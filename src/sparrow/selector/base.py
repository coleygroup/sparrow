from abc import ABC, abstractmethod
from typing import Dict, Union, List
from rdkit import Chem
from tqdm import tqdm
from pathlib import Path
from datetime import datetime
import warnings
import csv 
import json 
import numpy as np

from sparrow.scorer import Scorer
from sparrow.condition_recommender import Recommender
from sparrow.route_graph import RouteGraph
from sparrow.coster import Coster, LookupCoster
from sparrow.nodes import ReactionNode

from pulp import LpProblem
import gurobipy as gp 

Problem = Union[LpProblem, gp.Model]

class Selector(ABC):
    """ 
    A Selector performs the selection of molecules and their synthetic routes. 
    The selection is performed on a RouteGraph using PuLP to set up the 
    optimization problem and Gurobi to solve it. 
    This class also collects all information required for the RouteGraph, 
    such as reaction scoring, condition recommendation, and assigning compound
    costs. 
    """
    def __init__(self, 
                 route_graph: RouteGraph, 
                 target_dict: Dict[str, float],                 
                 rxn_scorer: Scorer = None, 
                 condition_recommender: Recommender = None,
                 constrain_all_targets: bool = False, 
                 max_targets: int = None,
                 coster: Coster = None, 
                 cost_per_rxn: float = 100,
                 weights: List = [1,1,1,0,0],
                 output_dir: str = 'debug',
                 remove_dummy_rxns_first: bool = False,
                 clusters: dict = None,
                 N_per_cluster: int = 1, 
                 rxn_classes: dict = None,
                 rxn_classifier_dir: str = None,
                 max_rxn_classes: int = None,
                 max_rxns: int = None,
                 sm_budget: float = None,
                 dont_buy_targets: bool = False, 
                 cycle_constraints: bool = False,
                 max_seconds: int = None,
                 extract_routes: bool = True, 
                 ) -> None:

        self.dir = Path(output_dir)
        self.cost_per_rxn = cost_per_rxn
        self.max_rxns = max_rxns
        self.sm_budget = sm_budget
        self.runtime = None
        self.cycle_constraints = cycle_constraints
        self.max_seconds = max_seconds
        self.extract_routes = extract_routes

        self.graph = route_graph  
        if remove_dummy_rxns_first: 
            self.graph.remove_dummy_rxns()
        else: 
            self.graph.prune_dummy_rxns()

        Path(self.dir/'chkpts').mkdir(parents=True, exist_ok=True)
        save_freq = 1e6 if type(coster) == LookupCoster else 500
        self.graph.set_buyable_compounds_and_costs(coster, save_json_dir=self.dir/'chkpts', save_freq=save_freq)

        self.add_dummy_starting_rxn_nodes()

        self.graph.id_nodes()

        self.clean_target_dict(target_dict)
        self.targets = list(self.target_dict.keys())

        if clusters: 
            self.clusters, self.id_to_clusters = self.clean_clusters(clusters) 
        else: 
            self.clusters = self.id_to_clusters = None
        self.N_per_cluster = N_per_cluster

        if rxn_classes:
            self.rxn_classes = rxn_classes.keys()
            self.rxn_class_dict, self.id_to_classes = self.clean_rxn_classes(rxn_classes)
        else:
            self.rxn_classes = self.rxn_class_dict = None
        self.rxn_classifier_dir = rxn_classifier_dir

        self.target_dict = self.graph.set_compound_types(self.target_dict, coster=coster, save_dir=self.dir/'chkpts')

        if dont_buy_targets: 
            self.graph.remove_dummies_to(id_list=self.targets)
        
        self.rxn_scorer = rxn_scorer
        self.condition_recommender = condition_recommender
              
        self.constrain_all_targets = constrain_all_targets
        self.weights = weights
        self.max_targets = max_targets
        self.max_rxn_classes = max_rxn_classes
        
        if self.condition_recommender is not None: 
            self.get_recommendations()

        if self.rxn_scorer is not None: 
            self.get_rxn_scores()

        self.problem = self.initialize_problem()

    @abstractmethod
    def initialize_problem(self) -> Problem:
        """ Initialize a PuLP LpProblem or Gurobi Model """
    
    def add_dummy_starting_rxn_nodes(self): 
        """ Adds reaction nodes that form all starting materials """
        for start_node in self.graph.buyable_nodes(): 
            dummy_rxn_smiles = f">>{start_node.smiles}"
            self.graph.add_reaction_node(
                dummy_rxn_smiles, 
                children=[start_node.smiles], 
                dummy=True, 
                penalty=0, 
                score = 1,
            )  

    def clean_target_dict(self, target_dict: Dict[str, float]) -> Dict[str, float]:
        """ Converts target dict from Dict[smiles, reward] to Dict[id, reward] """
        self.target_dict = {}
        
        c=0
        old_smis = []
        for old_smi, reward in tqdm(target_dict.items(), desc='Checking targets and rewards'):
             
            try: 
                float(reward)
            except: 
                # warnings.warn(f'Target {old_smi} has an invalid reward ({reward}) and is being removed from target set')
                c+=1      
                continue     
            
            node = self.graph.compound_nodes.get(old_smi, None)
            if node:
                self.target_dict[node.id] = reward
                old_smis.append(old_smi)
            else: 
                clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
                node = self.graph.compound_nodes.get(clean_smi, None)
                
                if node: 
                    self.target_dict[node.id] = reward    
                    old_smis.append(old_smi)       
                else:       
                    # warnings.warn(f'Target {old_smi} is not in routes and is being removed from target set')
                    c+=1

        if len(self.target_dict) < len(target_dict):
            p = self.dir / 'cleaned_tar_dict.csv'

            warnings.warn(f'{c} targets are not in routes and are being removed from the candidate set')
            print(f'Saving {len(self.target_dict)} remaining targets, ids, and rewards to {p}')

            save_list = [
                {'SMILES': self.graph.smiles_from_id(id), 'Old SMILES': old_smi, 'ID': id, 'Reward': reward}
                for (id, reward), old_smi in zip(self.target_dict.items(), old_smis)
            ]
            
            with open(p, 'w') as csvfile: 
                writer = csv.DictWriter(csvfile, fieldnames=['SMILES', 'Old SMILES', 'ID', 'Reward'])
                writer.writeheader()
                writer.writerows(save_list)

        return self.target_dict
    
    def clean_rxn_classes(self, classes): 
        clean_clusters = {}
        for c_name, rxnsmis in classes.items(): 
            clean_smis = []
            for old_smi in rxnsmis: 
                if old_smi in self.graph.reaction_nodes:
                    clean_smis.append(self.graph.reaction_nodes[old_smi].id)
                else: 
                    reacs, prods = old_smi.split('>>')
                    clean_smi = f'{Chem.MolToSmiles(Chem.MolFromSmiles(reacs))}>>{Chem.MolToSmiles(Chem.MolFromSmiles(prods))}'
                    if clean_smi in self.graph.reaction_nodes: 
                        clean_smis.append(self.graph.reaction_nodes[clean_smi].id)       
            
            clean_clusters[c_name] = clean_smis  

        id_to_classes = {}
        for class_name, rxn_list in clean_clusters.items():
            for rxn in rxn_list:
                if rxn not in id_to_classes.keys():
                    id_to_classes[rxn] = [class_name]
                else:
                    id_to_classes[rxn].append(class_name)
        return clean_clusters, id_to_classes
        
    def clean_clusters(self, clusters): 
        clean_clusters = {}
        for c_name, smis in clusters.items(): 
            clean_smis = []
            for old_smi in smis: 
                if old_smi in self.graph.compound_nodes:
                    clean_smis.append(self.graph.compound_nodes[old_smi].id)
                else: 
                    clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
                    if clean_smi in self.graph.compound_nodes: 
                        clean_smis.append(self.graph.compound_nodes[clean_smi].id)       
            
            clean_clusters[c_name] = clean_smis     

        id_to_clusters = {tid: [] for tid in self.targets}
        for c_name, tids in clean_clusters.items():
            for tid in tids: 
                id_to_clusters[tid].append(c_name)

        return clean_clusters, id_to_clusters

    def get_recommendations(self): 
        """ Completes condition recommendation for any reaction node that does not have conditions """
        count = 0
        for node in tqdm(self.graph.non_dummy_nodes(), 'Recommending Conditions'):
            if node.condition_set: 
                continue
            
            condition = self.condition_recommender(node.smiles)
            node.update_condition(condition)
            count += 1

            if count % 1000 == 0: 
                time = datetime.now().strftime("%H-%M-%S")
                self.graph.to_json(self.dir / 'chkpts' / f'trees_w_conditions_{time}.json')
            
        self.graph.to_json(self.dir / 'chkpts' / f'trees_w_conditions.json')

    def get_rxn_scores(self): 
        """ Scores all reactions in the graph that are not already scored """
        count = 0
        for node in tqdm(self.graph.reaction_nodes_only(), 'Scoring reactions'): 
            if (node.score_set and node.score > 0) or node.dummy: 
                continue 

            try:    
                score = self.rxn_scorer(rxn_smi=node.smiles, condition=node.condition)
            except: 
                try: 
                    print(f'Could not score {node.smiles} with conditions {node.condition}, trying again without conditions')
                    score = self.rxn_scorer(rxn_smi=node.smiles, condition=[[]])
                except:
                    print(f'Could not score {node.smiles}, returning score of 0')
                    score = 0 

            node.update(score=score)
            count += 1
            if count % 100 == 0: 
                time = datetime.now().strftime("%H-%M-%S")
                self.graph.to_json(self.dir / 'chkpts' / f'trees_w_scores_{time}.json')
        
        self.graph.to_json(self.dir / 'chkpts' / f'trees_w_scores.json')
    
    def cost_of_dummy(self, dummy_node: ReactionNode = None, dummy_id: str = None) -> float:
        if dummy_node is None: 
            dummy_node = self.graph.node_from_id(dummy_id)

        start_node = list(dummy_node.children.values())[0]
        return start_node.cost_per_g 
    
    @abstractmethod
    def define_variables(self): 
        """ Defines decision variables for self.problem """
    
    @abstractmethod
    def set_constraints(self): 
        """ Sets constraints for self.problem """
        
    @abstractmethod
    def set_objective(self): 
        """ Sets objective function for self.problem """

    @abstractmethod
    def optimize(self): 
        """ Optimizes self.problem, with a time limit of self.max_seconds """

    def formulate_and_optimize(self, extract_routes: bool = True): 
        """ Formulates the problem, solves it, and outputs the solution """
        self.define_variables()
        self.set_objective()
        self.set_constraints()
        self.optimize()
        self.post_processing(extract_routes=extract_routes)
        return 

    def extract_vars(self, output_dir=None, extract_routes=True): 
        mol_ids, rxn_ids, class_ids = self.extract_selected_ids()

        dummy_ids = [rxn for rxn in rxn_ids if self.graph.node_from_id(rxn).dummy]
        non_dummy_ids = [rxn for rxn in rxn_ids if self.graph.node_from_id(rxn).dummy == 0]

        # DFS to check that selected nodes for rxns and compounds are connected
        selected_targets = set(mol_ids) & set(self.targets)
        selected_starting = set([self.graph.child_of_dummy(dummy) for dummy in dummy_ids])
        print(f'{len(selected_targets)} targets selected using {len(non_dummy_ids)} reactions and {len(selected_starting)} starting materials')
        self.export_selected_nodes(rxn_ids, selected_starting, selected_targets, output_dir)
        
        avg_rxn_score = np.mean([self.graph.node_from_id(rxn).score for rxn in non_dummy_ids]) if len(non_dummy_ids) > 0 else None

        summary = {
            'Number targets': len(selected_targets), 
            'Fraction targets selected': len(selected_targets)/len(self.targets),
            'Cumulative reward of selected compounds': sum([self.target_dict[tar] for tar in selected_targets]),
            'Possible reward': sum(self.target_dict.values()),
            'Number starting materials (some may be unused)': len(selected_starting),
            'Cost starting materials (some may be unused)': sum([self.cost_of_dummy(dummy_id=d_id) for d_id in dummy_ids]),
            'Number reaction steps (some may be unused)': len(non_dummy_ids),
            'Average reaction score (some may be unused)': avg_rxn_score,
            'Run time': self.runtime,
            'Number of variables': self.get_num_variables(),
            'Number of constraints': self.get_num_constraints(),
        }

        if self.rxn_classes: 
            summary['Number of reaction classes selected'] = len(class_ids)

        if extract_routes:
            summary['Expected Reward'] = 0
            storage = {}
            for target in tqdm(selected_targets, desc='Extracting routes'): 
                store_dict = {'Compounds':[], 'Reactions':[]}
                smi = self.graph.smiles_from_id(target)
                # recursive search between rxn and compound parent funcs
                storage[smi] = self.find_mol_parents(store_dict, target, mol_ids, rxn_ids, target=target)
                storage[smi]['Reward'] = self.target_dict[target]
                p_success = np.prod(np.array([entry['score'] for entry in storage[smi]['Reactions'] if 'score' in entry]))
                er = self.target_dict[target]*p_success
                storage[smi]['Expected Reward'] = er
                summary['Expected Reward'] += er

            with open(output_dir/f'routes.json','w') as f: 
                json.dump(storage, f, indent='\t')

            print(f'Total expected reward: {summary["Expected Reward"]:0.2f}')

            # calculate actually used starting materials and reactions 
            all_rxns = [[rentry['smiles'] for rentry in entry['Reactions']] for entry in storage.values()]
            all_rxns = set([rxn for entry in all_rxns for rxn in entry])
            rxn_steps = [rxn for rxn in all_rxns if not rxn.startswith('>>')]
            starting_materials = [rxn.split('>>')[1] for rxn in all_rxns if rxn.startswith('>>')]
            summary['Number reaction steps (actual)'] = len(rxn_steps)
            summary['Number starting materials (actual)'] = len(starting_materials)
            summary['Average reaction score (actual)'] =  np.mean([self.graph.node_from_smiles(rxn).score for rxn in rxn_steps]) if len(rxn_steps) > 0 else None

        return summary 
    
    def rxn_class_analysis(self, output_dir):

        with open(output_dir/f'routes.json', 'r') as f: 
            data = json.load(f)

        dict = {} # target to classes in order
        class_count = {}
        score_dict = {} # general shared reaction count
        
        # need rxn classes in order that were used to get to target
        for target in data:
            reactions = data[target]['Reactions']
            dict[target] = []
            for reaction in reactions:
                if 'class' in reaction: 
                    rxn_class = reaction['class'][0]
                    dict[target].append(rxn_class)
                    class_count[rxn_class] = class_count.get(rxn_class, 0) + 1

        # this just counts the shared reactions between the path to a target compound and the paths to all other targets
        for target in dict:
            score_dict[target] = 0
            for reaction in dict[target]:
                for other_target in dict:
                    if other_target != target:
                        if reaction in dict[other_target]:
                            score_dict[target] += 1

        results = {}
        results["Frequency of selected reaction classes"] = class_count
        return results

    def post_processing(self, 
                        output_dir: str = None, 
                        extract_routes: bool = True): 
        
        if output_dir is None: 
            output_dir = self.dir / 'solution'

        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True) 

        summary = self.extract_vars(output_dir=output_dir, extract_routes=extract_routes)

        if self.rxn_classes:
            rxn_class_summary = self.rxn_class_analysis(output_dir)
            with open(output_dir/'rxn_class_summary.json', 'w') as f:
                json.dump(rxn_class_summary, f, indent='\t')
                
        with open(output_dir/'summary.json', 'w') as f:
            json.dump(summary, f, indent='\t')
        
        return summary

    @abstractmethod
    def extract_selected_ids(self):
        """ Returns a set of node IDs corresponding to selected compounds and reactions """

    def export_selected_nodes(self, rxn_list, starting_list, target_list, output_dir):
        storage = {'Starting Materials': [], 'Reactions': [], 'Targets': []}
        graph = self.graph

        for rxn_id in rxn_list:
            node = graph.node_from_id(rxn_id)
            if node.dummy: 
                storage['Reactions'].append({
                    'smiles': node.smiles,
                    'starting material cost ($/g)': self.cost_of_dummy(node), 
                })
            else: 
                new_rxn_entry = {
                    'smiles': node.smiles,
                    'conditions': node.get_condition(1)[0], 
                    'score': node.score,
                }
                if self.rxn_classes:
                    new_rxn_entry['class'] = self.id_to_classes[rxn_id] 
                storage['Reactions'].append(new_rxn_entry)


        for cpd_id in starting_list: 
            node = graph.node_from_id(cpd_id)
            storage['Starting Materials'].append({
                'smiles': node.smiles,
                'cost': node.cost_per_g
            })

        for cpd_id in target_list: 
            node = graph.node_from_id(cpd_id)
            if self.clusters: 
                storage['Targets'].append({
                    'smiles': node.smiles,
                    'reward': node.reward,
                    'clusters': self.id_to_clusters[cpd_id]
                })
            else:
                storage['Targets'].append({
                    'smiles': node.smiles,
                    'reward': node.reward,
                })
        
        with open(Path(output_dir)/'solution_list_format.json','w') as f:
            json.dump(storage, f, indent='\t')

        return storage 

    def find_rxn_parents(self, store_dict, rxn_id, selected_mols, selected_rxns, target=None):

        par_ids = [n.id for n in self.graph.node_from_id(rxn_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_mols)
        for par in selected_pars: 
            store_dict['Compounds'].append(self.graph.smiles_from_id(par))
            store_dict = self.find_mol_parents(store_dict, par, selected_mols, selected_rxns)
        return store_dict

    def find_mol_parents(self, store_dict, mol_id, selected_mols, selected_rxns, target=None): 

        par_ids = [n.id for n in self.graph.node_from_id(mol_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_rxns)
        for par in selected_pars: 
            node = self.graph.node_from_id(par)
            if node.dummy: 
                store_dict['Reactions'].append({
                    'smiles': node.smiles,
                    'starting material cost ($/g)': self.cost_of_dummy(node), 
                })
            else: 
                new_rxn_entry = {
                    'smiles': node.smiles,
                    'conditions': node.get_condition(1)[0], 
                    'score': node.score,
                }
                if self.rxn_classes: 
                    new_rxn_entry['class'] = self.id_to_classes[par]
                store_dict['Reactions'].append(new_rxn_entry)  

            store_dict = self.find_rxn_parents(store_dict, par, selected_mols, selected_rxns)
        return store_dict

    @abstractmethod
    def get_num_variables(self):
        """ Returns the number of variables in the optimization problem """

    @abstractmethod
    def get_num_constraints(self): 
        """ Returns the number of constraints in the optimization problem """
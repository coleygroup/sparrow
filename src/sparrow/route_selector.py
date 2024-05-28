from sparrow.scorer import Scorer
from sparrow.condition_recommender import Recommender
from sparrow.route_graph import RouteGraph
from sparrow.coster import Coster, ChemSpaceCoster
from sparrow.nodes import ReactionNode
from sparrow.utils.cluster_utils import cluster_smiles
from typing import Dict, Union, List
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI
from pulp.apis import PULP_CBC_CMD
from rdkit import Chem
from tqdm import tqdm
from pathlib import Path
from datetime import datetime
import time
import warnings
import csv 
import json 

reward_type = Union[int, float]

class RouteSelector: 
    """ 
    RouteSelector performs the selection of molecules and their synthetic routes. 
    The selection is performed on a RouteGraph using PuLP to set up the 
    optimization problem and Gurobi to solve it. 
    """
    def __init__(self, 
                 route_graph: RouteGraph, 
                 target_dict: Dict[str, reward_type],                 
                 rxn_scorer: Scorer = None, 
                 condition_recommender: Recommender = None,
                 constrain_all_targets: bool = False, 
                 max_targets: int = None,
                 coster: Coster = None, 
                 weights: List = [1,1,1,1],
                 output_dir: str = 'debug',
                 remove_dummy_rxns_first: bool = False,
                 cluster_cutoff: float = 0.7,
                 dont_buy_targets: bool = False,
                 ) -> None:

        self.dir = Path(output_dir)

        self.graph = route_graph  
        if remove_dummy_rxns_first: 
            self.graph.remove_dummy_rxns()
        else: 
            self.graph.prune_dummy_rxns()

        Path(self.dir/'chkpts').mkdir(parents=True, exist_ok=True)
        if type(coster) == ChemSpaceCoster: 
            self.graph.set_buyable_compounds_and_costs(coster, save_json_dir=self.dir/'chkpts')
        else: 
            self.graph.set_buyable_compounds_and_costs(coster, save_json_dir=self.dir/'chkpts', save_freq=1e6)

        self.add_dummy_starting_rxn_nodes()

        self.graph.id_nodes()

        self.clean_target_dict(target_dict)
        self.targets = list(self.target_dict.keys())

        self.target_dict = self.graph.set_compound_types(self.target_dict, coster=coster, save_dir=self.dir/'chkpts')

        if dont_buy_targets: 
            self.graph.remove_dummies_to(id_list=self.targets)

        self.rxn_scorer = rxn_scorer
        self.condition_recommender = condition_recommender
              
        self.constrain_all_targets = constrain_all_targets
        self.weights = weights
        self.cluster_cutoff = cluster_cutoff
        self.max_targets = max_targets
        
        if self.condition_recommender is not None: 
            self.get_recommendations()

        if self.rxn_scorer is not None: 
            self.get_rxn_scores()

        self.problem = LpProblem("Route_Selection", LpMinimize)

    def clean_target_dict(self, target_dict: Dict[str, float]) -> Dict[str, float]:
        """ Converts target dict from Dict[smiles, reward] to Dict[id, reward] """
        self.target_dict = {}
        
        c=0
        for old_smi, reward in tqdm(target_dict.items(), desc='Checking targets and rewards'):
             
            try: 
                float(reward)
            except: 
                warnings.warn(f'Target {old_smi} has an invalid reward ({reward}) and is being removed from target set')
                c+=1      
                continue     
            
            node = self.graph.compound_nodes.get(old_smi, None)
            if node:
                self.target_dict[node.id] = reward
            else: 
                clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
                node = self.graph.compound_nodes.get(clean_smi, None)
                if node: 
                    self.target_dict[node.id] = reward          
                else:       
                    warnings.warn(f'Target {old_smi} is not in routes and is being removed from target set')
                    c+=1

        if len(self.target_dict) < len(target_dict):
            p = self.dir / 'cleaned_tar_dict.csv'
            print(f'Saving {len(self.target_dict)} remaining targets, ids, and rewards to {p}')

            save_list = [
                {'SMILES': self.graph.smiles_from_id(id), 'ID': id, 'Reward': reward,}
                for id, reward in self.target_dict.items()
            ]
            
            with open(p, 'w') as csvfile: 
                writer = csv.DictWriter(csvfile, fieldnames=['SMILES', 'ID', 'Reward'])
                writer.writeheader()
                writer.writerows(save_list)

        return self.target_dict

    def get_recommendations(self): 
        """ Completes condition recommendation for any reaction node that does not have conditions """
        count = 0
        for node in tqdm(self.graph.non_dummy_nodes(), 'Recommending Conditions'):
            if node.condition_set: 
                continue
            
            condition = self.condition_recommender(node.smiles)
            node.update_condition(condition)
            count += 1

            if count % 100 == 0: 
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

    def define_variables(self): 
        """ 
        TODO: explain in readme what variables mean, refer to that here 
        (currently in my thesis proposal)
        TODO: include conditions 
        """

        rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]
        self.r = LpVariable.dicts(
            "rxn", 
            indices=rxn_ids, 
            cat="Binary",
        )

        mol_ids = [node.id for node in self.graph.compound_nodes_only()]
        self.m = LpVariable.dicts(
            "mol", 
            mol_ids, 
            cat="Binary",
        )
        
        return 

    def set_constraints(self, set_cycle_constraints=True):
        """ Sets constraints defined in TODO: write in README all constraints """
        print('Setting constraints ...')
        # implement constrain_all_targets later

        self.set_rxn_constraints()
        self.set_mol_constraints()

        if set_cycle_constraints: 
            self.set_cycle_constraints()
        
        if self.max_targets:
            self.set_max_target_constraint()
        
        if self.constrain_all_targets:
            self.set_constraint_all_targets()

        return 
    
    def set_rxn_constraints(self): 

        for node in tqdm(self.graph.reaction_nodes_only(), desc='Reaction constraints'): 
            if node.dummy: 
                continue 
            par_ids = [par.id for par in node.parents.values()]
            for par_id in par_ids: 
                self.problem += (
                    self.m[par_id] >= self.r[node.id]
                )
        
        return 
    
    def set_mol_constraints(self): 

        for node in tqdm(self.graph.compound_nodes_only(), 'Compound constraints'): 
            parent_ids = [par.id for par in node.parents.values()]
            self.problem += (
                self.m[node.id] <= lpSum(self.r[par_id] for par_id in parent_ids)
            )
        
        return 

    def set_cycle_constraints(self): 

        cycles = self.graph.dfs_find_cycles_nx()
        for cyc in tqdm(cycles, desc='Cycle constraints'): 
            self.problem += (
                lpSum(self.r[rid] for rid in cyc) <= (len(cyc) - 1)
            )

        return 
    
    def set_max_target_constraint(self):
        """ Sets constraint on the maximum number of selected targets. """
        self.problem += (
            lpSum(self.m[target] for target in self.targets) <= self.max_targets
        )
    
    def set_constraint_all_targets(self): 
        """ Constrains that all targets must be synthesized """
        for target in self.targets: 
            self.problem += (
                self.m[target] == 1
            )
    
    def get_child_and_parent_ids(self, smi: str = None, id: str = None): 
        """ Returns list of child node smiles and parent node smiles for a given
        compound smiles """
        if smi is None and id is None: 
            print('No node information given')
            return None 
        
        if id is not None: 
            child_ids = [child.id for child in self.graph.node_from_id(id).children.values()] 
            parent_ids = [parent.id for parent in self.graph.node_from_id(id).parents.values()]
        elif smi is not None: 
            child_ids = [child.id for child in self.graph.node_from_smiles(smi).children.values()] 
            parent_ids = [parent.id for parent in self.graph.node_from_smiles(smi).parents.values()]

        return parent_ids, child_ids
    
    def add_dummy_starting_rxn_nodes(self): 
        """ Adds reaction nodes that form all starting materials, as described in 
        TODO: describe this in README """
        for start_node in self.graph.buyable_nodes(): 
            dummy_rxn_smiles = f">>{start_node.smiles}"
            self.graph.add_reaction_node(
                dummy_rxn_smiles, 
                children=[start_node.smiles], 
                dummy=True, 
                penalty=0, 
                score = 10**6,
            )
    
    def set_objective(self): 
        # TODO: Add consideration of conditions 
        print('Setting objective function ...')

        reward_mult = self.weights[0] # / ( len(self.target_dict)) # *max(self.target_dict.values()) )
        cost_mult = self.weights[1] # / (len(self.graph.dummy_nodes_only())) # * max([node.cost_per_g for node in self.graph.buyable_nodes()]) ) 
        pen_mult = self.weights[2] # / (len(self.graph.non_dummy_nodes())) # * max([node.penalty for node in self.graph.non_dummy_nodes()]) )

        self.problem += -1*reward_mult*lpSum([float(self.target_dict[target])*self.m[target] for target in self.targets]) \
        + cost_mult*lpSum([self.cost_of_dummy(dummy)*self.r[dummy.id] for dummy in self.graph.dummy_nodes_only()]) \
        + pen_mult*lpSum([self.r[node.id]*float(node.penalty) for node in self.graph.non_dummy_nodes()])
            # reaction penalties, implement CSR later 
        
        if self.weights[3]>0: 
            self.add_diversity_objective()

        return 

    def add_diversity_objective(self): 
        """ Adds scalarization objective to increase the number of clusters represented, requires defining new variable """
        print('Clustering molecules for diversity objective')
        cs_ind = cluster_smiles([self.graph.smiles_from_id(id) for id in self.targets], cutoff=self.cluster_cutoff)
        cs = [[self.targets[ind] for ind in cluster] for cluster in cs_ind]
        
        cs_file = self.dir / 'clusters.json'
        print(f'Saving list of {len(cs)} clusters to {cs_file}')
        
        with open(cs_file,'w') as f: 
            json.dump(cs, f, indent='\t')

        # d_i : whether cluster i is represented by the selected set 
        self.d = LpVariable.dicts(
            "cluster", 
            indices=range(len(cs)), 
            cat="Binary",
        )

        # constraint: d_i <= sum(c_j) for j in cluster i
        for i, ids_in_cluster in enumerate(cs): 
            self.problem += (
                self.d[i] <= lpSum(self.m[cpd_id] for cpd_id in ids_in_cluster)
            )
        
        # add objective 
        print(f'adding objective with {self.weights[3]}')
        self.problem += self.problem.objective - self.weights[3]*lpSum(self.d)

        return 

    
    def cost_of_dummy(self, dummy_node: ReactionNode = None, dummy_id: str = None) -> float:
        if dummy_node is None: 
            dummy_node = self.graph.node_from_id(dummy_id)

        start_node = list(dummy_node.children.values())[0]
        return start_node.cost_per_g 

    def optimize(self, solver=None):

        # self.problem.writeLP("RouteSelector.lp", max_length=300)
        print("Solving optimization problem...")
        opt_start = time.time()
        if solver == 'GUROBI': 
            self.problem.solve(GUROBI(timeLimit=86400))
        else: 
            self.problem.solve(PULP_CBC_CMD(gapRel=1e-7, gapAbs=1e-9, msg=False))

        print(f"Optimization problem completed. Took {time.time()-opt_start:0.2f} seconds to solve")
        
        return 
    
    def optimal_variables(self):
        """ Returns nonzero variables """
        nonzero_vars = [
            var for var in self.problem.variables() if var.varValue > 0.01
        ]

        return nonzero_vars

     
    
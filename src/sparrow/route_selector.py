from sparrow.scorer import Scorer
from sparrow.condition_recommender import Recommender
from sparrow.route_graph import RouteGraph
from sparrow.coster import Coster, ChemSpaceCoster
from sparrow.nodes import ReactionNode
from sparrow.utils.cluster_utils import cluster_smiles
from typing import Dict, Union, List
from rdkit import Chem
from tqdm import tqdm
from pathlib import Path
from datetime import datetime
from math import log 
import time
import warnings
import csv 
import re
import json 
import numpy as np 

import gurobipy as gp 
from gurobipy import GRB

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
                 cost_per_rxn: float = 100,
                 output_dir: str = 'debug',
                 remove_dummy_rxns_first: bool = False,
                 cluster_cutoff: float = 0.7,
                 custom_clusters: dict = None,
                 max_rxns: int = None,
                 sm_budget: float = None,
                 ) -> None:

        self.dir = Path(output_dir)
        self.cost_per_rxn = cost_per_rxn
        self.max_rxns = max_rxns
        self.sm_budget = sm_budget

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

        if custom_clusters: 
            self.clusters = self.clean_clusters(custom_clusters) 
        else: 
            self.clusters = None

        self.target_dict = self.graph.set_compound_types(self.target_dict, coster=coster, save_dir=self.dir/'chkpts')

        self.rxn_scorer = rxn_scorer
        self.condition_recommender = condition_recommender
              
        self.constrain_all_targets = constrain_all_targets
        self.cluster_cutoff = cluster_cutoff
        self.max_targets = max_targets
        
        if self.condition_recommender is not None: 
            self.get_recommendations()

        if self.rxn_scorer is not None: 
            self.get_rxn_scores()

        self.model = gp.Model("sparrow")

    def clean_target_dict(self, target_dict: Dict[str, float]) -> Dict[str, float]:
        """ Converts target dict from Dict[smiles, reward] to Dict[id, reward] """
        self.target_dict = {}
        
        c=0
        old_smis = []
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
                old_smis.append(old_smi)
            else: 
                clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
                node = self.graph.compound_nodes.get(clean_smi, None)
                
                if node: 
                    self.target_dict[node.id] = reward
                    old_smis.append(old_smi)          
                else:       
                    warnings.warn(f'Target {old_smi} is not in routes and is being removed from target set')
                    c+=1

        if len(self.target_dict) < len(target_dict):
            p = self.dir / 'cleaned_tar_dict.csv'
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

    def clean_clusters(self, clusters): 

        clean_clusters = {}
        for c_name, smis in clusters.items(): 
            clean_smis = []
            for old_smi in smis: 
                node = self.graph.compound_nodes.get(old_smi, None)
                if node:
                    clean_smis.append(node.id)
                else: 
                    clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
                    node = self.graph.compound_nodes.get(clean_smi, None)
                    if node: 
                        clean_smis.append(node.id)       
                    else:       
                        warnings.warn(f'Target {old_smi} is not in routes and is being removed from clusters')
            
            clean_clusters[c_name] = clean_smis     

        return clean_clusters

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

        # whether rxn is used 
        self.rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]
        self.r = self.model.addVars(
            self.rxn_ids,
            vtype=GRB.BINARY,
            name="rxn",
        )

        # whether molecule is selected 
        mol_ids = [node.id for node in self.graph.compound_nodes_only()]
        self.m = self.model.addVars(
            mol_ids,
            vtype=GRB.BINARY,
            name="mol"
        )

        # whether rxn i is used for target j 
        self.u = self.model.addVars(
            self.rxn_ids, self.targets, 
            vtype=GRB.BINARY,
            name="rxnfor"
        )

        # whether compound j is used int he route to target j 
        self.cpdfor = self.model.addVars(
            mol_ids, self.targets,
            vtype=GRB.BINARY,
            name="cpdfor"
        )

        # additional variables for nonlinearity 
        sum_log_scores = np.log([node.score if node.score>0 else 1e-10 for node in self.graph.non_dummy_nodes()]).sum()
        self.usum = self.model.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,  
            name="usum",
            ub=0,
            lb=sum_log_scores,
        )
        self.expusum = self.model.addVars(
            self.targets,
            vtype=GRB.CONTINUOUS,
            name="expusum",
            ub=1,
            lb=0
        )
        self.probsuccess = self.model.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,
            name="probsuccess",
            ub=1,
            lb=0,
        )
        max_n_rxns = len(self.graph.non_dummy_nodes())
        max_sm_costs = np.sum([self.cost_of_dummy(dummy) for dummy in self.graph.dummy_nodes_only()])
        max_costs = max_n_rxns*self.cost_per_rxn + max_sm_costs
        self.allcosts = self.model.addVar(
            name="allcosts",
            vtype=GRB.CONTINUOUS,
            ub=max_costs,
            lb=0,
        )
        max_er = sum(self.target_dict.values())
        self.sumer = self.model.addVar(
            name="sumer",
            vtype=GRB.CONTINUOUS,
            ub=max_er,
            lb=0,
        )
        # self.logsumer = self.model.addVar(
        #     name="logsumer",
        #     vtype=GRB.CONTINUOUS,
        #     ub=1e20,
        #     lb=0,
        # )
        # self.logallcosts = self.model.addVar(
        #     name='logallcosts',
        #     vtype=GRB.CONTINUOUS,
        #     ub=1e20,
        #     lb=0,
        # )
        
        return 

    def set_constraints(self, set_cycle_constraints=True):
        """ Sets constraints """
        print('Setting constraints ...')
        # implement constrain_all_targets later

        self.set_rxn_constraints()
        self.set_mol_constraints()

        if set_cycle_constraints: 
            self.set_cycle_constraints()
        
        if self.max_targets:
            self.set_max_target_constraint()
        
        if self.constrain_all_targets:
            self.set_constraint_n_targets(N=len(self.targets))

        if self.clusters: 
            self.set_cluster_constraints()

        if self.sm_budget: 
            self.set_budget_constraint()

        self.set_nonlinear_constraints()

        return 
    
    def set_nonlinear_constraints(self): 
        # sets exponential constraints for usum 
        # usum = sum_i u_i,j * log (Li)
        log_scores = np.log([node.score if node.score>0 else 1e-10 for node in self.graph.non_dummy_nodes()])
        # log_scores = [log(node.score) if node.score>0 else -1e20 for node in self.graph.non_dummy_nodes()]
        for t in self.targets: 
            self.model.addConstr(
                self.usum[t] == gp.quicksum(
                    self.u[node.id, t]*logscore 
                    for logscore, node in zip(log_scores, self.graph.non_dummy_nodes())
                ),
                name=f'usumConstr_{t}',
            )
            self.model.addGenConstrExp(self.usum[t], self.expusum[t], name=f'expusumConstr_{t}',options="FuncPieces=-2 FuncPieceError=0.00001")
            self.model.addQConstr(
                self.probsuccess[t], GRB.LESS_EQUAL, self.expusum[t]*self.cpdfor[t, t], 
                name=f'probsuccessConstr_{t}',
            )
        
        # self.model.addConstr(
        #     self.allcosts == gp.quicksum(
        #         [self.cost_per_rxn*gp.quicksum(self.r[node.id]*float(node.penalty) for node in self.graph.non_dummy_nodes()),
        #         gp.quicksum(self.cost_of_dummy(dummy)*self.r[dummy.id] for dummy in self.graph.dummy_nodes_only()),]
        #     ),
        #     name='allcosts_constr'
        # )
        
        self.model.addQConstr(
            self.sumer, GRB.LESS_EQUAL, gp.quicksum(
                self.target_dict[t]*self.probsuccess[t] for t in self.targets
            ),
            name='sumER_constr',
        )
        # self.model.addGenConstrLog(
        #     self.allcosts,
        #     self.logallcosts,
        #     name=''
        # )

        # self.model.addGenConstrLog(
        #     self.sumer,
        #     self.logsumer,            
        # )

        return 
    
    def set_rxn_constraints(self): 
        # TODO: consider redudancy in these constraints and simplify problem 
        
        if self.max_rxns: 
            self.model.addConstr(
                self.max_rxns >= gp.quicksum(self.r[node.id] for node in self.graph.non_dummy_nodes())
            )

        for node in tqdm(self.graph.reaction_nodes_only(), desc='Reaction constraints'): 
            if node.dummy: 
                continue 
            par_ids = [par.id for par in node.parents.values()]
            for par_id in par_ids:
                # whether reaction selected at all
                self.model.addConstr(
                    self.m[par_id] >= self.r[node.id],
                    name=f'rxnConstr_{node.id}_{par_id}'
                )

                # # reaction for specific targets
                # self.model.addConstrs(
                #     (
                #         self.m[par_id] >= self.u[node.id, target]
                #         for target in self.targets 
                #     ),
                # )
                
                # if a reaction is selected, every parent (reactant) must also be selected 
                self.model.addConstrs(
                    self.cpdfor[par_id, t] >= self.u[node.id, t]
                    for t in self.targets
                )

            # if a reaction is used for a target, it is also used in general 
            self.model.addConstr(
                len(self.targets)*self.r[node.id] >= gp.quicksum(self.u[node.id, tar] for tar in self.targets),
                name=f'rxnusedConstr_{node.id}'
            )

            self.model.addConstr(
                self.r[node.id] <= gp.quicksum(self.u[node.id, tar] for tar in self.targets),
                name=f'rxnusedConstr2_{node.id}'
            )

        return 
    
    def set_mol_constraints(self): 
        
        for t in self.targets: 
            self.model.addConstr(
                self.cpdfor[t, t] == self.m[t]
            )

        for node in tqdm(self.graph.compound_nodes_only(), 'Compound constraints'): 
            par_ids = [par.id for par in node.parents.values()]
            self.model.addConstr(
                self.m[node.id] <= gp.quicksum(self.r[par_id] for par_id in par_ids),
                name=f'cpdConstr_{node.id}'
            )

            # if a compound is used for a target, at least one parent reaction must also be used for that target 
            self.model.addConstrs(
                self.cpdfor[node.id, t] <= gp.quicksum(self.u[par_id, t] for par_id in par_ids) 
                for t in self.targets
            )

            # if a compound is used for a target, it is also used in general 
            self.model.addConstr(
                len(self.targets)*self.m[node.id] >= gp.quicksum(self.cpdfor[node.id, tar] for tar in self.targets),
                name=f'molusedConstr_{node.id}'
            )     

            self.model.addConstr(
                self.m[node.id] <= gp.quicksum(self.cpdfor[node.id, tar] for tar in self.targets),
                name=f'molusedConstr_{node.id}'
            )      
        
        return 

    def set_cycle_constraints(self): 

        cycles = self.graph.dfs_find_cycles_nx()
        c = 0
        for cyc in tqdm(cycles, desc='Cycle constraints'): 
            self.model.addConstr(
                gp.quicksum(self.r[rid] for rid in cyc) <= (len(cyc) - 1),
                name=f'cycConstr_{c}'
            )
            c += 1

        return 
    
    def set_max_target_constraint(self):
        """ Sets constraint on the maximum number of selected targets. """
        self.model.addConstr(
            gp.quicksum(self.m[target] for target in self.targets) <= self.max_targets,
            name='max_targets'
        )
    
    def set_constraint_n_targets(self, N): 
        """ Constrains that n targets must be synthesized """
        
        self.model.addConstr(
            gp.quicksum(self.m[target] for target in self.targets) == N,
            name='constrain_n'
        )

        return 'constrain_n'
    
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
        """ Adds reaction nodes that form all starting materials """
        for start_node in self.graph.buyable_nodes(): 
            dummy_rxn_smiles = f">>{start_node.smiles}"
            self.graph.add_reaction_node(
                dummy_rxn_smiles, 
                children=[start_node.smiles], 
                dummy=True, 
                penalty=0, 
                score = 10**6,
            )
    
    def set_budget_constraint(self): 
        self.model.addConstr(
            gp.quicksum(self.cost_of_dummy(dummy)*self.r[dummy.id] for dummy in self.graph.dummy_nodes_only()) <= self.sm_budget
        )

    def set_objective(self, cost_weight=0): 
        # TODO: Add consideration of conditions 
        print('Setting objective function ...')

        self.model.ModelSense = GRB.MAXIMIZE

        if cost_weight > 0:
            self.model.setObjective(self.sumer - cost_weight*self.allcosts)
        else: 
            self.model.setObjective(self.sumer)

        return 
    
    def set_cluster_constraints(self): 
        """ 
        Constrains that every cluster must be represented by the selected set of candidates 
        clusters: a list of lists
                  each list represents a cluster, and each element of each list is a CompoundNode ID 
                  corresponding to the target in that cluster 
        """
        
        # self.crep = self.model.addVars(
        #     list(self.clusters.keys()),
        #     vtype=GRB.BINARY,
        #     name="cluster"
        # )

        for name, ids in self.clusters.items(): 
            self.model.addConstr(
                gp.quicksum(self.m[nid] for nid in ids) >= 1,
                name='cluster_represented'
            )

            # self.model.addConstr(
            #     self.crep[name] == 1
            # )

    def set_total_cost_constraints(self, cost_min, cost_max):

        self.model.addConstr(self.allcosts >= cost_min, name='min_cost')
        self.model.addConstr(self.allcosts <= cost_max, name='max_cost')

        return 'min_cost', 'max_cost'

    def add_diversity_objective(self): 
        """ Adds scalarization objective to increase the number of clusters represented, requires defining new variable """
        # NOT SUPPORTED WITH GUROBI 

        print("not supported with gurobi yet! ")

        return 
    
    def cost_of_dummy(self, dummy_node: ReactionNode = None, dummy_id: str = None) -> float:
        if dummy_node is None: 
            dummy_node = self.graph.node_from_id(dummy_id)

        start_node = list(dummy_node.children.values())[0]
        return start_node.cost_per_g 

    def optimize(self, savedir=None, solver=None):

        # self.problem.writeLP("RouteSelector.lp", max_length=300)
        if not savedir: 
            savedir = self.dir 

        print("Solving optimization problem...")
        opt_start = time.time()

        self.model.params.NonConvex = 2
        self.model.Params.TIME_LIMIT = 3*3600
        self.model.optimize()
        print(f"Optimization problem completed. Took {time.time()-opt_start:0.2f} seconds.")
        
        return 
    
    def optimize_MO(self, solver=None, set_cycle_constraints=True):

        self.model.params.outputflag = 0

        # first solve constraining all targets
        solution_dir = self.dir / 'solution_constrain_all'
        solution_dir.mkdir(exist_ok=True, parents=True)
        status = 'failed'
        N = len(self.targets)

        while status == 'failed': 
            constr_n = self.set_constraint_n_targets(N)
            self.model.write(str(solution_dir/"multiobj.lp"))
            self.set_objective(cost_weight=0.1)
            self.optimize()
            if self.model.SolCount > 0: 
                status = 'solved'
                cost_max = self.model.getVarByName('allcosts').X
                print(f'solved, max reasonable cost is {self.model.getVarByName("allcosts").X:0.2f}')
                self.extract_vars(solution_dir)
            else:  
                print(f'Could not solve constraining to {N} targets, reducing by 10%')
                N = min(int(np.ceil(0.9*N)), N-1)
            
            self.model.remove(self.model.getConstrByName(constr_n))


        # then solve constraining cost 
        self.set_objective()
        ranges = [[0.3, 0.5], [0.15, 0.3], [0.05, 0.15], [0, 0.05]]
        names = ['very high', 'high', 'medium', 'low']
        for name, range in zip(names, ranges): 
            solution_dir = self.dir / f'solution_{name}'
            solution_dir.mkdir(exist_ok=True, parents=True)
            min_constr, max_constr = self.set_total_cost_constraints(
                cost_min=range[0]*cost_max, 
                cost_max=range[1]*cost_max
            )
            self.model.write(str(solution_dir/"multiobj.lp"))
            self.optimize()
            print(f'solved, cost is {self.model.getVarByName("allcosts").X:0.2f}')
            self.extract_vars(solution_dir)
            self.model.remove(self.model.getConstrByName(min_constr))
            self.model.remove(self.model.getConstrByName(max_constr))
    
    def extract_vars(self, output_dir=None, extract_routes=True):
        if output_dir is None: 
            output_dir = self.dir / 'solution'

        output_dir.mkdir(exist_ok=True, parents=True)    
        nonzero_varnames = [
                var.VarName for var in self.model.getVars() if var.X > 0.01
            ]
        rxn_ids = re.findall(r'rxn\[(.*?)\]', ' '.join(nonzero_varnames))
        mol_ids = re.findall(r'mol\[(.*?)\]', ' '.join(nonzero_varnames))

        dummy_ids = [rxn for rxn in rxn_ids if self.graph.node_from_id(rxn).dummy]
        non_dummy_ids = [rxn for rxn in rxn_ids if self.graph.node_from_id(rxn).dummy == 0]

        selected_targets = set(mol_ids) & set(self.targets)
        selected_starting = set([self.graph.child_of_dummy(dummy) for dummy in dummy_ids])
        print(f'{len(selected_targets)} targets selected using {len(non_dummy_ids)} reactions and {len(selected_starting)} starting materials')
        self.export_selected_nodes(rxn_ids, selected_starting, selected_targets, output_dir)
        
        avg_rxn_score = np.mean([self.graph.node_from_id(rxn).score for rxn in non_dummy_ids]) if len(non_dummy_ids) > 0 else None

        summary = {
            'Cost/reaction weighting factor': self.cost_per_rxn,
            'Number targets': len(selected_targets), 
            'Fraction targets': len(selected_targets)/len(self.targets),
            'Total reward': sum([self.target_dict[tar] for tar in selected_targets]),
            'Possible reward': sum(self.target_dict.values()),
            'Number starting materials': len(selected_starting),
            'Cost starting materials': sum([self.cost_of_dummy(dummy_id=d_id) for d_id in dummy_ids]),
            'Number reaction steps': len(non_dummy_ids),
            'Average reaction score': avg_rxn_score,
        }

        if extract_routes:
            storage = {}
            for target in tqdm(selected_targets, desc='Extracting routes'): 
                store_dict = {'Compounds':[], 'Reactions':[]}
                smi = self.graph.smiles_from_id(target)
                storage[smi] = self.find_mol_parents(store_dict, target, mol_ids, rxn_ids)
                storage[smi]['Reward'] = self.target_dict[target]

            with open(output_dir/f'routes.json','w') as f: 
                json.dump(storage, f, indent='\t')

        with open(output_dir/'summary.json', 'w') as f:
            json.dump(summary, f, indent='\t')

        return summary 

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
                storage['Reactions'].append({
                    'smiles': node.smiles,
                    'conditions': node.get_condition(1)[0], 
                    'score': node.score,
                })

        for cpd_id in starting_list: 
            node = graph.node_from_id(cpd_id)
            storage['Starting Materials'].append({
                'smiles': node.smiles,
                'cost': node.cost_per_g
            })

        for cpd_id in target_list: 
            node = graph.node_from_id(cpd_id)
            storage['Targets'].append({
                'smiles': node.smiles,
                'reward': node.reward,
            })
        
        with open(output_dir/'solution_list_format.json','w') as f:
            json.dump(storage, f, indent='\t')

        return storage 
    
    def find_rxn_parents(self, store_dict, rxn_id, selected_mols, selected_rxns):

        par_ids = [n.id for n in self.graph.node_from_id(rxn_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_mols)
        for par in selected_pars: 
            store_dict['Compounds'].append(self.graph.smiles_from_id(par))
            store_dict = self.find_mol_parents(store_dict, par, selected_mols, selected_rxns)
        return store_dict

    def find_mol_parents(self, store_dict, mol_id, selected_mols, selected_rxns): 

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
                store_dict['Reactions'].append({
                    'smiles': node.smiles,
                    'conditions': node.get_condition(1)[0], 
                    'score': node.score,
                })
            store_dict = self.find_rxn_parents(store_dict, par, selected_mols, selected_rxns)
        return store_dict

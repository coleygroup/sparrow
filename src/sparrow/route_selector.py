from sparrow import Scorer, Recommender
from sparrow.route_graph import RouteGraph
from sparrow.coster import Coster
from typing import Dict, Union, List
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI, LpInteger
from rdkit import Chem
from tqdm import tqdm
from pathlib import Path
from datetime import datetime

import csv 

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
                 coster: Coster = None, 
                 weights: List = [1,1,1,1],
                 output_dir: str = 'debug'
                 ) -> None:
        
        self.graph = route_graph  
        self.graph.id_nodes()
        self.dir = Path(output_dir)

        self.target_dict = self.clean_target_dict(target_dict)
        self.targets = list(self.target_dict.keys())

        self.rxn_scorer = rxn_scorer
        self.condition_recommender = condition_recommender
              
        self.constrain_all_targets = constrain_all_targets
        self.weights = weights
        
        self.target_dict = self.graph.set_compound_types(self.target_dict, coster=coster)

        if self.condition_recommender is not None: 
            self.get_recommendations()

        if self.rxn_scorer is not None: 
            self.get_rxn_scores()

        self.add_dummy_starting_rxn_nodes()

        self.problem = LpProblem("Route_Selection", LpMinimize)

    def clean_target_dict(self, target_dict: Dict[str, float]) -> Dict[str, float]:
        """ Converts target dict from Dict[smiles, reward] to Dict[id, reward] """
        new_target_dict = {}
        
        c=0
        for old_smi, reward in target_dict.items():
            clean_smi = Chem.MolToSmiles(Chem.MolFromSmiles(old_smi))
            if clean_smi in self.graph.compound_nodes.keys():
                id = self.graph.id_from_smiles(clean_smi)
                new_target_dict[id] = reward
            elif old_smi in self.graph.compound_nodes.keys(): 
                id = self.graph.id_from_smiles(old_smi)
                new_target_dict[id] = reward                
            else: 
                print(f'Target {old_smi} not in routes! Being removed from target set!!!')
                c+=1

        p = self.dir / 'cleaned_tar_dict.csv'
        print(f'Saving remaining targets, ids, and rewards to {p}')

        save_list = [
            {'SMILES': self.graph.smiles_from_id(id), 'ID': id, 'Reward': reward,}
            for id, reward in new_target_dict.items()
        ]
        
        with open(p, 'w') as csvfile: 
            writer = csv.DictWriter(csvfile, fieldnames=['SMILES', 'ID', 'Reward'])
            writer.writeheader()
            writer.writerows(save_list)

        return new_target_dict

    def get_recommendations(self): 
        """ Completes condition recommendation for any reaction node that does not have conditions """
        for node in tqdm(self.graph.reaction_nodes_only(), 'Recommending Conditions'):
            if node.condition_set: 
                continue
            
            condition = self.condition_recommender(node.smiles)
            node.update_condition(condition)
    
    def get_rxn_scores(self): 
        """ Scores all reactions in the graph that are not already scored """
        count = 0
        for node in tqdm(self.graph.reaction_nodes_only(), 'Scoring reactions'): 
            if node.score_set or node.dummy: 
                continue 

            try:    
                score = self.rxn_scorer(rxn_smi=node.smiles, condition=node.condition)
            except: 
                print(f'Reaction {node.smiles} could not be scored, setting score=0')
                score = 0 
                
            node.update(score=score)
            count += 1
            if count % 100 == 0: 
                time = datetime.now().strftime("%H-%M-%S")
                self.graph.to_json(self.dir / 'chkpts' / f'trees_w_scores_{time}.json')
        


    def define_variables(self): 
        """ 
        TODO: explain in readme what variables mean, refer to that here 
        (currently in my thesis proposal)
        TODO: include conditions 
        """

        self.a = LpVariable.dicts(
            "target",
            self.targets,
            cat="Binary",
        ) # whether each target is selected for synthesis 
        
        starting_nodes = self.graph.buyable_nodes()
        self.s = LpVariable.dicts(
            "start",
            [node.id for node in starting_nodes],
            cat="Binary",
        ) # whether each starting material is used in the selected routes 
        
        rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]
        self.o = LpVariable.dicts(
            "rxnfortarget",
            (rxn_ids, self.targets),
            cat="Binary"
        ) # whether each reaction is used to synthesize a given target 

        self.f = LpVariable.dicts(
            "rxnflow",
            rxn_ids,
            lowBound=0,
            upBound=len(self.targets),
            cat=LpInteger,
        )  # flow through reaction node 

        self.r = LpVariable.dicts(
            "rxnused", 
            rxn_ids, 
            cat="Binary",
        )

        return 

    def set_constraints(self):
        """ Sets constraints defined in TODO: write in README all constraints """

        if self.constrain_all_targets: 
            self.set_all_targets_selected_constraint()
        else: 
            self.set_target_flow_constraint()

        self.set_intermediate_flow_constraint()
        self.set_starting_material_constraint()
        self.set_reaction_used_constraint()
        self.set_overall_flow_constraint()

        return 
    
    def set_target_flow_constraint(self):
        """ Sets constraint on flow through a target node """
        # flow into target node - flow out of target node = 0 (target not selected)
        # or 1 (target is selected )
        for target in self.targets: 
            parent_ids, child_ids = self.get_child_and_parent_ids(id=target)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == self.a[target],
                f"Flow_Through_Target_{target}"
            )

        return 
    
    def set_intermediate_flow_constraint(self): 
        """ Sets constraint on flow through intermediate nodes: net flow must be zero """
        intermediate_nodes = [*self.graph.intermediate_nodes(), *self.graph.buyable_nodes()]
        
        for inter in intermediate_nodes: 
            parent_ids, child_ids = self.get_child_and_parent_ids(id=inter.id)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == 0,
                f"Flow_Through_Inter_{inter.id}"
            )
           
        return 
    
    def set_starting_material_constraint(self):
        """ Sets constraint on 's' variables and dummy reaction nodes """

        N = len(self.targets)

        for start in self.graph.buyable_nodes(): 
            parent_ids, _ = self.get_child_and_parent_ids(id=start.id) 
            # ^ parent_smis should only have one dummy rxn node if this is done correctly 
            self.problem += (
                N*self.s[start.id] >= lpSum([self.f[rxn] for rxn in parent_ids]),
                f"Start_{start.id}_from_dummy_flow"
            )

        return 
    
    def set_reaction_used_constraint(self):
        """ Sets constraint on 'r' abd 'f' variables, so if f_rxn > 0, r_rxn = 1"""

        N = len(self.targets)
        for rxn_smi, node in self.graph.reaction_nodes.items(): 
            # ^ parent_smis should only have one dummy rxn node if this is done correctly 
            self.problem += (
                N*self.r[node.id] >= self.f[node.id],
                f"Rxnused_flow_{node.id}"
            )

        return 

    def set_overall_flow_constraint(self):
        """ Sets constraint between o_mn and f_m for reaction node m """

        for rxn in self.graph.reaction_nodes_only(): 
            self.problem += (
                self.f[rxn.id] == lpSum([self.o[rxn.id][target] for target in self.targets]),
                f"Total_flow_through_{rxn.id}"
            )
        
        return 
    
    def set_all_targets_selected_constraint(self):
        """ Sets constraint that all targets are selected """
        for target in self.targets: 
            parent_ids, child_ids = self.get_child_and_parent_ids(target)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == 1,
                f"Flow_Through_Target_{target}"
            )

        return 
    
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
        # FIX THIS 
        rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]

        if self.constrain_all_targets: #TODO: fix this 
            self.problem += self.weights[0]
        else:
            self.problem += -1*self.weights[0]*lpSum([self.target_dict[target]*self.a[target] for target in self.targets]) # reward
            self.problem += self.problem.objective + self.weights[1]*lpSum([s for s in self.s.values()]) # starting materials 
            self.problem += self.problem.objective + self.weights[2]*0 # not considering conditions yet 
            self.problem += self.problem.objective + self.weights[3]*lpSum([self.r[rxn]*self.graph.node_from_id(rxn).penalty for rxn in rxn_ids])
                # reaction penalties 
        return 
    
    def optimize(self, solver=None):

        # self.problem.writeLP("RouteSelector.lp", max_length=300)

        if solver == 'GUROBI': 
            self.problem.solve(GUROBI(timeLimit=86400))
        else: 
            self.problem.solve()

        print("Optimization problem completed...")

        return 
    
    def optimal_variables(self):
        """ Returns nonzero variables """
        nonzero_vars = [
            var for var in self.problem.variables() if var.varValue > 0.01
        ]

        return nonzero_vars

     
    
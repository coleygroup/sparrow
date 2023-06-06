from mars.route_graph import RouteGraph
from typing import Dict, Union, List, Optional
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI, LpInteger
from askcos.synthetic.evaluation.evaluator import Evaluator
from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
import askcos.global_config as gc
from pathlib import Path
import pickle 

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
                 constrain_all_targets: bool = False, 
                 weights: List = [1,1,1,1],
                 calc_reaction_scores: Optional[bool] = True,
                 ) -> None:
        
        self.graph = route_graph        
        self.target_dict = target_dict

        self.rewards = list(target_dict.values())

        self.constrain_all_targets = constrain_all_targets

        self.graph.set_compound_types(target_dict)

        if calc_reaction_scores: 
            context_recommender = NeuralNetContextRecommender()
            context_recommender.load_nn_model(
                model_path=gc.NEURALNET_CONTEXT_REC['model_path'], 
                info_path=gc.NEURALNET_CONTEXT_REC['info_path'], 
                weights_path=gc.NEURALNET_CONTEXT_REC['weights_path']
            )
            self.graph.calc_reaction_scores(
                context_recommender,
                evaluator=Evaluator(),
            )
        
        self.add_dummy_starting_rxn_nodes()

        self.graph.id_nodes()
        self.targets = [self.graph.compound_nodes[tar].id for tar in target_dict.keys()]
        self.target_dict = {self.graph.compound_nodes[tar].id: score for tar, score in target_dict.items()}

        self.problem = LpProblem("Route_Selection", LpMinimize)

        self.weights = weights

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
        
        starting_nodes = self.graph.starting_material_nodes()
        self.s = LpVariable.dicts(
            "start",
            [node.id for node in starting_nodes],
            cat="Binary",
        ) # whether each starting material is used in the selected routes 
        
        rxn_ids = [node.id for node in self.graph.reaction_nodes.values()]
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

        return 

    def set_constraints(self):
        """ Sets constraints defined in TODO: write in README all constraints """

        if self.constrain_all_targets: 
            self.set_all_targets_selected_constraint()
        else: 
            self.set_target_flow_constraint()

        self.set_intermediate_flow_constraint()
        self.set_starting_material_constraint()
        self.set_overall_flow_constraint()

        return 
    
    def set_target_flow_constraint(self):
        """ Sets constraint on flow through a target node """
        # flow into target node - flow out of target node = 0 (target not selected)
        # or 1 (target is selected )
        for target in self.targets: 
            target_smi = self.graph.ids[target].smiles
            parent_ids, child_ids = self.get_compound_child_and_parent_ids(target_smi)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == self.a[target],
                f"Flow_Through_Target_{target}"
            )

        return 
    
    def set_intermediate_flow_constraint(self): 
        """ Sets constraint on flow through intermediate nodes: net flow must be zero """
        intermediate_nodes = [*self.graph.intermediate_nodes(), *self.graph.starting_material_nodes()]
        inter_smiles = [node.smiles for node in intermediate_nodes]
        
        for inter in inter_smiles: 
            parent_ids, child_ids = self.get_compound_child_and_parent_ids(inter)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == 0,
                f"Flow_Through_Inter_{self.graph.compound_nodes[inter].id}"
            )
           
        return 
    
    def set_starting_material_constraint(self):
        """ Sets constraint on 's' variables and dummy reaction nodes """

        starting_smis = [node.smiles for node in self.graph.starting_material_nodes()]
        N = len(self.targets)
        for start in starting_smis: 
            parent_ids, _ = self.get_compound_child_and_parent_ids(start) 
            # ^ parent_smis should only have one dummy rxn node if this is done correctly 
            self.problem += (
                N*self.s[self.graph.compound_nodes[start].id] >= lpSum([self.f[rxn] for rxn in parent_ids]),
                f"Start_{self.graph.compound_nodes[start].id}_from_dummy_flow"
            )

        return 

    def set_overall_flow_constraint(self):
        """ Sets constraint between o_mn and f_m for reaction node m """

        rxn_ids = [node.id for node in self.graph.reaction_nodes.values()]
        for rxn in rxn_ids: 
            self.problem += (
                self.f[rxn] == lpSum([self.o[rxn][target] for target in self.targets]),
                f"Total_flow_through_{rxn}"
            )
        
        return 
    
    def set_all_targets_selected_constraint(self):
        """ Sets constraint that all targets are selected """
        for target in self.targets: 
            target_smi = self.graph.ids[target].smiles
            parent_ids, child_ids = self.get_compound_child_and_parent_ids(target_smi)
            self.problem += (
                lpSum([self.f[parent_id] for parent_id in parent_ids])
                    - lpSum([self.f[child_id] for child_id in child_ids])  
                    == 1,
                f"Flow_Through_Target_{target}"
            )

        return 
    
    def get_compound_child_and_parent_ids(self, node_smis: str): 
        """ Returns list of child node smiles and parent node smiles for a given
        compound smiles """
        child_ids = [child.id for child in self.graph.compound_nodes[node_smis].children.values()] 
        parent_ids = [parent.id for parent in self.graph.compound_nodes[node_smis].parents.values()]

        return parent_ids, child_ids
    
    def add_dummy_starting_rxn_nodes(self): 
        """ Adds reaction nodes that form all starting materials, as described in 
        TODO: describe this in README """
        for start_node in self.graph.starting_material_nodes(): 
            dummy_rxn_smiles = f">>{start_node.smiles}"
            self.graph.add_reaction_node(dummy_rxn_smiles, children=[start_node.smiles], dummy=True)

            # make penalty of dummy reactions zero 
            self.graph.reaction_nodes[dummy_rxn_smiles].penalty = 0
    
    def set_objective(self): 
        # FIX THIS 
        rxn_ids = [node.id for node in self.graph.reaction_nodes.values()]

        if self.constrain_all_targets: #TODO: fix this 
            self.problem += self.weights[0]
        else:
            self.problem += -1*self.weights[0]*lpSum([self.target_dict[target]*self.a[target] for target in self.targets]) 
            self.problem += self.problem.objective + self.weights[1]*lpSum([s for s in self.s.values()])
            self.problem += self.problem.objective + self.weights[2]*0 # not considering conditions yet 
            self.problem += self.problem.objective + self.weights[3]*lpSum([self.o[rxn][target]*self.graph.id_nodes()[rxn].penalty for rxn in rxn_ids for target in self.targets])
                
        return 
    
    def optimize(self):

        self.problem.writeLP("RouteSelector.lp", max_length=300)
        self.problem.solve(GUROBI(timeLimit=86400))
        print("Optimization problem completed...")

        return 
    
    def routes_from_optimal_variables(self):
        """ Takes optimal variables from problem solution and converts it to a set of routes """
        # TODO: implement this 
        nonzero_vars = [
            var for var in self.problem.variables() if var.varValue > 0.01
        ]
        selected_targets = { var.name.split("target_")[1]: self.graph.ids[var.name.split("target_")[1]]
            for var in nonzero_vars
            if var.name.find("target_") == 0 
        } # dict {id: smiles} for selected targets 

        print('Selected targets: ')
        for tar in selected_targets: 
            print(tar)

        return nonzero_vars

     
    
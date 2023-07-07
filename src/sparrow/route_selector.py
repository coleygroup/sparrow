from sparrow.route_graph import RouteGraph
from sparrow import Scorer, Recommender
from typing import Dict, Union, List
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI, LpInteger
from rdkit import Chem

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
                 weights: List = [1,1,1,1],
                 ) -> None:
        
        self.graph = route_graph  
        self.graph.id_nodes()

        self.target_dict = self.clean_target_dict(target_dict)
        self.targets = list(target_dict.keys())
        self.rewards = list(self.target_dict.values())        

        self.rxn_scorer = rxn_scorer
        self.condition_recommender = condition_recommender
              
        self.constrain_all_targets = constrain_all_targets
        self.weights = weights

        if self.condition_recommender is not None: 
            self.get_recommendations()

        if self.rxn_scorer is not None: 
            self.get_rxn_scores()

        self.graph.set_compound_types(self.target_dict)
        self.add_dummy_starting_rxn_nodes()

        self.problem = LpProblem("Route_Selection", LpMinimize)

    def clean_target_dict(self, target_dict: Dict[str, float]) -> Dict[str, float]:
        """ Converts target dict from Dict[smiles, reward] to Dict[id, reward] """
        new_target_dict = {
            Chem.MolToSmiles(Chem.MolFromSmiles(old_smi), isomericSmiles=False) : reward
            for old_smi, reward in target_dict
        }

        return new_target_dict

    def get_recommendations(self): 
        """ Completes condition recommendation for any reaction node that does not have conditions """
        for node in self.graph.reaction_nodes_only():
            if node.conditions_set: 
                continue
            
            conditions = self.condition_recommender(node.smiles)
            node.update_conditions(conditions)
    
    def get_rxn_scores(self): 
        """ Scores all reactions in the graph that are not already scored """
        for node in self.graph.reaction_nodes_only(): 
            if node.score_set: 
                continue 
                
            score = self.rxn_scorer(rxn_smi=node.smiles, conditions=node.conditions)

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
        intermediate_nodes = [*self.graph.intermediate_nodes(), *self.graph.buyable_nodes()]
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

        starting_smis = [node.smiles for node in self.graph.buyable_nodes()]
        N = len(self.targets)
        for start in starting_smis: 
            parent_ids, _ = self.get_compound_child_and_parent_ids(start) 
            # ^ parent_smis should only have one dummy rxn node if this is done correctly 
            self.problem += (
                N*self.s[self.graph.compound_nodes[start].id] >= lpSum([self.f[rxn] for rxn in parent_ids]),
                f"Start_{self.graph.compound_nodes[start].id}_from_dummy_flow"
            )

        return 
    
    def set_reaction_used_constraint(self):
        """ Sets constraint on 'r' abd 'f' variables, so if f_rxn > 0, r_rxn = 1"""

        N = len(self.targets)
        for rxn_smi, node in self.graph.reaction_nodes.items(): 
            # parent_ids, _ = self.get_compound_child_and_parent_ids(rxn_smi) 
            # ^ parent_smis should only have one dummy rxn node if this is done correctly 
            self.problem += (
                N*self.r[node.id] >= self.f[node.id],
                f"Rxnused_flow_{node.id}"
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
        for start_node in self.graph.buyable_nodes(): 
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
            self.problem += -1*self.weights[0]*lpSum([self.target_dict[target]*self.a[target] for target in self.targets]) # reward
            self.problem += self.problem.objective + self.weights[1]*lpSum([s for s in self.s.values()]) # starting materials 
            self.problem += self.problem.objective + self.weights[2]*0 # not considering conditions yet 
            self.problem += self.problem.objective + self.weights[3]*lpSum([self.r[rxn]*self.graph.id_nodes()[rxn].penalty for rxn in rxn_ids])
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
        """ Takes optimal variables from problem solution and converts it to a set of routes """
        nonzero_vars = [
            var for var in self.problem.variables() if var.varValue > 0.01
        ]
        selected_targets = { var.name.split("target_")[1]: self.graph.ids[var.name.split("target_")[1]]
            for var in nonzero_vars
            if var.name.find("target_") == 0 
        } # dict {id: smiles} for selected targets 

        print('Selected targets: ')
        for tar in selected_targets: 
            print(self.graph.ids[tar].smiles)

        return nonzero_vars

     
    
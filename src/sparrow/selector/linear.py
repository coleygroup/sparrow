from typing import Dict, List
from pulp import LpProblem, LpMinimize, LpVariable, lpSum, GUROBI, LpStatus
from pulp.apis import PULP_CBC_CMD
from tqdm import tqdm 
import time 

from sparrow.condition_recommender import Recommender
from sparrow.coster import Coster
from sparrow.route_graph import RouteGraph
from sparrow.scorer import Scorer
from sparrow.selector.base import Selector
from sparrow.utils.cluster_utils import cluster_smiles

class LinearSelector(Selector):
    """ A route selector that uses pulp to formulate and solve the optimization for downselection """
    def __init__(self,
                 route_graph: RouteGraph, 
                 target_dict: Dict[str, float], 
                 rxn_scorer: Scorer = None, 
                 condition_recommender: Recommender = None, 
                 constrain_all_targets: bool = False, 
                 max_targets: int = None, 
                 coster: Coster = None, 
                 weights: List = [1, 1, 1, 0, 0], 
                 output_dir: str = 'debug', 
                 remove_dummy_rxns_first: bool = False, 
                 clusters: dict = None, 
                 N_per_cluster: int = 0,
                 rxn_classes: dict = None,
                 rxn_classifier_dir: str = None,
                 max_rxn_classes: int = None,
                 dont_buy_targets: bool = False,
                 solver: str = 'pulp',
                 max_rxns: int = None, 
                 sm_budget: float = None, 
                 cycle_constraints: bool = False,
                 max_seconds: int = None,
                 extract_routes: bool = True, 
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, weights=weights, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first, max_seconds=max_seconds,
            clusters=clusters, N_per_cluster=N_per_cluster, rxn_classes=rxn_classes,
            rxn_classifier_dir=rxn_classifier_dir, max_rxn_classes=max_rxn_classes, 
            dont_buy_targets=dont_buy_targets, max_rxns=max_rxns, 
            sm_budget=sm_budget, extract_routes=extract_routes,
            cycle_constraints=cycle_constraints
            )
        self.solver = solver 
        
    def initialize_problem(self) -> LpProblem:
        return LpProblem("Route_Selection", LpMinimize)
    
    def define_variables(self):

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

        self.c = None
        if self.rxn_class_dict != None:
            class_ids = self.rxn_classes
            self.c = LpVariable.dicts(
                "class",
                class_ids,
                cat="Binary", # decision var for including a class
            )
        
        return 

    def set_constraints(self):

        print('Setting constraints ...')

        self.set_rxn_constraints()
        self.set_mol_constraints()

        if self.cycle_constraints: 
            self.set_cycle_constraints()
        
        if self.max_targets:
            self.set_max_target_constraint()
        
        if self.constrain_all_targets:
            self.set_constraint_all_targets()
        
        if self.N_per_cluster is not None and self.N_per_cluster > 0:
            self.set_cluster_constraints()

        if self.rxn_class_dict:
            self.set_class_constraints()

        if self.max_rxn_classes:
            print('Setting maximum number of reaction classes')
            self.set_max_classes_constraint()

        if self.sm_budget: 
            self.set_budget_constraint()

        return 
    
    def set_rxn_constraints(self): 

        if self.max_rxns is not None: 
            self.problem += (
                self.max_rxns >= lpSum(self.r[rid] for rid in self.r.keys() if not self.graph.smiles_from_id(rid).startswith('>>'))
            )

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

    def set_max_classes_constraint(self):
        """ Sets constraint on the maximum number of rxn classes used to access selected targets. """
        self.problem += (
            lpSum(self.c) <= self.max_rxn_classes
        )
    
    def set_constraint_all_targets(self): 
        """ Constrains that all targets must be synthesized """
        for target in self.targets: 
            self.problem += (
                self.m[target] == 1
            )

    def set_cluster_constraints(self): 
        # check that self.N_per_cluster >= minimum cluster size 
        min_clus_size = min([len(clus) for clus in self.clusters.values()])
        if self.N_per_cluster > min_clus_size: 
            print(f'Smallest cluster has {min_clus_size} compounds, but N_per_cluster is {self.N_per_cluster}')
            print(f'Switching N_per_cluster to {min_clus_size}')
            self.N_per_cluster = min_clus_size
        
        for cname, ids in self.clusters.items(): 
            self.problem += (
                lpSum(self.m[nid] for nid in ids) >= self.N_per_cluster
            )

    def set_budget_constraint(self):
        dummies_in_r = [dummy.id for dummy in self.graph.dummy_nodes_only() if dummy.id in self.r]
        self.problem += (
            lpSum(self.cost_of_dummy(dummy_id=dummy)*self.r[dummy] for dummy in dummies_in_r) <= self.sm_budget
        )

    def set_class_constraints(self):

        for id in self.rxn_classes:
            rxns = self.rxn_class_dict[id]
            N = len(rxns)
        
            self.problem += (
                N * self.c[id] >= lpSum(self.r[rxn] for rxn in rxns)
            )

            self.problem += (
                self.c[id] <= lpSum(self.r[rxn] for rxn in rxns)
            )

    def set_objective(self, weights = None): 
        print('Setting objective function ...')
        weights = weights or self.weights

        reward_mult = weights[0] # / ( len(self.target_dict)) # *max(self.target_dict.values()) )
        cost_mult = weights[1] # / (len(self.graph.dummy_nodes_only())) # * max([node.cost_per_g for node in self.graph.buyable_nodes()]) ) 
        pen_mult = weights[2] # / (len(self.graph.non_dummy_nodes())) # * max([node.penalty for node in self.graph.non_dummy_nodes()]) )

        self.problem += -1*reward_mult*lpSum([float(self.target_dict[target])*self.m[target] for target in self.targets]) \
        + cost_mult*lpSum([self.cost_of_dummy(dummy)*self.r[dummy.id] for dummy in self.graph.dummy_nodes_only()]) \
        + pen_mult*lpSum([self.r[node.id]*float(node.penalty) for node in self.graph.non_dummy_nodes()])
            # reaction penalties, implement CSR later 
        
        if weights[3]>0: 
            self.add_diversity_objective()

        if self.c != None and len(weights) > 4 and weights[4] > 0: 
            self.add_rxn_class_objective()

        return 

    def add_diversity_objective(self): 
        """ Adds scalarization objective to increase the number of clusters represented, requires defining new variable """

        # d_i : whether cluster i is represented by the selected set 
        self.d = LpVariable.dicts(
            "cluster", 
            indices=range(len(self.clusters)), 
            cat="Binary",
        )

        # constraint: d_i <= sum(c_j) for j in cluster i
        for i, ids_in_cluster in enumerate(self.clusters.values()): 
            self.problem += (
                self.d[i] <= lpSum(self.m[cpd_id] for cpd_id in ids_in_cluster)
            )
        
        # add objective 
        print(f'adding diversity objective with weight {self.weights[3]}')
        self.problem += self.problem.objective - self.weights[3]*lpSum(self.d)

        return 
    
    def add_rxn_class_objective(self): 
        """ Adds objective to decrease the number of reaction classes used """
        
        # add objective 
        print(f'Adding reaction class objective with weight {self.weights[4]}')
        self.problem += self.problem.objective + self.weights[4]*lpSum(self.c)

        return
    
    def optimize(self, output_dir=None):

        # self.problem.writeMPS(self.dir/'chkpts'/'problem.mps')

        print("Solving optimization problem ...")
        opt_start = time.time()
        if self.solver == 'GUROBI' or self.solver == 'gurobi': 
            self.problem.solve(GUROBI(timeLimit=self.max_seconds))
        else: 
            self.problem.solve(PULP_CBC_CMD(gapRel=1e-7, gapAbs=1e-9, msg=False, timeLimit=self.max_seconds))

        if self.problem.status == -1: 
            raise RuntimeError('Problem is infeasible. To fix, try relaxing constraints.')
        
        self.runtime = time.time()-opt_start

        print(f"Optimization problem completed. Took {self.runtime:0.2f} seconds to solve")
        
        return self
    
    def extract_selected_ids(self):
        """ Returns nonzero variables names """
        
        nonzero_vars = [
            var for var in self.problem.variables() if var.varValue > 0.01
        ]

        rxn_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('rxn')]
        mol_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('mol')]
        class_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('class')]

        return mol_ids, rxn_ids, class_ids
    
    def get_num_variables(self): 
        return self.problem.numVariables()
    
    def get_num_constraints(self):
        return self.problem.numConstraints()
         
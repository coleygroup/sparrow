import numpy as np 
import gurobipy as gp 
import time 
import re 
from gurobipy import GRB
from typing import Dict
from tqdm import tqdm 

from sparrow.selector.base import Selector
from sparrow.scorer import Scorer
from sparrow.condition_recommender import Recommender
from sparrow.route_graph import RouteGraph
from sparrow.coster import Coster



class ExpectedRewardSelector(Selector):
    """ 
    A route selector that uses Gurobi to optimize expected reward given reaction 
    and starting material budget constraints 
    
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
                 output_dir: str = 'debug', 
                 remove_dummy_rxns_first: bool = False, 
                 clusters: dict = None, 
                 N_per_cluster: int = 1, 
                 max_rxns: int = None, 
                 sm_budget: float = None, 
                 dont_buy_targets: bool = False
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, cost_per_rxn=cost_per_rxn, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first,
            clusters=clusters, max_rxns=max_rxns, sm_budget=sm_budget, 
            dont_buy_targets=dont_buy_targets, N_per_cluster=N_per_cluster,
            )
        
    def initialize_problem(self) -> gp.Model:
        return gp.Model("sparrow")

    def define_variables(self): 

        # whether rxn is used 
        self.rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]
        self.r = self.problem.addVars(
            self.rxn_ids,
            vtype=GRB.BINARY,
            name="rxn",
        )

        # whether molecule is selected 
        mol_ids = [node.id for node in self.graph.compound_nodes_only()]
        self.m = self.problem.addVars(
            mol_ids,
            vtype=GRB.BINARY,
            name="mol"
        )

        # whether rxn i is used for target j 
        self.u = self.problem.addVars(
            self.rxn_ids, self.targets, 
            vtype=GRB.BINARY,
            name="rxnfor"
        )

        # whether compound j is used int he route to target j 
        self.cpdfor = self.problem.addVars(
            mol_ids, self.targets,
            vtype=GRB.BINARY,
            name="cpdfor"
        )

        # additional variables for nonlinearity 
        sum_log_scores = np.log([node.score if node.score>0 else 1e-10 for node in self.graph.non_dummy_nodes()]).sum()
        self.usum = self.problem.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,  
            name="usum",
            ub=0,
            lb=sum_log_scores,
        )
        self.expusum = self.problem.addVars(
            self.targets,
            vtype=GRB.CONTINUOUS,
            name="expusum",
            ub=1,
            lb=0
        )
        self.probsuccess = self.problem.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,
            name="probsuccess",
            ub=1,
            lb=0,
        )
        max_n_rxns = len(self.graph.non_dummy_nodes())
        max_sm_costs = np.sum([self.cost_of_dummy(dummy) for dummy in self.graph.dummy_nodes_only()])
        max_costs = max_n_rxns*self.cost_per_rxn + max_sm_costs
        self.allcosts = self.problem.addVar(
            name="allcosts",
            vtype=GRB.CONTINUOUS,
            ub=max_costs,
            lb=0,
        )
        max_er = sum(self.target_dict.values())
        self.sumer = self.problem.addVar(
            name="sumer",
            vtype=GRB.CONTINUOUS,
            ub=max_er,
            lb=0,
        )

        return 
        
    def set_constraints(self, set_cycle_constraints=True):
        """ Sets constraints """

        print('Setting constraints ...')

        self.set_rxn_constraints()
        self.set_mol_constraints()

        if set_cycle_constraints: 
            self.set_cycle_constraints()
        
        if self.max_targets:
            self.set_max_target_constraint()
        
        if self.constrain_all_targets:
            self.set_constraint_n_targets(N=len(self.targets))

        if self.N_per_cluster is not None and self.N_per_cluster > 0:
            self.set_cluster_constraints()

        if self.sm_budget: 
            self.set_budget_constraint()

        self.set_nonlinear_constraints()

        return 
    
    def set_nonlinear_constraints(self): 

        log_scores = np.log([node.score if node.score>0 else 1e-10 for node in self.graph.non_dummy_nodes()])
        
        for t in self.targets: 
            self.problem.addConstr(
                self.usum[t] == gp.quicksum(
                    self.u[node.id, t]*logscore 
                    for logscore, node in zip(log_scores, self.graph.non_dummy_nodes())
                ),
                name=f'usumConstr_{t}',
            )
            self.problem.addGenConstrExp(self.usum[t], self.expusum[t], name=f'expusumConstr_{t}',options="FuncPieces=-2 FuncPieceError=0.00001")
            self.problem.addQConstr(
                self.probsuccess[t], GRB.LESS_EQUAL, self.expusum[t]*self.cpdfor[t, t], 
                name=f'probsuccessConstr_{t}',
            )
        
        self.problem.addQConstr(
            self.sumer, GRB.LESS_EQUAL, gp.quicksum(
                self.target_dict[t]*self.probsuccess[t] for t in self.targets
            ),
            name='sumER_constr',
        )

        return 
    
    def set_rxn_constraints(self): 
        # TODO: consider redudancy in these constraints and simplify problem 
        
        if self.max_rxns: 
            self.problem.addConstr(
                self.max_rxns >= gp.quicksum(self.r[node.id] for node in self.graph.non_dummy_nodes())
            )

        for node in tqdm(self.graph.reaction_nodes_only(), desc='Reaction constraints'):  
            par_ids = [par.id for par in node.parents.values()]
            for par_id in par_ids:
                # whether reaction selected at all
                self.problem.addConstr(
                    self.m[par_id] >= self.r[node.id],
                    name=f'rxnConstr_{node.id}_{par_id}'
                )
                
                # if a reaction is selected, every parent (reactant) must also be selected 
                self.problem.addConstrs(
                    self.cpdfor[par_id, t] >= self.u[node.id, t]
                    for t in self.targets
                )

            # if a reaction is used for a target, it is also used in general 
            self.problem.addConstr(
                len(self.targets)*self.r[node.id] >= gp.quicksum(self.u[node.id, tar] for tar in self.targets),
                name=f'rxnusedConstr_{node.id}'
            )

            self.problem.addConstr(
                self.r[node.id] <= gp.quicksum(self.u[node.id, tar] for tar in self.targets),
                name=f'rxnusedConstr2_{node.id}'
            )

        return 
    
    def set_mol_constraints(self): 
        
        for t in self.targets: 
            self.problem.addConstr(
                self.cpdfor[t, t] == self.m[t]
            )

        for node in tqdm(self.graph.compound_nodes_only(), 'Compound constraints'): 
            par_ids = [par.id for par in node.parents.values()]
            self.problem.addConstr(
                self.m[node.id] <= gp.quicksum(self.r[par_id] for par_id in par_ids),
                name=f'cpdConstr_{node.id}'
            )

            # if a compound is used for a target, at least one parent reaction must also be used for that target 
            self.problem.addConstrs(
                self.cpdfor[node.id, t] <= gp.quicksum(self.u[par_id, t] for par_id in par_ids) 
                for t in self.targets
            )

            # if a compound is used for a target, it is also used in general 
            self.problem.addConstr(
                len(self.targets)*self.m[node.id] >= gp.quicksum(self.cpdfor[node.id, tar] for tar in self.targets),
                name=f'molusedConstr_{node.id}'
            )     

            self.problem.addConstr(
                self.m[node.id] <= gp.quicksum(self.cpdfor[node.id, tar] for tar in self.targets),
                name=f'molusedConstr_{node.id}'
            )      
        
        return 

    def set_cycle_constraints(self): 

        cycles = self.graph.dfs_find_cycles_nx()
        c = 0
        for cyc in tqdm(cycles, desc='Cycle constraints'): 
            self.problem.addConstr(
                gp.quicksum(self.r[rid] for rid in cyc) <= (len(cyc) - 1),
                name=f'cycConstr_{c}'
            )
            c += 1

        return 
    
    def set_max_target_constraint(self):
        """ Sets constraint on the maximum number of selected targets. """
        self.problem.addConstr(
            gp.quicksum(self.m[target] for target in self.targets) <= self.max_targets,
            name='max_targets'
        )
    
    def set_constraint_n_targets(self, N): 
        """ Constrains that n targets must be synthesized """
        
        self.problem.addConstr(
            gp.quicksum(self.m[target] for target in self.targets) == N,
            name='constrain_n'
        )

        return 'constrain_n'

    def set_budget_constraint(self): 
        self.problem.addConstr(
            gp.quicksum(self.cost_of_dummy(dummy)*self.r[dummy.id] for dummy in self.graph.dummy_nodes_only()) <= self.sm_budget
        )

    def set_cluster_constraints(self): 
        """ 
        Constrains that every cluster must be represented by the selected set of candidates 
        clusters: a list of lists
                  each list represents a cluster, and each element of each list is a CompoundNode ID 
                  corresponding to the target in that cluster 
        """
        # check that self.N_per_cluster >= minimum cluster size 
        min_clus_size = min([len(clus) for clus in self.clusters.values()])
        if self.N_per_cluster > min_clus_size: 
            print(f'Smallest cluster has {min_clus_size} compounds, but N_per_cluster is {self.N_per_cluster}')
            print(f'Switching N_per_cluster to {min_clus_size}')
            self.N_per_cluster = min_clus_size
        i=0
        for cname, ids in self.clusters.items(): 
            if i > 12: 
                continue 
            self.problem.addConstr(
                gp.quicksum(self.m[nid] for nid in ids) >= self.N_per_cluster,
                name=f'cluster_represented_{cname}'
            )
            i+=1
    
    def set_objective(self, cost_weight=0): 

        print('Setting objective function ...')

        self.problem.ModelSense = GRB.MAXIMIZE

        if cost_weight > 0:
            self.problem.setObjective(self.sumer - cost_weight*self.allcosts)
        else: 
            self.problem.setObjective(self.sumer)

        return 
    
    def optimize(self, savedir=None, solver=None):

        if not savedir: 
            savedir = self.dir 

        print("Solving optimization problem...")
        opt_start = time.time()

        self.problem.params.NonConvex = 2
        self.problem.Params.TIME_LIMIT = 3*3600
        self.problem.optimize()

        if self.problem.status == 3: 
            raise RuntimeError('Problem is infeasible. To fix, try relaxing constraints.')
    
        print(f"Optimization problem completed. Took {time.time()-opt_start:0.2f} seconds.")
        
        return 
    
    def extract_selected_ids(self):
        nonzero_varnames = [
            var.VarName for var in self.problem.getVars() if var.X > 0.01
        ]
        rxn_ids = re.findall(r'rxn\[(.*?)\]', ' '.join(nonzero_varnames))
        mol_ids = re.findall(r'mol\[(.*?)\]', ' '.join(nonzero_varnames))
        return mol_ids, rxn_ids
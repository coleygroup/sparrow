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
                 N_per_cluster: int = 0, 
                 max_rxns: int = None, 
                 sm_budget: float = None, 
                 dont_buy_targets: bool = False,
                 cycle_constraints: bool = True,
                 max_seconds: int = None,
                 rxn_classes: dict = None,
                 max_rxn_classes: int = None,
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, cost_per_rxn=cost_per_rxn, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first,
            clusters=clusters, max_rxns=max_rxns, sm_budget=sm_budget, 
            dont_buy_targets=dont_buy_targets, N_per_cluster=N_per_cluster,
            cycle_constraints=cycle_constraints, max_seconds=max_seconds,
            rxn_classes=rxn_classes, max_rxn_classes=max_rxn_classes,
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
        self.mol_ids = [node.id for node in self.graph.compound_nodes_only()]
        self.m = self.problem.addVars(
            self.mol_ids,
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
            self.mol_ids, self.targets,
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
        
        self.c = None
        if self.rxn_class_dict != None:
            self.c = self.problem.addVars(
                self.rxn_classes,
                vtype=GRB.BINARY,
                name="class", # decision var for including a class
            )

        return 
        
    def set_constraints(self):
        """ Sets constraints """

        print('Setting constraints ...')

        self.set_rxn_constraints()
        self.set_mol_constraints()

        if self.cycle_constraints: 
            self.set_cycle_constraints()
        
        if self.max_targets:
            self.set_max_target_constraint()
        
        if self.constrain_all_targets:
            self.set_constraint_n_targets(N=len(self.targets))

        if self.N_per_cluster is not None and self.N_per_cluster > 0:
            self.set_cluster_constraints()

        if self.sm_budget: 
            self.set_budget_constraint()

        if self.rxn_class_dict:
            self.set_class_constraints()

        if self.max_rxn_classes:
            self.set_max_classes_constraint()

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
        
        if self.max_rxns is not None: 
            self.problem.addConstr(
                self.max_rxns >= gp.quicksum(self.r[rid] for rid in self.r.keys() if not self.graph.smiles_from_id(rid).startswith('>>'))
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

        for t in self.targets: 
            self.problem.addConstr(
                gp.quicksum(self.u[r, t] for r in self.rxn_ids) <= len(self.rxn_ids)*self.m[t]
            )
            self.problem.addConstr(
                gp.quicksum(self.cpdfor[c, t] for c in self.mol_ids) <= len(self.mol_ids)*self.m[t]
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
            if all([rid in self.r for rid in cyc]):
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
        dummies_in_r = [dummy.id for dummy in self.graph.dummy_nodes_only() if dummy.id in self.r]
        self.problem.addConstr(
            gp.quicksum(self.cost_of_dummy(dummy_id=dummy)*self.r[dummy] for dummy in dummies_in_r) <= self.sm_budget
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

        for cname, ids in self.clusters.items(): 
            self.problem.addConstr(
                gp.quicksum(self.m[nid] for nid in ids) >= self.N_per_cluster,
                name=f'cluster_represented_{cname}'
            )
    
    def set_class_constraints(self):
        
        for id in tqdm(self.rxn_classes, 'Reaction class constraints'):
            rxns = self.rxn_class_dict[id]
            N = len(rxns)
        
            self.problem.addConstr(
                N * self.c[id] >= gp.quicksum(self.r[rxn] for rxn in rxns if rxn in self.r)
            )

            self.problem.addConstr(
                self.c[id] <= gp.quicksum(self.r[rxn] for rxn in rxns if rxn in self.r)
            )

        return 

    def set_max_classes_constraint(self): 
        self.problem.addConstr(
            gp.quicksum(self.c) <= self.max_rxn_classes
        )

    def set_objective(self, cost_weight=0): 

        print('Setting objective function ...')

        self.problem.ModelSense = GRB.MAXIMIZE

        if cost_weight > 0:
            self.problem.setObjective(self.sumer - cost_weight*self.allcosts)
        else: 
            self.problem.setObjective(self.sumer)

        return 
    
    def optimize(self, savedir=None):

        if not savedir: 
            savedir = self.dir 

        print("Solving optimization problem ...")
        opt_start = time.time()

        self.problem.Params.NonConvex = 2
        self.problem.Params.TIME_LIMIT = self.max_seconds
        self.problem.optimize()

        self.runtime = time.time()-opt_start

        if self.problem.status == 3: 
            raise RuntimeError('Problem is infeasible. To fix, try relaxing constraints.')
    
        print(f"Optimization problem completed. Took {self.runtime:0.2f} seconds.")
        
        return 
    
    def extract_selected_ids(self):
        nonzero_varnames = [
            var.VarName for var in self.problem.getVars() if var.X > 0.01
        ]
        rxn_ids = re.findall(r'rxn\[(.*?)\]', ' '.join(nonzero_varnames))
        mol_ids = re.findall(r'mol\[(.*?)\]', ' '.join(nonzero_varnames))
        return mol_ids, rxn_ids, []

    def find_rxn_parents(self, store_dict, rxn_id, selected_mols, selected_rxns, target=None):

        par_ids = [n.id for n in self.graph.node_from_id(rxn_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_mols)
        for par in selected_pars: 
            if self.cpdfor[par, target].X == 1:
                store_dict['Compounds'].append(self.graph.smiles_from_id(par))
                store_dict = self.find_mol_parents(store_dict, par, selected_mols, selected_rxns, target=target)
        return store_dict

    def find_mol_parents(self, store_dict, mol_id, selected_mols, selected_rxns, target=None): 

        par_ids = [n.id for n in self.graph.node_from_id(mol_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_rxns)
        for par in selected_pars: 
            if self.u[par, target].X == 1:
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

                store_dict = self.find_rxn_parents(store_dict, par, selected_mols, selected_rxns, target=target)
        return store_dict

    def get_num_variables(self): 
        return len(self.problem.getVars())
    
    def get_num_constraints(self):
        return len(self.problem.getConstrs())

class PrunedERSelector(ExpectedRewardSelector): 
    """ 
    This is an ExpectedRewardSelector that reduces the total number of decision variables
    that are defined in the optimization. Specifically, the ExpectedRewardSelector defines 
    the decision variable self.u to have Nr (number of reactions in the graph) * Nt (number of targets) 
    elements and the decision array self.cpdfor to have Nc (number of compounds in the graph) * Nt 
    elements. For large graphs this is prohibitive. A PrunedERSelector instead identifies which 
    nodes are within some distance (prune_distance) away from each target and only defines the 
    relevant decision variables for those connected nodes. For low prune_distance, this reduces the 
    time to define the problem and to optimize. For large prune_distance, the time to define the 
    problem may increase, but solving time will either remain the same or decrease. 
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
                 N_per_cluster: int = 0, 
                 max_rxns: int = None, 
                 sm_budget: float = None, 
                 dont_buy_targets: bool = False,
                 prune_distance: int = 10, 
                 cycle_constraints: bool = True,
                 max_seconds: int = None,
                 rxn_classes: dict = None,
                 max_rxn_classes: int = None,
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, cost_per_rxn=cost_per_rxn, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first, N_per_cluster=N_per_cluster,
            clusters=clusters, max_rxns=max_rxns, sm_budget=sm_budget, 
            dont_buy_targets=dont_buy_targets, cycle_constraints=cycle_constraints,
            max_seconds=max_seconds, rxn_classes=rxn_classes, max_rxn_classes=max_rxn_classes,
            )  

        self.prune_distance = prune_distance 
        
    def define_variables(self): 
     
        self.rxn_ids = [node.id for node in self.graph.reaction_nodes_only()]
        self.mol_ids = [node.id for node in self.graph.compound_nodes_only()]

        # whether rxn i is used for target j 
        # whether compound j is used in the route to target j 
        # only define this variable for reactions that are within some directed distance of targets
        self.u = {}
        self.cpdfor = {}

        self.u_by_targets = {t_id: [] for t_id in self.targets}
        self.cpdfor_by_targets = {t_id: [] for t_id in self.targets}

        for t_id in tqdm(self.targets, desc='"for" constraints'): 
            connected_nodes = self.graph.get_connected_nodes(t_id, max_distance=self.prune_distance)
            for n in connected_nodes: 
                n_id = n.id
                if n_id.startswith('C'):
                    if n_id not in self.cpdfor: 
                        self.cpdfor[n_id] = {}
                    self.cpdfor[n_id][t_id] = self.problem.addVar(
                        vtype=GRB.BINARY,
                        name=f'cpdfor[{n_id}][{t_id}]'
                    )
                    if n_id not in self.cpdfor_by_targets: 
                        self.cpdfor_by_targets[t_id].append(n_id)
                elif n_id.startswith('R'): 
                    if n_id not in self.u: 
                        self.u[n_id] = {}
                    self.u[n_id][t_id] = self.problem.addVar(
                        vtype=GRB.BINARY,
                        name=f'rxnfor[{n_id}][{t_id}]'
                    )
                    if n_id not in self.u_by_targets[t_id]:
                        self.u_by_targets[t_id].append(n_id)
            
        # whether rxn is used, only worth defining if in self.u
        self.r = self.problem.addVars(
            self.u.keys(),
            vtype=GRB.BINARY,
            name="rxn",
        )

        # whether molecule is selected, only worth defining if in self.cpdfor
        self.m = self.problem.addVars(
            self.cpdfor.keys(),
            vtype=GRB.BINARY,
            name="mol"
        )

        # additional variables for nonlinearity 
        usum_lower_bounds = [
            sum(sorted(np.log([
                self.graph.node_from_id(r_id).score for r_id in self.u_by_targets[t_id]
            ]))[:self.max_rxns])
            for t_id in self.targets
        ]
        self.usum = self.problem.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,  
            name="usum",
            ub=0,
            lb=usum_lower_bounds,
        )
        self.expusum = self.problem.addVars(
            self.targets,
            vtype=GRB.CONTINUOUS,
            name="expusum",
            ub=1,
            lb=0,
        )
        self.probsuccess = self.problem.addVars(
            self.targets, 
            vtype=GRB.CONTINUOUS,
            name="probsuccess",
            ub=1,
            lb=0,
        )

        if self.sm_budget: 
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

        if self.rxn_class_dict != None:
            self.c = self.problem.addVars(
                self.rxn_classes,
                vtype=GRB.BINARY,
                name="class", # decision var for including a class
            )
        
        return 

    def set_rxn_constraints(self):
        # TODO: consider redudancy in these constraints and simplify problem 
        
        if self.max_rxns is not None: 
            self.problem.addConstr(
                self.max_rxns >= gp.quicksum(self.r[rid] for rid in self.r.keys() if not self.graph.smiles_from_id(rid).startswith('>>'))
            )

        for rid in tqdm(self.u.keys(), desc='Reaction constraints'): 
            
            node = self.graph.node_from_id(rid)
            par_ids = [par.id for par in node.parents.values()]
            for par_id in par_ids:
                # whether reaction selected at all
                if par_id in self.m:
                    self.problem.addConstr(
                        self.m[par_id] >= self.r[node.id],
                        name=f'rxnConstr_{node.id}_{par_id}'
                    )
                
                # if a reaction is selected, every parent (reactant) must also be selected 
                if rid in self.u and par_id in self.cpdfor:
                    ts = self.cpdfor[par_id].keys() & self.u[node.id].keys()
                    self.problem.addConstrs(
                        self.cpdfor[par_id][t] >= self.u[node.id][t] for t in ts
                    )

            # if a reaction is used, its child compound(s) must also be selected 
            child_ids = [child.id for child in node.children.values()]
            self.problem.addConstrs(
                self.r[rid] <= self.m[c_id] for c_id in child_ids
            )

            # if a reaction is used for a target, it is also used in general 
            if node.id in self.u: 
                
                self.problem.addConstr(
                    len(self.u[node.id])*self.r[node.id] >= gp.quicksum(self.u[node.id].values()),
                    name=f'rxnusedConstr_{node.id}'
                )

                self.problem.addConstr(
                    self.r[node.id] <= gp.quicksum(self.u[node.id].values()),
                    name=f'rxnusedConstr2_{node.id}'
                )

        for t in tqdm(self.targets, desc='Reaction constraints for targets'): 
            # rs = [r for r in self.rxn_ids if r in self.u and t in self.u[r]]
            rs = self.u_by_targets[t]
            self.problem.addConstr(
                gp.quicksum(self.u[r][t] for r in rs) <= len(rs)*self.m[t]
            )

            # cs = [c for c in self.mol_ids if c in self.cpdfor and t in self.cpdfor[c]]
            cs = self.cpdfor_by_targets[t]
            self.problem.addConstr(
                gp.quicksum(self.cpdfor[c][t] for c in cs) <= len(cs)*self.m[t]
            )

        return 
    
    def set_mol_constraints(self): 
        
        for t in self.targets: 
            self.problem.addConstr(
                self.cpdfor[t][t] == self.m[t]
            )

        for cid in tqdm(self.m.keys(), 'Compound constraints'): 
            node = self.graph.node_from_id(cid)
            par_ids = [par.id for par in node.parents.values() if par.id in self.r]
            self.problem.addConstr(
                self.m[cid] <= gp.quicksum(self.r[par_id] for par_id in par_ids),
                name=f'cpdConstr_{node.id}'
            )

            # if a compound is used for a target, at least one parent reaction must also be used for that target 
            if cid in self.cpdfor:
                for t in self.cpdfor[cid].keys():
                    par_idss = [par_id for par_id in par_ids if par_id in self.u and t in self.u[par_id]]
                    self.problem.addConstr(
                        self.cpdfor[cid][t] <= gp.quicksum(self.u[par_id][t] for par_id in par_idss) 
                    )

                # if a compound is used for a target, it is also used in general 
                self.problem.addConstr(
                    len(self.cpdfor[cid])*self.m[cid] >= gp.quicksum(self.cpdfor[cid].values()),
                    name=f'molusedConstr_{cid}'
                )     

                self.problem.addConstr(
                    self.m[cid] <= gp.quicksum(self.cpdfor[cid].values()),
                    name=f'molusedConstr_{cid}'
                )      

            # if a non-target compound is selected, it must have at least one child reaction be selected
            if cid not in self.targets: 
                child_ids = [child.id for child in node.children.values() if child.id in self.r]
                self.problem.addConstr(
                    self.m[cid] <= gp.quicksum(self.r[r_id] for r_id in child_ids)
                )

        return 
    
    def set_nonlinear_constraints(self):

        log_scores = np.log([max(node.score, 1e-3) for node in self.graph.non_dummy_nodes()])
        log_scores = {node.id: score for node, score in zip(self.graph.reaction_nodes_only(), log_scores)}

        for t in tqdm(self.targets, 'Nonlinear constraints'): 
            nids = [nid for nid in self.u_by_targets[t] if nid in log_scores]
                    #[node.id for node in self.graph.non_dummy_nodes()
                   # if node.id in self.u and t in self.u[node.id]]
            self.problem.addConstr(
                self.usum[t] == gp.quicksum(
                    self.u[nid][t]*log_scores[nid] 
                    for nid in nids
                ),
                name=f'usumConstr_{t}',
            )
            self.problem.addGenConstrExp(self.usum[t], self.expusum[t], name=f'expusumConstr_{t}')
            self.problem.addQConstr(
                self.probsuccess[t], GRB.EQUAL, self.expusum[t]*self.cpdfor[t][t], 
                name=f'probsuccessConstr_{t}',
            )
        
        self.problem.addQConstr(
            self.sumer, GRB.EQUAL, gp.quicksum(
                self.target_dict[t]*self.probsuccess[t] for t in self.targets
            ),
            name='sumER_constr',
        )

        return 
    
    def find_rxn_parents(self, store_dict, rxn_id, selected_mols, selected_rxns, target=None):

        par_ids = [n.id for n in self.graph.node_from_id(rxn_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_mols)
        for par in selected_pars: 
            if target in self.cpdfor[par] and self.cpdfor[par][target].X>0: # TODO FIX 
                store_dict['Compounds'].append(self.graph.smiles_from_id(par))
                store_dict = self.find_mol_parents(store_dict, par, selected_mols, selected_rxns, target=target)
        return store_dict

    def find_mol_parents(self, store_dict, mol_id, selected_mols, selected_rxns, target=None): 

        par_ids = [n.id for n in self.graph.node_from_id(mol_id).parents.values()]
        selected_pars = set(par_ids) & set(selected_rxns)
        for par in selected_pars: 
            if target in self.u[par] and self.u[par][target].X>0: # TODO FIX
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

                store_dict = self.find_rxn_parents(store_dict, par, selected_mols, selected_rxns, target=target)
        return store_dict

class QuantityAwareERSelector(Selector): 
    """ 
    A selector that optimizes expected reward, but approximates the amount of 
    SM needed in the cost function. Assumes that the minimum required quanttiy of starting material j 
    equals T*m, where T is the number of targets whose routes use the starting material j and m is 
    a constant parameter with units mass/target. This class also requires a cost function for 
    starting materials, i.e., a list of quantities and corresponding costs. This can be 
    addressed using the QuantityLookupCoster, or by including the cost function in the 
    provided graph. 
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
                 cluster_cutoff: float = 0.7, 
                 custom_clusters: dict = None, 
                 max_rxns: int = None, 
                 sm_budget: float = None, 
                 dont_buy_targets: bool = False
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, cost_per_rxn=cost_per_rxn, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first, cluster_cutoff=cluster_cutoff, 
            custom_clusters=custom_clusters, max_rxns=max_rxns, sm_budget=sm_budget, 
            dont_buy_targets=dont_buy_targets
            )
        

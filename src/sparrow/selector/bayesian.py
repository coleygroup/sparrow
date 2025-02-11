from typing import Dict, List
from skopt import gp_minimize
from pathlib import Path
import json
import shutil
import math 
import time 

from sparrow.condition_recommender import Recommender
from sparrow.coster import Coster
from sparrow.route_graph import RouteGraph
from sparrow.scorer import Scorer
from sparrow.selector.linear import LinearSelector

class BOLinearSelector(LinearSelector):
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
                 bayes_iters: int = 20,
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, weights=weights, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first, max_seconds=max_seconds, 
            clusters=clusters, N_per_cluster=N_per_cluster, rxn_classes=rxn_classes,
            rxn_classifier_dir=rxn_classifier_dir, max_rxn_classes=max_rxn_classes, 
            dont_buy_targets=dont_buy_targets, cycle_constraints=cycle_constraints,
            max_rxns=max_rxns, sm_budget=sm_budget, extract_routes=extract_routes,
            )

        self.bayes_iters = bayes_iters
        self.solver = solver 
    
    def formulate_and_optimize(self, **kwargs):
        self.define_variables()
        self.set_constraints()
        opt_start_time = time.time()
        res = self.optimize_weights()
        bo_run_time = time.time() - opt_start_time
        r_weight = res.x[0]
        print(f'Optimal weighting factors:\n\tReward weight: {r_weight:0.3e}\n\tReaction weight: {1-r_weight:0.3e}')
        print(f'Yielding expected reward of {res.fun * -1:0.3f}')

        best_output_dir = self.dir / 'tuning_results' / f'lambda_{r_weight:0.3f}'
        dest = self.dir/'BEST_SOLUTION'
        if dest.exists(): 
            print(f'Overwriting directory: {dest}')
            shutil.rmtree(dest)
        shutil.copytree(best_output_dir, self.dir/'BEST_SOLUTION')  
        
        with open(self.dir/'BEST_SOLUTION'/'summary.json', 'r') as f:
            final_summary = json.load(f) 
        
        final_summary['Total Tuning Run Time'] = bo_run_time
        final_summary['Reward weighting factor in each tuning iteration'] = res.x_iters
        final_summary['Function values in each tuning iteration'] = list(-1*res.func_vals) 
        with open(self.dir/'BEST_SOLUTION'/'summary.json', 'w') as f:
            final_summary = json.dump(final_summary, f, indent='\t') 

        return self
    
    def run_opt_vary_weights(self, weights: list, output_dir: str, extract_routes: bool = True):
        self.set_objective(weights=weights)
        self.optimize()
        summary = self.post_processing(extract_routes=extract_routes, output_dir=output_dir)
        summary['Weights'] = weights
        with open(output_dir/'summary.json', 'w') as f:
            json.dump(summary, f, indent='\t') 

        return summary

    def expected_reward(self, reward_weight):
        all_weights = [reward_weight[0], 0, 1-reward_weight[0], 0, 0]
        print(f'Optimizing with weights: {all_weights}')
        (self.dir / 'tuning_results').mkdir(exist_ok=True, parents=True)
        output_dir = self.dir / 'tuning_results' / f'lambda_{reward_weight[0]:0.3f}'
        output_dir.mkdir(exist_ok=True)
        summary = self.run_opt_vary_weights(weights=all_weights, output_dir=output_dir)
        return -1 * summary['Expected Reward']
        
    def optimize_weights(self, random_state=0):
        res = gp_minimize(
            self.expected_reward,    # the function to minimize
            [(1e-5,1-1e-5)],      # the bounds on each dimension of x
            acq_func="EI",      # the acquisition function
            n_calls=self.bayes_iters,         # the number of evaluations of f
            n_initial_points=2, # min(10, math.floor(self.bayes_iters/5)),  # the number of random initialization points
            random_state=random_state, # random state for reproducibility 
        )

        return res
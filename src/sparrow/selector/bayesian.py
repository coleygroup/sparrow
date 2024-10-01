from typing import Dict, List
from pulp import lpSum
from skopt import gp_minimize
from pathlib import Path
import json

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
                 post_opt_class_score: str = None,
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
            post_opt_class_score=post_opt_class_score,
            )

        self.bayes_iters = bayes_iters
        self.solver = solver 
    
    def optimize(self):
        res = self.optimize_weights()
        print("BAYESIAN RESULT:" + str(res.fun * -1) + " at " + str(res.x))
        final_summary = {}
        result_weights = [res.x[0], 0, res.x[1], self.weights[3], self.weights[4]] # reconstructing all_weights

        summary = json.load(open(Path(self.dir, "BO" + str(result_weights), "summary.json"), 'r'))

        final_summary["Highest Expected Reward"] = res.fun * -1
        final_summary["Optimal Rxn Utility and Penalty Params"] = res.x
        final_summary["Iterations of Bayesian Opt"] = self.bayes_iters
        final_summary["Solution summary"] = summary
        final_summary["Full Bayesian Opt Output"] = str(res)

        with open(self.dir/f'bayesian_summary.json','w') as f: 
                json.dump(final_summary, f, indent='\t')
        return self

    def expected_reward(self, weights):
        # full re-initialization of the problem
        all_weights = self.weights # make sure copy not reference
        all_weights[0] = weights[0] # rxn utility
        all_weights[1] = 0 # start material cost - clearing just in case
        all_weights[2] = weights[1] # rxn penalty
        print(all_weights)
        self.problem = self.initialize_problem()

        output_path = "BO" + str(all_weights) 
        Path(self.dir, output_path).mkdir(exist_ok=True) #TODO: what's wrong with self.output_dir?
        super().optimize(weights=all_weights, output_dir=Path(self.dir, output_path))
        summary = json.load(open(Path(self.dir, output_path, "summary.json"), 'r'))
        return -1 * summary['Expected Reward']
    
        #TODO: routes.json blank in bayesian directories [is this wrong, the routes in the solution directory WAS non-empty?]
    
    def optimize_weights(self):
        res = gp_minimize(self.expected_reward,    # the function to minimize
        [(0.0,1.0), (0.0,1.0)],      # the bounds on each dimension of x
        acq_func="EI",      # the acquisition function
        n_calls=self.bayes_iters,         # the number of evaluations of f
        n_initial_points=10,  # the number of random initialization points
        # random_state=[]   # the random seed
        )

        return res
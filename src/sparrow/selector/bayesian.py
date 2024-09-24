from typing import Dict, List
from pulp import lpSum
from skopt import gp_minimize

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
                 set_cycle_constraints: bool = False,
                 max_seconds: int = None,
                 extract_routes: bool = True, 
                 post_opt_class_score: str = None,
                 ) -> None:
        
        super().__init__(
            route_graph=route_graph, target_dict=target_dict, 
            rxn_scorer=rxn_scorer, condition_recommender=condition_recommender, 
            constrain_all_targets=constrain_all_targets, max_targets=max_targets, 
            coster=coster, weights=weights, output_dir=output_dir, 
            remove_dummy_rxns_first=remove_dummy_rxns_first, max_seconds=max_seconds, 
            clusters=clusters, N_per_cluster=N_per_cluster, rxn_classes=rxn_classes,
            rxn_classifier_dir=rxn_classifier_dir, max_rxn_classes=max_rxn_classes, 
            dont_buy_targets=dont_buy_targets, set_cycle_constraints=set_cycle_constraints,
            max_rxns=max_rxns, sm_budget=sm_budget, extract_routes=extract_routes,
            post_opt_class_score=post_opt_class_score,
            )
        self.solver = solver 
    
    def optimize(self):
        self.optimize_weights()
        return self

    # very basic reward metric
    def expected_reward_basic(self, weights):
        super().optimize(weights=weights)
        summary = self.extract_vars()
        return summary['Expected Reward']
    
    def optimize_weights(self):
        res = gp_minimize(-1 * self.expected_reward_basic,                  # the function to minimize
        [(0,1), (0,1)],      # the bounds on each dimension of x
        acq_func="EI",      # the acquisition function
        n_calls=15,         # the number of evaluations of f
        n_initial_points=10,  # the number of random initialization points
        random_state=[])   # the random seed

        return res
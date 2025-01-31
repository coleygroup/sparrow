from configargparse import ArgumentParser
from pathlib import Path

def get_args(args: str = None):
    parser = ArgumentParser()

    parser.add_argument('--config', is_config_file=True, help="the filepath of the configuration file")
    parser.add_argument('--target-csv', action='store', type=str, required=True, help="the filepath of the target csv file")
    parser.add_argument('--output-dir', action='store', type=str, default='sparrow_results')
    parser.add_argument('--no-routes', action='store_true', default=False)
    
    # path to graph object if just rerunning optimization 
    parser.add_argument('--graph', type=str, action='store', default=None, 
                        help='path to route graph json file. If provided, no route planning is performed')

    parser = add_tree_build_args(parser)
    parser = add_condition_rec_args(parser)
    parser = add_coster_args(parser)
    parser = add_scorer_args(parser)
    parser = add_cluster_args(parser)
    parser = add_rxn_class_args(parser)
    parser = add_optimization_args(parser)

    return parser.parse_args(args)

def add_tree_build_args(parser: ArgumentParser): 
    parser.add_argument('--path-finder', default=None, action='store', type=str,
                        choices=['lookup', 'api', 'apiv1'], 
                        help='type of tree builder to use')
    
    # lookup tree builder 
    parser.add_argument('--tree-lookup-dir', action='store', default=None, type=str,
                        help='path of lookup json file with combined retrosynthesis tree')

    # ASKCOS API tree builder  
    parser.add_argument('--time-per-target', default=30, action='store', type=int,
                        help='expansion time in seconds for each target')
    parser.add_argument('--max-ppg', default=10000, action='store', type=int, 
                        help='maximum price per gram in dollars for starting materials for ASKCOS MCTS tree search')
    parser.add_argument('--max-branching', default=20, action='store', type=int, 
                        help='maximum branch factor for ASKCOS MCTS tree search')
    parser.add_argument('--tree-host', default=None, action='store', type=str, 
                    help='host address for tree builder, if using ASKCOS API path finder')
    return parser

def add_condition_rec_args(parser: ArgumentParser):    
    parser.add_argument('--recommender', default=None, action='store', type=str,
                        choices=['lookup', 'api', 'apiv1'], 
                        help='type of context recommender to use')
    
    # API Recommender 
    parser.add_argument('--context-host', action='store', default=None, type=str,
                        help='host address for context recommender, if using API recommender')
    
    # Lookup Recommender
    parser.add_argument('--context-lookup', action='store', default=None, type=str,
                        help='path of lookup csv file for lookup context recommender (not implemented yet)')
    return parser

def add_coster_args(parser: ArgumentParser): 
    parser.add_argument('--coster', default=None, action='store', type=str,
                        choices=['naive', 'chemspace', 'lookup'], 
                        help='type of compound coster to use')
    parser.add_argument('--dont-buy-targets', action='store_true', default=False, 
                        help="ensures that the solution does not propose directly buying a target compound")
        
    # Naive Coster - no arguments 

    # Chemspace coster
    parser.add_argument('--key-path', action='store', type=str, default=str(Path.cwd()),
                        help='path that includes the file keys.py with chemspace api key')    
    
    # lookup coster 
    parser.add_argument('--inventory', action='store', type=str, default=None,
                        help='path of lookup file for lookup cost and buyability')
    parser.add_argument('--skip-canon', action='store_true', default=False,
                        help='whether to skip canonicalization of smiles in the inventory set')
    
    # optimization side of costing 
    parser.add_argument('--variable-costs', action='store_true', default=False,
                        help='whether to use a cost function of quantity instead of constant cost for buyables')
    
    return parser

def add_scorer_args(parser: ArgumentParser): 
    parser.add_argument('--scorer', default=None, action='store', type=str,
                        choices=['lookup', 'api', 'apiv1'], 
                        help='type of scorer to use')    
    
    # API Scorer 
    parser.add_argument('--scorer-host', action='store', default=None, type=str,
                        help='host address for reaction scorer, if using API recommender')
    
    # Lookup Recommender
    parser.add_argument('--scorer-lookup', action='store', default=None, type=str,
                        help='path of reaction scorer csv file for lookup reaction scorer (not implemented yet)')
    
    return parser 

def add_optimization_args(parser: ArgumentParser): 

    parser.add_argument('--formulation', action='store', default='linear', type=str, choices=['linear', 'expected_reward'],
                        help='whether to optimize the linear or expected reward formulation, latter requirs Gurobi license')

    # Specific to linear formulation 
    parser.add_argument('--reward-weight', action='store', type=float, default=1,
                        help='weighting factor for reward objective')
    parser.add_argument('--start-cost-weight', action='store', type=float, default=1,
                        help='weighting factor for starting material cost objective')
    parser.add_argument('--reaction-weight', action='store', type=float, default=1,
                        help='weighting factor for reaction objective')
    parser.add_argument('--solver', action='store', type=str, choices=['pulp', 'gurobi'],
                        default='pulp', help='solver to use for linear optimization')

    # specific to expected reward maximization
    parser.add_argument('--cost_of_rxn_weight', action='store', type=float, default=100,
                        help='weighting factor for reaction objective')
    parser.add_argument('--max-rxns', action='store', type=int, default=None,
                        help='maximum number of reaction steps to select')
    parser.add_argument('--starting-material-budget', action='store', type=float, default=None,
                        help='maximum budget on starting material costs (all on per g basis, does not consider amount needed!)')
    parser.add_argument('--prune-distance', '--prune', action='store', type=int, default=None,
                        help='To reduce the number of variables that is defined (and speed up solving), set this to a nonzero integer (~2X max route length)')

    # general settings and constraints 
    parser.add_argument('--constrain-all', action='store_true', default=False,
                        help='whether to constrain that all candidates are selection')
    parser.add_argument('--max-targets', action='store', default=None, type=int,
                        help='maximum number of selected targets (useful if testing is a bottleneck)')
    parser.add_argument('--acyclic', action='store_true', default=False, 
                        help='if the reaction network graph is known to be acyclic')
    parser.add_argument('--time-limit', action='store', default=12, type=float,
                        help='time limit on solving the optimization, in hours')
    
    # bayesian optimization 
    parser.add_argument('--bayes-iters', action='store', type=int, default=None,
                        help='By setting number of iterations, enable bayesian optimization of weights')
    
    return parser

def add_cluster_args(parser: ArgumentParser): 
    # currently specific to linear optimization 
    parser.add_argument('--diversity-weight', action='store', type=float, default=0,
                    help='weighting factor for diversity, encourages many clusters to be represented')

    # not specific to any formulation 
    parser.add_argument('--cluster', action='store', default=None, choices=[None, 'custom', 'similarity'],
                        help='How to define clusters. If "custom", should be included in targets.csv file ')
    parser.add_argument('--cluster-cutoff', action='store', type=float, default=0.7,
                        help='if using automatic clustering, cutoff for Butina clustering algorithm (lower cutoff -> more small clusters)') 
    parser.add_argument('--N-per-cluster', action='store', type=int, default=None,
                        help='To constrain that N compounds per cluster will be selected')
    return parser

def add_rxn_class_args(parser: ArgumentParser):
    parser.add_argument('--rxn-class-weight', action='store', type=float, default=0,
                        help='weighting factor for shared rxn classes')
    parser.add_argument('--rxn-classifier-path', action='store', type=str, default=None, 
                        help='Enter path to NameRXN directory or the full path to a csv with a custom rxn class mapping')
    parser.add_argument('--max-rxn-classes', action='store', type=int, default=None,
                        help='Constrain maximum number of reaction classes allowed to get to target set')
    return parser
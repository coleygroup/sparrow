from configargparse import ArgumentParser
from pathlib import Path

def get_args():
    parser = ArgumentParser()

    parser.add_argument('--config', is_config_file=True, help="the filepath of the configuration file")
    parser.add_argument('--target-csv', action='store', type=str, help="the filepath of the target csv file")
    parser.add_argument('--output-dir', action='store', type=str)

    # local askcos implementation (applies to local scorer and local context recommender)
    parser.add_argument('--askcos-path', action='store', default=None, type=str,
                    help='path where askcos exists if using askcos local modules')
    
    # path to graph object if just rerunning optimization 
    parser.add_argument('--graph', type=str, action='store', default=None, 
                        help='path to route graph json file. If provided, no route planning is performed')

    parser = add_tree_build_args(parser)
    parser = add_condition_rec_args(parser)
    parser = add_coster_args(parser)
    parser = add_scorer_args(parser)
    parser = add_optimization_args(parser)

    return parser.parse_args()

def add_tree_build_args(parser): 
    parser.add_argument('--path-finder', default=None, action='store', type=str,
                        choices=['lookup', 'api'], 
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

def add_condition_rec_args(parser):    
    parser.add_argument('--recommender', default=None, action='store', type=str,
                        choices=['lookup', 'local', 'api'], 
                        help='type of context recommender to use')
    
    # API Recommender 
    parser.add_argument('--context-host', action='store', default=None, type=str,
                        help='host address for context recommender, if using API recommender')

    # Local Recommender (--askcos-path in main function)
    
    # Lookup Recommender
    parser.add_argument('--context-lookup', action='store', default=None, type=str,
                        help='path of lookup csv file for lookup context recommender (not implemented yet)')
    return parser

def add_coster_args(parser): 
    parser.add_argument('--coster', default=None, action='store', type=str,
                        choices=['naive', 'chemspace'], 
                        help='type of compound coster to use')
    
    # Naive Coster - no arguments 

    # Chemspace coster
    parser.add_argument('--key-path', action='store', type=str, default=str(Path.cwd()),
                        help='path that includes the file keys.py with chemspace api key')    
    
    # lookup coster 
    parser.add_argument('--coster-lookup', action='store', type=str, default=None,
                        help='path of lookup file for lookup cost and buyability (not implemented yet)')
    
    return parser

def add_scorer_args(parser): 
    parser.add_argument('--scorer', default=None, action='store', type=str,
                        choices=['lookup', 'local', 'api'], 
                        help='type of scorer to use')    
    
    # API Scorer 
    parser.add_argument('--scorer-host', action='store', default=None, type=str,
                        help='host address for reaction scorer, if using API recommender')

    # Local Scorer (--askcos-path in main function)
    
    # Lookup Recommender
    parser.add_argument('--scorer-lookup', action='store', default=None, type=str,
                        help='path of reaction scorer csv file for lookup reaction scorer (not implemented yet)')
    
    return parser 

def add_optimization_args(parser): 

    parser.add_argument('--constrain-all', action='store_true', default=False,
                        help='whether to constrain that all candidates are selection (not implemented yet)')
    parser.add_argument('--reward-weight', action='store', type=float, default=1,
                        help='weighting factor for reward objective')
    parser.add_argument('--start-cost-weight', action='store', type=float, default=1,
                        help='weighting factor for starting material cost objective')
    parser.add_argument('--reaction-weight', action='store', type=float, default=1,
                        help='weighting factor for reaction objective')
    parser.add_argument('--solver', action='store', type=str, choices=['pulp, gurobi'],
                        default='pulp', help='solver to use for linear optimization')
    return parser


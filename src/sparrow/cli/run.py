from pathlib import Path
import pandas as pd 
import sys 
import numpy as np 

from sparrow.path_finder import AskcosAPIPlanner, LookupPlanner
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.condition_recommender import AskcosRecommender, AskcosAPIRecommender
from sparrow.scorer import AskcosScorer, AskcosAPIScorer
from sparrow.coster import ChemSpaceCoster, NaiveCoster, LookupCoster
from sparrow.cli.args import get_args

def get_target_dict(filepath: Path, clusters: bool = False): 
    df = pd.read_csv(filepath)
    target_dict = {
        smiles: reward
        for smiles, reward in zip(df['SMILES'], df['Reward'])
    }

    if not clusters: 
        return target_dict, list(df['SMILES']), None
    
    if 'Cluster' not in df: 
        print(f'No "Cluster" column in {filepath}, but custom cluster argument was specified')
        print('Continuing without clusters')
        return target_dict, list(df['SMILES']), None

    c_names = set([c for c in df['Cluster'] if not np.isnan(c)])
    cluster_smis = {
        c: list(df.loc[df.Cluster==c]['SMILES'])
        for c in c_names
    }
    return target_dict, list(df['SMILES']), cluster_smis

def optimize(selector, params):

    selector.define_variables()
    selector.set_objective()
    selector.set_constraints(set_cycle_constraints=not params['acyclic'])
    
    # solver = {'pulp': None, 'gurobi': 'GUROBI'}[params['solver']]
    selector.optimize() # solver='GUROBI' for GUROBI (license needed)
    # selector.optimize_MO()
    return selector 


def get_path_storage(params, targets): 
    if params['graph'] is not None: 
        return None 
    
    if params['path_finder'] == 'lookup': 
        planner = LookupPlanner(
            json_dir=Path(params['tree_lookup_dir']),
            output_dir=Path(params['output_dir']),
        )
    elif params['path_finder'] == 'api': 
        planner = AskcosAPIPlanner( 
            host=params['tree_host'],
            output_dir=Path(params['output_dir']),
            time_per_target=params['time_per_target'], 
            max_ppg=params['max_ppg'],
            max_branching=params['max_branching'],
        )

    return planner.get_save_trees(targets)

def build_recommender(params): 
    rec = params['recommender']
    if rec == 'api': 
        return AskcosAPIRecommender(host=params['context_host'])
    elif rec == 'local': 
        return AskcosRecommender(askcos_path=params['askcos_path'])
    elif rec is None: 
        return None    
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Context recommender {rec} not implemented')

def build_scorer(params): 
    rec = params['scorer']
    if rec == 'api': 
        return AskcosAPIScorer(host=params['scorer_host'])
    elif rec == 'local': 
        return AskcosScorer(askcos_path=params['askcos_path'])
    elif rec is None: 
        return None  
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Scorer {rec} not implemented')
    
def build_coster(params): 
    rec = params['coster']
    if rec == 'chemspace': 
        sys.path.append(params['key_path'])
        from keys import chemspace_api_key
        return ChemSpaceCoster(api_key=chemspace_api_key)
    elif rec == 'naive': 
        return NaiveCoster()
    elif rec == 'lookup': 
        return LookupCoster(lookup_file=params['inventory'], canonicalize= not params['skip_canon'])
    elif rec is None: 
        return None  
    else:
        raise NotImplementedError(f'Scorer {rec} not implemented')
     
def build_selector(params, target_dict, storage_path, clusters):
    if storage_path is None: 
        graph = RouteGraph(node_filename=params['graph'])
    else: 
        graph = RouteGraph(node_filename=storage_path)

    # weights = [params['reward_weight'], params['start_cost_weight'], params['reaction_weight'], params['diversity_weight']]

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=build_recommender(params), 
        output_dir=Path(params['output_dir']),
        rxn_scorer=build_scorer(params),
        coster=build_coster(params),
        cost_per_rxn=params['cost_of_rxn_weight'],
        constrain_all_targets=params['constrain_all'],
        max_targets=params['max_targets'],
        custom_clusters=clusters,
        max_rxns=params['max_rxns'],
        sm_budget=params['starting_material_budget'],
        variable_costs=params['variable_costs'],
    )

    if storage_path is not None: 
        filepath = Path(params['output_dir'])/'trees_w_info.json'
        print(f'Saving route graph with contexts, reaction scores, and costs to {filepath}')
        selector.graph.to_json(filepath)
    
    return selector

def save_args(params): 
    filename = Path(params['output_dir'])/'params.ini'
    with open(filename, 'w') as f:  
        for k, v in sorted(params.items()):
            f.write(f'{k}: {v}\n')

    return 

def run():
    args = get_args()
    params = vars(args)

    print('SPARROW will be run with the following parameters:')    
    for k, v in sorted(params.items()):
        if v is not None: 
            print(f"  {k}: {v}")
    print(flush=True)
    output_dir = Path(params['output_dir'])
    output_dir.mkdir(exist_ok=True, parents=True)
    
    save_args(params)

    target_dict, targets, clusters = get_target_dict(params['target_csv'], clusters=params['custom_cluster']) 
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path, clusters)
    selector = optimize(selector, params)

if __name__ == '__main__':
    run()
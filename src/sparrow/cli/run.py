from pathlib import Path
import pandas as pd 
import sys 
import numpy as np 
import json 

from sparrow.path_finder import AskcosAPIPlanner, LookupPlanner
from sparrow.route_graph import RouteGraph
from sparrow.selector.linear import LinearSelector
from sparrow.selector.nonlinear import ExpectedRewardSelector, PrunedERSelector
from sparrow.condition_recommender import AskcosAPIRecommender
from sparrow.scorer import AskcosAPIScorer
from sparrow.coster import ChemSpaceCoster, NaiveCoster, LookupCoster
from sparrow.cli.args import get_args
from sparrow.utils import cluster_utils

def get_target_dict(filepath: Path): 
    df = pd.read_csv(filepath)
    target_dict = {
        smiles: reward
        for smiles, reward in zip(df['SMILES'], df['Reward'])
    }
    return target_dict, list(df['SMILES'])

def get_clusters(cluster_type, filepath, cutoff, outdir=None): 
    if cluster_type is None: 
        return None 
    elif cluster_type == 'custom':
        df = pd.read_csv(filepath)
        if 'Cluster' not in df: 
            print(f'No "Cluster" column in {filepath}, treating each additional column as a cluster group')
            c_names = set(df.columns)
            c_names.remove('SMILES')
            c_names.remove('Reward')
            cluster_smis = {
                c_name: list(df.loc[df[c_name]==True]['SMILES'])
                for c_name in c_names
            }
            return cluster_smis

        c_names = set([c for c in df['Cluster'] if not np.isnan(c)])
        cluster_smis = {
            c: list(df.loc[df.Cluster==c]['SMILES'])
            for c in c_names
        }
        return cluster_smis
    elif cluster_type == 'similarity': 
        df = pd.read_csv(filepath)
        smis = list(df['SMILES'])
        clusters = cluster_utils.cluster_smiles(smis, cutoff=cutoff)
        cluster_smis = {
            f'SimCluster_{i}': [smis[j] for j in clusters[i]]
            for i in range(len(clusters))
        }
        
        cs_file = Path(outdir) / 'clusters.json'
        print(f'Saving list of {len(cluster_smis)} clusters to {cs_file}')
        
        with open(cs_file,'w') as f: 
            json.dump(cluster_smis, f, indent='\t')

        return cluster_smis

def optimize(selector, params):

    selector.define_variables()
    selector.set_objective()
    selector.set_constraints(set_cycle_constraints=not params['acyclic'])
    selector.optimize() 
    selector.extract_vars()
    
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
    if params['formulation'] == 'expected_reward' and params['prune_distance'] is None:
        selector = ExpectedRewardSelector(
            route_graph=graph,
            target_dict=target_dict,
            rxn_scorer=build_scorer(params),
            condition_recommender=build_recommender(params),
            constrain_all_targets=params['constrain_all'],
            max_targets=params['max_targets'],
            coster=build_coster(params),
            cost_per_rxn=params['cost_of_rxn_weight'],
            output_dir=Path(params['output_dir']),
            clusters=clusters,
            max_rxns=params['max_rxns'],
            sm_budget=params['starting_material_budget'],
            dont_buy_targets=params['dont_buy_targets']
        )
    elif params['formulation'] == 'expected_reward': 
        selector = PrunedERSelector(
            route_graph=graph,
            target_dict=target_dict,
            rxn_scorer=build_scorer(params),
            condition_recommender=build_recommender(params),
            constrain_all_targets=params['constrain_all'],
            max_targets=params['max_targets'],
            coster=build_coster(params),
            cost_per_rxn=params['cost_of_rxn_weight'],
            output_dir=Path(params['output_dir']),
            clusters=clusters,
            max_rxns=params['max_rxns'],
            sm_budget=params['starting_material_budget'],
            dont_buy_targets=params['dont_buy_targets'],
            N_per_cluster=params['N_per_cluster']
        )
    elif params['formulation'] == 'expected_reward': 
        selector = PrunedERSelector(
            route_graph=graph,
            target_dict=target_dict,
            rxn_scorer=build_scorer(params),
            condition_recommender=build_recommender(params),
            constrain_all_targets=params['constrain_all'],
            max_targets=params['max_targets'],
            coster=build_coster(params),
            cost_per_rxn=params['cost_of_rxn_weight'],
            output_dir=Path(params['output_dir']),
            cluster_cutoff=params['cluster_cutoff'],
            custom_clusters=clusters,
            max_rxns=params['max_rxns'],
            sm_budget=params['starting_material_budget'],
            dont_buy_targets=params['dont_buy_targets'],
            prune_distance=params['prune_distance']
        )        
    else: 
        weights = [params['reward_weight'], params['start_cost_weight'], params['reaction_weight'], params['diversity_weight']]
        
        selector = LinearSelector(
            route_graph=graph,
            target_dict=target_dict,
            rxn_scorer=build_scorer(params),
            condition_recommender=build_recommender(params),
            constrain_all_targets=params['constrain_all'],
            max_targets=params['max_targets'],
            coster=build_coster(params),
            weights=weights,
            output_dir=Path(params['output_dir']),
            cluster_cutoff=params['cluster_cutoff'],
            custom_clusters=clusters,
            dont_buy_targets=params['dont_buy_targets'],
            solver=params['solver'],
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

    target_dict, targets = get_target_dict(params['target_csv']) 
    if params['diversity_weight'] > 0: 
        print(f'Overwriting parameter "--cluster" to be similarity for diversity objective')
        params['cluster'] = 'similarity'
    clusters = get_clusters(
        cluster_type=params['cluster'], 
        filepath=params['target_csv'], 
        cutoff=params['cluster_cutoff'], 
        outdir=params['output_dir']
    )
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path, clusters)
    selector = optimize(selector, params)

if __name__ == '__main__':
    run()
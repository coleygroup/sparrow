from pathlib import Path
import pandas as pd 
import sys 
import numpy as np 
import json 
from tqdm import tqdm

from sparrow.path_finder import AskcosAPIPlanner, LookupPlanner, AskcosV1APIPlanner
from sparrow.route_graph import RouteGraph
from sparrow.selector.linear import LinearSelector
from sparrow.selector.bayesian import BOLinearSelector
from sparrow.selector.nonlinear import ExpectedRewardSelector, PrunedERSelector
from sparrow.condition_recommender import AskcosAPIRecommender, AskcosV1APIRecommender
from sparrow.scorer import AskcosAPIScorer, AskcosV1APIScorer
from sparrow.coster import ChemSpaceCoster, NaiveCoster, LookupCoster
from sparrow.rxn_classifier import NameRxnClass, LookupClass
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

        c_names = set([c for c in df['Cluster'] if not pd.isnull(c)])
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
    elif params['path_finder'] == 'apiv1':
        planner = AskcosV1APIPlanner( 
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
    elif rec == 'apiv1': 
        return AskcosV1APIRecommender(host=params['context_host'])
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
    elif rec == 'apiv1': 
        return AskcosV1APIScorer(host=params['scorer_host'])
    elif rec is None: 
        return None  
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Scorer {rec} not implemented')
    
def build_coster(params): 
    rec = params['coster']
    if rec == 'chemspace': 
        sys.path.append(str(Path(params['key_path']).parent))
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
    
def build_rxn_classes(params, graph: RouteGraph):
    cls_path = params['rxn_classifier_path']
    if cls_path == None:
        return None
    elif not Path(cls_path).is_dir(): 
        classifier = LookupClass(cls_path)
    else:
        classifier = NameRxnClass(cls_path)

    rxn_smis = [r.smiles for r in graph.non_dummy_nodes()]
    size = np.ceil(len(rxn_smis) / 5).astype(int) # currently splits into 5 batches (hard-coded)
    batch_rxns = list(
        map(lambda x: rxn_smis[x * size:x * size + size],
        list(range(5)))
    )
    batch_classes = [
        classifier.get_rxn_classes(batch) 
        for batch in tqdm(batch_rxns, desc='Classifying reactions')
    ]

    class_nums = [c for batch in batch_classes for c in batch]

    rxn_classes = {}
    for c, rxn in zip(class_nums, rxn_smis): 
        if c in rxn_classes.keys():
            rxn_classes[c].append(rxn)
        else:
            rxn_classes[c] = [rxn]

    if isinstance(classifier, NameRxnClass):
        rxn_df = pd.DataFrame({'SMILES': rxn_smis, 'Class': class_nums})
        rxn_df.to_csv(Path(params['output_dir'])/'reaction_classes.csv', index=False)
        
    return rxn_classes
     
def build_selector(params, target_dict, storage_path, clusters):
    # storage_path is a dict of SMILES to reward
    if storage_path is None: 
        graph = RouteGraph(node_filename=params['graph'])
    else: 
        graph = RouteGraph(node_filename=storage_path)

    weights = [params['reward_weight'], params['start_cost_weight'], params['reaction_weight'], params['diversity_weight'], params['rxn_class_weight']]
    args = {
        'route_graph': graph,
        'target_dict': target_dict,
        'rxn_scorer': build_scorer(params),
        'condition_recommender': build_recommender(params),
        'constrain_all_targets': params['constrain_all'],
        'max_targets': params['max_targets'],
        'coster': build_coster(params),
        'output_dir': Path(params['output_dir']),
        'clusters': clusters,
        'max_rxns': params['max_rxns'],
        'sm_budget': params['starting_material_budget'], 
        'dont_buy_targets': params['dont_buy_targets'],
        'N_per_cluster': params['N_per_cluster'],
        'cycle_constraints': not params['acyclic'],
        'max_seconds': params['time_limit']*3600,
        'rxn_classes': build_rxn_classes(params, graph) if 'rxn_classifier_path' in params else None,
        'max_rxn_classes': params['max_rxn_classes'] if 'max_rxn_classes' in params else None,
    }
    
    if params['formulation'] == 'expected_reward' and params['prune_distance'] is None:
        selector = ExpectedRewardSelector(
            cost_per_rxn=params['cost_of_rxn_weight'],
            **args
        )
    elif params['formulation'] == 'expected_reward': 
        selector = PrunedERSelector(
            cost_per_rxn=params['cost_of_rxn_weight'],
            prune_distance=params['prune_distance'],
            **args
        )      
    elif params['bayes_iters']:
        selector = BOLinearSelector(
            weights=weights,
            rxn_classifier_dir = params['rxn_classifier_path'] if 'rxn_classifier_path' in params else None,
            solver=params['solver'],
            bayes_iters=params['bayes_iters'],
            **args
        ) 
    else: 
        selector = LinearSelector(
            weights=weights,
            rxn_classifier_dir = params['rxn_classifier_path'] if 'rxn_classifier_path' in params else None,
            solver=params['solver'],
            **args
        )

    # if no graph is given, creates graph from the info file path
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
    if params['diversity_weight'] > 0 and params['cluster'] is None: 
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
    selector = selector.formulate_and_optimize(extract_routes=not params['no_routes'])

if __name__ == '__main__':
    run()
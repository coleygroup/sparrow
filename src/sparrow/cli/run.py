from pathlib import Path
import pandas as pd 
import json 
import sys 
import numpy as np 
from tqdm import tqdm 

from sparrow.path_finder import AskcosAPIPlanner, LookupPlanner
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.condition_recommender import AskcosAPIRecommender
from sparrow.scorer import AskcosAPIScorer
from sparrow.coster import ChemSpaceCoster, NaiveCoster, LookupCoster
from sparrow.cli.args import get_args



def get_target_dict(filepath: Path): 
    df = pd.read_csv(filepath)
    target_dict = {
        smiles: reward
        for smiles, reward in zip(df['SMILES'], df['Reward'])
    }
    return target_dict, list(df['SMILES'])

def optimize(selector, params):

    selector.define_variables()
    selector.set_objective()
    selector.set_constraints(set_cycle_constraints=not params['acyclic'])
    
    solver = {'pulp': None, 'gurobi': 'GUROBI'}[params['solver']]
    selector.optimize(solver=solver) # solver='GUROBI' for GUROBI (license needed)

    return selector 

def export_selected_nodes(selector: RouteSelector, rxn_list, starting_list, target_list, output_dir): 
    storage = {'Starting Materials': [], 'Reactions': [], 'Targets': []}
    graph = selector.graph

    for rxn_id in rxn_list:
        node = graph.node_from_id(rxn_id)
        if node.dummy: 
            storage['Reactions'].append({
                'smiles': node.smiles,
                'starting material cost ($/g)': selector.cost_of_dummy(node), 
            })
        else: 
            storage['Reactions'].append({
                'smiles': node.smiles,
                'conditions': node.get_condition(1)[0], 
                'score': node.score,
            })

    for cpd_id in starting_list: 
        node = graph.node_from_id(cpd_id)
        storage['Starting Materials'].append({
            'smiles': node.smiles,
            'cost': node.cost_per_g
        })

    for cpd_id in target_list: 
        node = graph.node_from_id(cpd_id)
        storage['Targets'].append({
            'smiles': node.smiles,
            'reward': node.reward,
        })
    
    with open(output_dir/'solution_list_format.json','w') as f:
        json.dump(storage, f, indent='\t')

    return storage 

def extract_vars(selector: RouteSelector, output_dir, extract_routes=True): 

    nonzero_vars = [
            var for var in selector.problem.variables() if var.varValue > 0.01
        ]
    rxn_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('rxn')]
    mol_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('mol')]
    dummy_ids = [rxn for rxn in rxn_ids if selector.graph.node_from_id(rxn).dummy]
    non_dummy_ids = [rxn for rxn in rxn_ids if selector.graph.node_from_id(rxn).dummy == 0]

    selected_targets = set(mol_ids) & set(selector.targets)
    selected_starting = set([selector.graph.child_of_dummy(dummy) for dummy in dummy_ids])
    print(f'{len(selected_targets)} targets selected using {len(non_dummy_ids)} reactions and {len(selected_starting)} starting materials')
    export_selected_nodes(selector, rxn_ids, selected_starting, selected_targets, output_dir)
    
    avg_rxn_score = np.mean([selector.graph.node_from_id(rxn).score for rxn in non_dummy_ids]) if len(non_dummy_ids) > 0 else None

    summary = {
        'Weights': selector.weights,
        'Number targets': len(selected_targets), 
        'Fraction targets': len(selected_targets)/len(selector.targets),
        'Total reward': sum([selector.target_dict[tar] for tar in selected_targets]),
        'Possible reward': sum(selector.target_dict.values()),
        'Number starting materials': len(selected_starting),
        'Cost starting materials': sum([selector.cost_of_dummy(dummy_id=d_id) for d_id in dummy_ids]),
        'Number reaction steps': len(non_dummy_ids),
        'Average reaction score': avg_rxn_score,
    }

    if extract_routes:
        storage = {}
        for target in tqdm(selected_targets, desc='Extracting routes'): 
            store_dict = {'Compounds':[], 'Reactions':[]}
            smi = selector.graph.smiles_from_id(target)
            storage[smi] = find_mol_parents(store_dict, target, mol_ids, rxn_ids, selector)
            storage[smi]['Reward'] = selector.target_dict[target]

        with open(output_dir/f'routes.json','w') as f: 
            json.dump(storage, f, indent='\t')

    return summary 

def find_rxn_parents(store_dict, rxn_id, selected_mols, selected_rxns, selector: RouteSelector):
    graph = selector.graph 
    par_ids = [n.id for n in graph.node_from_id(rxn_id).parents.values()]
    selected_pars = set(par_ids) & set(selected_mols)
    for par in selected_pars: 
        store_dict['Compounds'].append(graph.smiles_from_id(par))
        store_dict = find_mol_parents(store_dict, par, selected_mols, selected_rxns, selector)
    return store_dict

def find_mol_parents(store_dict, mol_id, selected_mols, selected_rxns, selector: RouteSelector): 
    graph = selector.graph 
    par_ids = [n.id for n in graph.node_from_id(mol_id).parents.values()]
    selected_pars = set(par_ids) & set(selected_rxns)
    for par in selected_pars: 
        node = graph.node_from_id(par)
        if node.dummy: 
            store_dict['Reactions'].append({
                'smiles': node.smiles,
                'starting material cost ($/g)': selector.cost_of_dummy(node), 
            })
        else: 
            store_dict['Reactions'].append({
                'smiles': node.smiles,
                'conditions': node.get_condition(1)[0], 
                'score': node.score,
            })
        store_dict = find_rxn_parents(store_dict, par, selected_mols, selected_rxns, selector)
    return store_dict

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
     
def build_selector(params, target_dict, storage_path):
    if storage_path is None: 
        graph = RouteGraph(node_filename=params['graph'])
    else: 
        graph = RouteGraph(node_filename=storage_path)

    weights = [params['reward_weight'], params['start_cost_weight'], params['reaction_weight'], params['diversity_weight']]

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=build_recommender(params), 
        output_dir=Path(params['output_dir']),
        rxn_scorer=build_scorer(params),
        coster=build_coster(params),
        weights=weights,
        constrain_all_targets=params['constrain_all'],
        max_targets=params['max_targets'],
        dont_buy_targets=params['dont_buy_targets'],
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
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path)
    selector = optimize(selector, params)
    summary = extract_vars(selector, output_dir, extract_routes=not params['no_routes'] )
    
    with open(output_dir/'summary.json', 'w') as f:
        json.dump(summary, f, indent='\t')

if __name__ == '__main__':
    run()
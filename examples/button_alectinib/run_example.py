import os, sys
sys.path.append('/home/jfromer/sparrow/askcos-base/askcos-core') # change this for your system to where askcos folder is located
sys.path.append('/home/jfromer/sparrow/')

from pathlib import Path
import pandas as pd 
import json 
from typing import List, Dict
from tqdm import tqdm
from keys import chemspace_api_key

from sparrow.tree_build_utils import get_paths
from sparrow.json_utils import storage_from_api_response, save_storage_dict
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.condition_recommender import AskcosRecommender
from sparrow.scorer import AskcosScorer
from sparrow.coster import ChemSpaceCoster

def get_target_dict(filepath: Path): 
    df = pd.read_csv(filepath)
    target_dict = {
        smiles: reward
        for smiles, reward in zip(df['SMILES'], df['Reward'])
    }
    return target_dict, list(df['SMILES'])

def get_trees(targets, base_dir, time_per_target=60, max_branching=str(20), host='https://18.4.94.12'):
    params = {
        'buyable_logic': 'or',
        'max_depth': '10',
        'expansion_time': str(time_per_target),
        'max_ppg': '100',
        'return_first': 'false', 
        'max_branching': max_branching,
    }

    results = get_paths(targets, host=host, store_dir=base_dir/'paths_cariprazine', params=params)

    return results 

def combine_outputs(result_ls: List[Path], path_dir: Path):

    trees = {}
    for p in tqdm(result_ls, desc='Combining ASKCOS outputs'): 
        with open(p,'r') as f: 
            entry = json.load(f)
        for smi, path in entry.items(): 
            if 'output' in path and len(path['output'])>0: 
                trees[smi] = {"result": path} 
    
    with open(path_dir/'askcos_outputs.json','w') as f: 
        json.dump(trees, f, indent='\t')

    return trees

def process_tree(trees: Dict = None, tree_file: Path = None) -> Dict: 
    if trees is None: 
        with open(tree_file, 'r') as f: 
            trees = json.load(f)
    

    storage = None
    for response in tqdm(trees.values(), desc='Reading ASKCOS API results'):
        storage = storage_from_api_response(response, storage)

    return storage 

def add_contexts_scores_to_trees(json_file: str, target_dict: str, base_dir) -> None: 
    graph = RouteGraph(node_filename=json_file)

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=AskcosRecommender(), 
        output_dir=base_dir,
        rxn_scorer=AskcosScorer(),
    )

    selector.graph.to_json(base_dir/'trees_w_context_scores.json')

    return selector


def cost_compounds(json_file: str, target_dict: Dict, base_dir):
    graph = RouteGraph(node_filename=json_file)

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=None, # already in graph   
        rxn_scorer=None,     # already in graph      
        coster=ChemSpaceCoster(api_key=chemspace_api_key),
        output_dir=base_dir,
    )
    selector.graph.to_json(base_dir/'trees_w_context_scores_costs.json')
    return selector

def optimize(json_file: str, target_dict: Dict, base_dir, weights=[1,1,1,1]):
    graph = RouteGraph(node_filename=json_file)
    
    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=None, # already in graph
        rxn_scorer=None, # already in graph       
        coster=None, # already in graph 
        output_dir=base_dir,
        remove_dummy_rxns_first=True, 
        weights=weights
    )
    selector.define_variables()
    selector.set_objective()
    selector.set_constraints()
    selector.optimize(solver=None) # solver='GUROBI' for GUROBI (license needed)

    return selector 

def export_selected_nodes(selector: RouteSelector, rxn_list, starting_list, target_list, base_dir): 
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
    
    with open(base_dir/'optimal_routes.json','w') as f:
        json.dump(storage, f, indent='\t')

    return storage 



def extract_vars(selector: RouteSelector, base_dir): 
    sys.setrecursionlimit(100000)
    nonzero_vars = [
            var for var in selector.problem.variables() if var.varValue > 0.01
        ]
    rxn_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('rxn')]
    mol_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('mol')]

    selected_targets = set(mol_ids) & set(selector.targets)
    starting_mats = set([node.id for node in selector.graph.buyable_nodes()])
    selected_starting = set(mol_ids) & starting_mats
    print(f'{len(selected_targets)} targets selected using {len(rxn_ids)} reactions and {len(selected_starting)} starting materials')
    export_selected_nodes(selector, rxn_ids, selected_starting, selected_targets, base_dir)

    storage = {}
    for target in selected_targets: 
        store_dict = {'Compounds':[], 'Reactions':[]}
        smi = selector.graph.smiles_from_id(target)
        storage[smi] = find_mol_parents(store_dict, target, mol_ids, rxn_ids, selector.graph)
        storage[smi]['Reward'] = selector.target_dict[target]

    with open(base_dir/f'routes_{len(selected_targets)}tars.json','w') as f: 
        json.dump(storage, f, indent='\t')

    return storage 

def find_rxn_parents(store_dict, rxn_id, selected_mols, selected_rxns, graph): 
    par_ids = [n.id for n in graph.node_from_id(rxn_id).parents.values()]
    selected_pars = set(par_ids) & set(selected_mols)
    for par in selected_pars: 
        store_dict['Compounds'].append(graph.smiles_from_id(par))
        store_dict = find_mol_parents(store_dict, par, selected_mols, selected_rxns, graph)
    return store_dict

def find_mol_parents(store_dict, mol_id, selected_mols, selected_rxns, graph): 
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
        store_dict = find_rxn_parents(store_dict, par, selected_mols, selected_rxns, graph)
    return store_dict

if __name__ == '__main__':
    drug = 'cariprazine'
    base_dir = Path('examples/button_analogs')
    host = 'http://molgpu01.mit.edu:9100' # 'https://3.139.77.247/'

    path_dir = base_dir/f'paths_{drug}'
    
    target_dict, targets = get_target_dict(base_dir/f'{drug}_rewards.csv') 
    # get_trees(targets, base_dir, time_per_target=60, host=host)

    paths_ls = list(path_dir.glob('paths*'))
    # trees = combine_outputs(paths_ls, path_dir)
    
    # tree_file = list(path_dir.glob('*outputs*'))[0]
    #storage = process_tree(tree_file)

    save_storage_dict(storage, base_dir/'trees.json')
    selector = add_contexts_scores_to_trees(base_dir/'trees.json', target_dict, base_dir)
    selector = cost_compounds(base_dir/'trees_w_context_scores.json', target_dict, base_dir)
    selector = optimize(base_dir/'trees_w_context_scores_costs.json', target_dict, base_dir, weights=[5,1,1,1])
    storage = extract_vars(selector, base_dir)

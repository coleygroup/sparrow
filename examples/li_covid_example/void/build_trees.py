import os, sys
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located

import json 
import pandas as pd 
from typing import Dict
from tqdm import tqdm
from pathlib import Path
from sparrow.json_utils import storage_from_api_response, save_storage_dict 
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.scorer import AskcosScorer
from sparrow.condition_recommender import AskcosLookupRecommender
from sparrow.coster import ChemSpaceCoster
from sparrow.visualizer import Visualizer
from sparrow.json_utils import build_retro_graph_local, build_retro_graph_api

def make_target_dict(input_csv: str) -> Dict: 
    df = pd.read_csv(input_csv)
    target_dict = {
        smiles: reward_from_rank(rank)
        for smiles, rank in zip(df['SMILES'], df['Reward'])
    }
    return target_dict
    
def reward_from_rank(rank: int) -> float: 
    return rank

def process_trees(tree_dir: str) -> Dict: 
    storage = None 
    for tree_file in tqdm(sorted(Path(tree_dir).glob('*.json')), desc='Reading tree files'):
        with open(tree_file, 'r') as f:
            trees = json.load(f)

        for response in trees.values():     
            storage = storage_from_api_response(response, storage)

    return storage 

def add_contexts_scores_to_trees(json_file: str, target_dict: str, ) -> None: 
    import askcos.utilities.contexts as context_cleaner
    graph = RouteGraph(node_filename=json_file)

    recommender = AskcosLookupRecommender(
        lookup_file='darpa_example/contexts/', 
        context_cleaner=context_cleaner,
    )
    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=recommender, 
        output_dir='darpa_example',
        rxn_scorer=AskcosScorer(),
    )

    selector.graph.to_json('darpa_example/trees_w_context_scores.json')

    return selector

def cost_compounds(json_file: str, target_dict: Dict):
    graph = RouteGraph(node_filename=json_file)

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=None, # already in graph   
        rxn_scorer=None,     # already in graph      
        coster=ChemSpaceCoster(api_key="dZH7vZYK2JDKWxgMSCKIBQZcKfteL395UuYtCuHoVk1WUcpq1MIeiPn95mBLsXOh"),
        output_dir='darpa_example',
    )
    selector.graph.to_json('darpa_example/trees_w_context_scores_costs.json')
    return 

def optimize(json_file: str, target_dict: Dict):
    graph = RouteGraph(node_filename=json_file)
    
    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=None, # already in graph
        rxn_scorer=None, # already in graph       
        coster=None, # already in graph 
        output_dir='darpa_example',
        remove_dummy_rxns_first=True, 
        weights=[1,1,1]
    )
    selector.define_variables()
    selector.set_objective()
    selector.set_constraints()
    selector.optimize(solver=None) # solver='GUROBI' for GUROBI (license needed)

    return selector 

def extract_vars(selector: RouteSelector): 
    nonzero_vars = [
            var for var in selector.problem.variables() if var.varValue > 0.01
        ]
    rxn_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('rxn')]
    mol_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('mol')]

    selected_targets = set(mol_ids) & set(selector.targets)
    starting_mats = set([node.id for node in selector.graph.buyable_nodes()])
    selected_starting = set(mol_ids) & starting_mats
    print(f'{len(selected_targets)} targets selected using {len(rxn_ids)} reactions and {len(selected_starting)} starting materials')

    storage = {}
    for target in selected_targets: 
        store_dict = {'Compounds':[], 'Reactions':[]}
        smi = selector.graph.smiles_from_id(target)
        storage[smi] = find_mol_parents(store_dict, target, mol_ids, rxn_ids, selector.graph)
        storage[smi]['Reward'] = selector.target_dict[target]

    with open(f'darpa_example/routes_{len(selected_targets)}tars.json','w') as f: 
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

if __name__=='__main__': 
    base_dir = Path('examples/li_covid_example')
    target_dict = make_target_dict(base_dir/'li_rewards_copy.csv')
    
    for tar in target_dict.keys(): 
        print(tar)
        storage = build_retro_graph_local(        
            target_smis=[tar],
            filename=base_dir/'trees.json',        
            time_per_target=30,
        )
        # storage = build_retro_graph_api(        
        #     target_smis=[tar],
        #     host='https://18.4.94.12',
        #     filename=base_dir/'trees.json',       
        #     time_per_target=30,
        # )
        # save_storage_dict(storage, base_dir/'askcos_trees.json')
    # selector = add_contexts_scores_to_trees('darpa_example/trees.json', target_dict)v
    # selector = optimize('darpa_example/trees_w_context_scores_costs.json', target_dict)
    # storage = extract_vars(selector)
    # print('done')

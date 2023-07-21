import os, sys
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located

import json 
import pandas as pd 
from typing import Dict
from tqdm import tqdm
from pathlib import Path
from sparrow.tree_build_utils import storage_from_api_response, save_storage_dict 
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.scorer import AskcosScorer
from sparrow.condition_recommender import AskcosLookupRecommender
import askcos.utilities.contexts as context_cleaner

def make_target_dict(input_csv: str) -> Dict: 
    df = pd.read_csv(input_csv)
    target_dict = {
        smiles: reward_from_rank(rank)
        for smiles, rank in zip(df['SMILES'], df['NDS Rank (0 is best)'])
    }
    return target_dict
    
def reward_from_rank(rank: int) -> float: 
    return 20 - rank 

def process_trees(tree_dir: str) -> Dict: 
    storage = None 
    for tree_file in tqdm(sorted(Path(tree_dir).glob('*.json')), desc='Reading tree files'):
        with open(tree_file, 'r') as f:
            trees = json.load(f)

        for response in trees.values():     
            storage = storage_from_api_response(response, storage)

    return storage 

def add_contexts_to_trees(json_file: str, target_dict: str, ) -> None: 
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
    return selector

if __name__=='__main__': 
    target_dict = make_target_dict('darpa_example/test_mols.csv')
    # storage = process_trees('darpa_example/trees/')
    # save_storage_dict(storage, 'darpa_example/trees.json')
    selector = add_contexts_to_trees('darpa_example/trees.json', target_dict)
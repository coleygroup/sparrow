import os
os.environ['KERAS_BACKEND'] = 'theano'

from rdkit import Chem

from typing import List, Dict, Union
import requests 
import json
import numpy as np 
from pathlib import Path
from tqdm import tqdm 


def build_retro_graph_local(
        target_smis: List[str],
        filename: Union[str, Path] = 'debug/askcos_paths.json',        
        n_cpus: int = 5,
        time_per_target: int = 15,
        storage = None, 
        ) -> Dict: 

    """ 
    Builds retrosynthesis trees for the targets 
    and outputs a json file with path information 
    to the specified filename. Returns the storage dictionary 
    that was saved to the json file. Assumes a local 
    version of ASKCOS exists and is in the path   
    """    
    import tensorflow as tf 
    from askcos.retrosynthetic.mcts.tree_builder import MCTS
    from askcos.utilities.io.logger import MyLogger

    tf.config.set_visible_devices([], 'GPU')
    MyLogger.initialize_logFile()
    celery = False 

    
    mols = [Chem.MolFromSmiles(smi) for smi in target_smis]
    
    is_first_target = True
    
    Tree = MCTS(nproc=n_cpus)

    for mol in mols:
        
        soft_reset = False if is_first_target else True 

        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)   
        print('expanding target {}'.format(smiles)) 
        paths, _, _ = Tree.get_buyable_paths(
            smiles,
            nproc=n_cpus,
            expansion_time=time_per_target, 
            termination_logic={'or': ['buyable']},
            # min_chemical_history_dict={'as_reactant':1, 'as_product':1,'logic':'and'},
            soft_reset=soft_reset,
            soft_stop=True
        )

        is_first_target = False

        print('done for target {}'.format(smiles))

        storage = storage_from_paths(paths, storage)

    Tree.stop()

    with open(filename, 'w') as f:
        json.dump(make_dict_jsonable(storage), f, indent="\t")
    
    print('saved storage to {}'.format(filename))
    
    for p in Tree.workers:
        if p and p.is_alive():
            p.terminate()

    return storage

def build_retro_graph_api(        
        target_smis: List[str],
        host: str,
        filename: Union[str, Path] = 'debug/askcos_paths.json',        
        time_per_target: int = 15,
    ) -> Dict: 
    
    storage = {} 
    mols = [Chem.MolFromSmiles(smi) for smi in target_smis]

    for mol in tqdm(mols): 
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)  
        params = {
            'smiles': smiles, # required
            'buyable_logic': 'or',
            'max_depth': '10',
            'expansion_time': str(time_per_target),
            'max_ppg': '100',
            'return_first': 'false'
        }
        
        resp = requests.get(host+'/api/treebuilder/', params=params, verify=False)
        trees = resp.json()['trees']
        storage = storage_from_paths(trees, storage)
        
        with open(filename, 'w') as f:
            json.dump(make_dict_jsonable(storage), f, indent="\t")
    
    return storage

def storage_from_api_response(response, storage=None) -> Dict: 
    tree = response['result']['output']
    storage = storage_from_paths(tree, storage)
    return storage 

def storage_from_paths(paths, storage=None) -> Dict:
    for path in paths: 
        storage=update_storage_from_path(path, storage)
    return storage 

def update_storage_from_path(path, storage=None) -> Dict:
    if storage is None: 
        storage = {'Compound Nodes': [],
                   'Reaction Nodes': []}
        stored_rsmis = []
        stored_smis = []
    else: 
        stored_rsmis = [r['smiles'] for r in storage['Reaction Nodes']]
        stored_smis = [c['smiles'] for c in storage['Compound Nodes']]

    if len(path['children']) > 0: 
        for p in path['children']:
            update_storage_from_path(p, storage)

    if ">>" in path['smiles']: 
        update_storage_element(storage['Reaction Nodes'], stored_rsmis, path)
    else:
        update_storage_element(storage['Compound Nodes'], stored_smis, path) 


    return storage 

def update_storage_element(node_list, stored_smiles, path): 
    parents = set([p['smiles'] for p in path['children']])
    # ^ confusing because ASKCOS take products of reaction to be parents

    if path['smiles'] in stored_smiles: 
        entry = np.where(np.array(stored_smiles)==path['smiles'])[0][0]
        node_list[entry]['parents'].update(parents)
    else:
        entry = {
            'smiles': path['smiles'], 
            'parents': parents,
        }
        node_list.append(entry)
    
    return node_list 

def make_dict_jsonable(storage): 
    """ Convert sets in storage dict to lists """
    
    for node_list in storage.values(): 
        for entry in node_list: 
            if 'parents' in entry:
                entry['parents'] = list(entry['parents'])
    
    return storage 

def save_storage_dict(storage, filename): 
    with open(filename, 'w') as f:
        json.dump(make_dict_jsonable(storage), f, indent="\t")
    return 


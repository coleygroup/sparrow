import os
os.environ['KERAS_BACKEND'] = 'theano'

from askcos.retrosynthetic.mcts.tree_builder import MCTS
from askcos.utilities.io.logger import MyLogger
from rdkit import Chem
import tensorflow as tf 
from typing import List, Dict, Union

from sparrow.route_graph import RouteGraph
import json
import numpy as np 
from pathlib import Path

tf.config.set_visible_devices([], 'GPU')
MyLogger.initialize_logFile()
celery = False 


def build_rxn_graph(
        target_smis: List[str],
        n_cpus: int = 5,
        time_per_target: int = 15,
        filename: Union[str, Path] = 'debug/askcos_paths.json'
        ) -> RouteGraph: 
    
    Tree = MCTS(nproc=n_cpus)
    mols = [Chem.MolFromSmiles(smi) for smi in target_smis]
    
    is_first_target = True
    storage = None
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

        for p in Tree.workers:
            if p and p.is_alive():
                p.terminate()
    
    with open(filename, 'w') as f:
        json.dump(make_dict_jsonable(storage), f, indent="\t")

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
            entry['parents'] = list(entry['parents'])
    
    return storage 




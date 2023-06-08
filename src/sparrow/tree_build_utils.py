import os
os.environ['KERAS_BACKEND'] = 'theano'

from askcos.retrosynthetic.mcts.tree_builder import MCTS
from askcos.utilities.io.logger import MyLogger
from rdkit import Chem
import tensorflow as tf 
from typing import List 

from sparrow.route_graph import RouteGraph
import pickle 

tf.config.set_visible_devices([], 'GPU')
MyLogger.initialize_logFile()
celery = False 

def build_rxn_graph(
        target_smis: List[str],
        n_cpus: int = 5,
        time_per_target: int = 15,
        ) -> RouteGraph: 
    
    Tree = MCTS(nproc=n_cpus)

    is_first_target = True

    mols = [Chem.MolFromSmiles(smi) for smi in target_smis]

    for mol in mols:
        
        soft_reset = True
        soft_stop = True 

        # if is_first_target:
        #     soft_reset = True # False
        #     is_first_target = False
        #     soft_stop = False # True
        # else:
        #     soft_reset = True
        #     soft_stop = False

        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)   
        print('expanding target {}'.format(smiles)) 
        paths, stats, graph = Tree.get_buyable_paths(
            smiles,
            nproc=n_cpus,
            expansion_time=time_per_target, 
            termination_logic={'or': ['buyable']},
            # min_chemical_history_dict={'as_reactant':1, 'as_product':1,'logic':'and'},
            soft_reset=soft_reset,
            soft_stop=soft_stop
        )

        with open('paths.pickle', 'wb') as f:
            pickle.dump(paths, f)

        print('done for target {}'.format(smiles))

        Tree.stop()

    
    route_graph = RouteGraph(askcos_MCTS_tree=Tree)

    return route_graph







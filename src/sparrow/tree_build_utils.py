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

    route_graph = RouteGraph()

    mols = [Chem.MolFromSmiles(smi) for smi in target_smis]
    
    is_first_target = True
    for mol in mols:
        
        if is_first_target:
            soft_reset = False 
            is_first_target = False
        else:
            soft_reset = True


        # if is_first_target:
        #     soft_reset = True # False
        #     is_first_target = False
        #     soft_stop = False # True
        # else:
        #     soft_reset = True
        #     soft_stop = False

        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)   
        print('expanding target {}'.format(smiles)) 
        Tree.get_buyable_paths(
            smiles,
            nproc=n_cpus,
            expansion_time=time_per_target, 
            termination_logic={'or': ['buyable']},
            # min_chemical_history_dict={'as_reactant':1, 'as_product':1,'logic':'and'},
            soft_reset=soft_reset,
            soft_stop=True
        )

        print('done for target {}'.format(smiles))

        route_graph.add_from_MCTS(MCTS_tree=Tree)
        
    Tree.stop()

    return route_graph







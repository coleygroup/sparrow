from distutils import core
import os, sys
os.environ['KERAS_BACKEND'] = 'theano'

from askcos.retrosynthetic.mcts.tree_builder import MCTS
from askcos.synthetic.evaluation.tree_evaluator import TreeEvaluator
from askcos.utilities.io.logger import MyLogger
import askcos.global_config as gc
from askcos.utilities.io import name_parser
import pandas as pd
from rdkit import Chem
import pickle
import pprint
import tensorflow as tf 

from src.mars.route_graph import RouteGraph


os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
tf.config.set_visible_devices([], 'GPU')

MyLogger.initialize_logFile()


n_paths = []

celery = False
NCPUS = 5

Tree = MCTS(nproc=NCPUS)


all_routes = []
all_trees = []
all_routes_contexts = []
index_list_all_routes = []
is_first_target = True
status = 0

reaction_dict = {}
for mol in [Chem.MolFromSmiles('O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1')]:
    if is_first_target:
        soft_reset = False
        is_first_target = False
        soft_stop = True
    else:
        soft_reset = True
        soft_stop = False

    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)   
    print('expanding target {}'.format(smiles)) 
    paths, status, graph = Tree.get_buyable_paths(smiles,
                                        nproc=NCPUS,
                                        expansion_time=10, 
                                        termination_logic={'or': ['buyable']},
                                        # min_chemical_history_dict={'as_reactant':1, 'as_product':1,'logic':'and'},
                                        soft_reset=soft_reset,
                                        soft_stop=soft_stop)
    print('done for target {}'.format(smiles))

    # route_graph = RouteGraph(askcos_MCTS_tree=Tree)

    for chemical in Tree.Chemicals:
        chemical_obj = Tree.Chemicals[chemical]
        for _id, result in chemical_obj.template_idx_results.items():
            if result.waiting:
                continue
            for rsmi,R in result.reactions.items():
                if (not R.valid) or R.price == -1:
                    continue
                rxn = '.'.join(sorted(R.reactant_smiles)) + '>>' + chemical
                if rxn not in reaction_dict:
                    reaction_dict[rxn]={}
                    for rct in rsmi.split('.'):
                        reaction_dict[rxn][rct]=-1
                    reaction_dict[rxn][chemical]=1
    print('number of reactions explored: {}'.format(len(reaction_dict)))
    feasible_trees = []
    n_paths.append(len(paths))






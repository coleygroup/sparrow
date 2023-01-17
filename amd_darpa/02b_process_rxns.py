import os, sys
os.environ['KERAS_BACKEND'] = 'theano'
os.environ["CUDA_VISIBLE_DEVICES"] = '-1'
sys.path.append('/home/jfromer/Molecule_library_synthesis/')

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1' 

import tensorflow as tf
tf.autograph.experimental.do_not_convert
tf.config.set_visible_devices([], 'GPU')
tf.keras.backend.set_floatx('float32')


import pickle
import pandas as pd
from rdkit import Chem
# from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if mo

import askcos.utilities.contexts as context_cleaner
from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
from askcos.synthetic.evaluation.evaluator import Evaluator
from pathlib import Path
from tqdm import tqdm
import askcos.global_config as gc
import collections
import sys
import json


rxns_path = Path('amd_darpa/second_sulfur_explore_scaffold_round4_reactions_for_grouping.json')
with open(rxns_path, 'r') as file: 
    rxns = json.load(file)

###label rxns
i=0
rxn_le = {}
encoded_rxn_dict = {}
count = 0
evaluator = Evaluator() 
cond_dict = collections.defaultdict(list)

for tree in tqdm(rxns.values()):
    for reactions in tree.values():
        for rxn in reactions.values(): 
            try:
                rxn_le[count] = rxn['reaction']
                stoich = { 
                    **{smi: -1 for smi in rxn['reactants']},
                    **{smi: 1 for smi in rxn['products']}
                }
                encoded_rxn_dict[count] = stoich
                
                # put together contexts
                contexts = []
                for ct in rxn['contexts'].values(): 
                    this_ct = [ct['temperature']]

                    if type(ct['solvent']) == list:
                        this_ct.append(('.'.join(ct['solvent'])))
                    elif ct['solvent']:
                        this_ct.append(ct['solvent'])
                    else: 
                        this_ct.append('')

                    if type(ct['reagent']) == list:
                        this_ct.append(('.'.join(ct['reagent'])))
                    elif ct['reagent']:
                        this_ct.append(ct['reagent'])
                    else: 
                        this_ct.append('')

                    if type(ct['catalyst']) == list:
                        this_ct.append(('.'.join(ct['reagent'])))
                    elif ct['catalyst']:
                        this_ct.append(ct['catalyst'])
                    else: 
                        this_ct.append('')

                    contexts.append(this_ct)

                rsmi, psmi = rxn['reaction'].split('>>')
                contexts = [context_cleaner.clean_context(context) for context in contexts]
                try:
                    eval_res = evaluator.evaluate(rsmi, psmi, contexts,  )
                except Exception as e:
                    print(e)
                    eval_res = [{'target':{'prob':0}}]*len(contexts)
                encoded_rxn_dict[count]['cond'] = {}
                encoded_rxn_dict[count]['score'] = {}
                for i in range(len(eval_res)):
                    eval_res[i]['context'] = contexts[i]
                    encoded_rxn_dict[count]['cond'][i] = eval_res[i]['context']
                    encoded_rxn_dict[count]['score'][i] = eval_res[i]['target']['prob']
                    for cond in contexts[i]:
                        if cond is not "":
                            cond_dict[cond].append((count,i))
                
                count = count + 1

            except BaseException as e:
                print(e)

                     
with open((rxns_path.parent / 'encoded_rxn_dict.pkl'),'wb') as RXN_DICT:
    pickle.dump(encoded_rxn_dict, RXN_DICT)

with open((rxns_path.parent / 'cond_dict.pkl'),'wb') as COND_DICT:
    pickle.dump(cond_dict, COND_DICT)

with open((rxns_path.parent / 'rxn_le.pkl'),'wb') as RXN_LE:
    pickle.dump(rxn_le, RXN_LE)
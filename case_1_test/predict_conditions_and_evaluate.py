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
import askcos.global_config as gc
import collections
import sys

case =  'case_1_test/'
with open(case+'reaction_dict.pickle','rb') as RD:
    rxn_dict = pickle.load(RD)

target_df = pd.read_csv(case+'target_list.csv')
target_list = list(target_df.loc[target_df['numberofpaths']>0]['SMILES'])
target_dict = set([Chem.MolToSmiles(Chem.MolFromSmiles(target)) for target in target_list])
small_target_dict = target_dict

###label rxns
i=0
rxn_le = {}
encoded_rxn_dict = {}
for key,value in rxn_dict.items():
    rxn_le[i]=key
    encoded_rxn_dict[i]=value
    i+=1

##load conditionprediction and forward evaluation models
cont = NeuralNetContextRecommender()
cont.load_nn_model(model_path=gc.NEURALNET_CONTEXT_REC['model_path'], info_path=gc.NEURALNET_CONTEXT_REC[
                       'info_path'], weights_path=gc.NEURALNET_CONTEXT_REC['weights_path'])
evaluator = Evaluator()

from tqdm import tqdm
cond_dict = collections.defaultdict(list)

ctr = 0

for key, value in tqdm(encoded_rxn_dict.items()):
    print(key)
    rxn = rxn_le[key]
    rsmi = rxn.split('>>')[0]
    psmi = rxn.split('>>')[1]
    uncleaned_contexts = cont.get_n_conditions(rxn, n=10, return_separate=True, return_scores=False) # from return_separate (?)
    contexts = cont.get_n_conditions(rxn, n=10, return_separate=False, return_scores=False)
    contexts = [context_cleaner.clean_context(context) for context in contexts]
    print(contexts[0])
    try:
        eval_res = evaluator.evaluate(rsmi, psmi, contexts,  )
    except Exception as e:
        print(e)
        eval_res = [{'target':{'prob':0}}]*len(uncleaned_contexts)
    encoded_rxn_dict[key]['cond'] = {}
    encoded_rxn_dict[key]['score'] = {}
    for i in range(len(eval_res)):
        eval_res[i]['context'] = uncleaned_contexts[i][:5]
        encoded_rxn_dict[key]['cond'][i] = eval_res[i]['context']
        encoded_rxn_dict[key]['score'][i] = eval_res[i]['target']['prob']
        for cond in uncleaned_contexts[i][:5]:
            if cond is not "":
                cond_dict[cond].append((key,i))
    # ctr+=1
    # if ctr>=10:
    #     break

                        
with open(case+'/encoded_rxn_dict_with_cond.pkl','wb') as RXN_DICT:
    pickle.dump(encoded_rxn_dict, RXN_DICT)

with open(case+'/cond_dict.pkl','wb') as COND_DICT:
    pickle.dump(cond_dict, COND_DICT)
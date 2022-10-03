import pickle
import pandas as pd
from rdkit import Chem
import collections
import urllib
import argparse
from tqdm import tqdm
from askcos.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI
from typing import List, bool, str
# from gurobipy import *

class Optimizer():
    def __init__(self,
            weights:List[int] = [1,1,1,10],
            constrain_all_targets: bool = 1, 
            case_name: str = 'case_1_Jenna',
            reaction_dict_filename: str = None,
            target_list_filename: str = None, 
            encoded_rxn_filename: str = None,
            condition_dict_filename: str = None,
        ) -> None:
        
        self.weights = weights
        self.constrain_all_targets = constrain_all_targets
        self.reaction_dict_filename = reaction_dict_filename or \
            case_name + '/reaction_dict.pickle'
        self.target_list_filename = target_list_filename or \
            case_name + '/target_list.csv'
        self.encoded_rxn_filename = encoded_rxn_filename or \
            case_name + '/encoded_rxn_dict_with_cond.pkl'
        self.condition_dict_filename = condition_dict_filename or \
            case_name + '/cond_dict.pkl'
        self.case_name = case_name 

        self.rxn_dict = self.load_reactions()
        self.target_dict = self.load_targets()
        self.encoded_rxn_dict = self.load_evaluated_rxns()
        self.cond_dict = self.load_cond_dict()


    def load_reactions(self):
        with open(self.reaction_dict_filename,'rb') as RD:
            rxn_dict = pickle.load(RD)
        return rxn_dict

    def load_targets(self):
        target_df = pd.read_csv(self.target_list_filename)
        target_list = list(target_df.loc[target_df['numberofpaths']>0]['SMILES'])
        target_dict = set([Chem.MolToSmiles(Chem.MolFromSmiles(target),False) for target in target_list])
        print('Loaded {} targets!'.format(len(target_list)))
        return target_dict 

    def load_evaluated_rxns(self):
        # load 
        with open(self.encoded_rxn_filename,'rb') as RXN_DICT:
            try:
                encoded_rxn_dict = pickle.load(RXN_DICT,encoding="byte")
            except: 
                encoded_rxn_dict = pickle.load(RXN_DICT,encoding="latin1")
        
        # add penalties
        for key,value in encoded_rxn_dict.items():
            try:
                encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
                    for key,score in encoded_rxn_dict[key]['score'].items()}
            except:
                encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
                    for key,score in encoded_rxn_dict[key]['score'.encode('utf-8')].items()}

        return encoded_rxn_dict  

    def load_cond_dict(self):
        with open(self.condition_dict_filename,'rb') as COND_DICT:
            cond_dict = pickle.load(COND_DICT,encoding="bytes")
        return cond_dict

    def encode_reactions(self):
        encoded_rxn_dict = {}
        for i, (key,value) in enumerate(self.rxn_dict.items()):
            encoded_rxn_dict[i]=value
        return encoded_rxn_dict

    

    def 
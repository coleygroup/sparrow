""" Creates dictionaries needed for optimization and outputs them in \dictionaries folder """
import argparse
import numpy as np 
from pathlib import Path 
import pickle 
import csv
from rdkit import Chem

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dict-folder", 
        help="folder containing encoded_rxn_dict and cond_dict", 
        default="amd_darpa", 
        action="store"
    )

    return parser.parse_args()

def load_targets_and_rewards(folder): 
    rewards = {}
    with open(dict_folder  / 'target_list.csv','r') as data:
        for line in csv.reader(data):
            if line[0] != 'smiles':
                rewards[line[0]] = float(line[1])

    target_list = list(rewards.keys())
    print(f"Loaded {len(target_list)} targets and rewards")
    target_dict = set(target_list)
    
    return rewards, target_list, target_dict

def calculate_penalties(encoded_rxn_dict):
    for key,value in encoded_rxn_dict.items():
        try:
            encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
                for key,score in encoded_rxn_dict[key]['score'].items()}
        except:
            encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
                for key,score in encoded_rxn_dict[key]['score'.encode('utf-8')].items()}

    return encoded_rxn_dict


def load_dicts(dict_folder): 

    with open(dict_folder / 'encoded_rxn_dict.pkl','rb') as RXN_DICT:
        try:
            encoded_rxn_dict = pickle.load(RXN_DICT,encoding="byte")
        except: 
            encoded_rxn_dict = pickle.load(RXN_DICT,encoding="latin1")
    
    encoded_rxn_dict = calculate_penalties(encoded_rxn_dict)
	
    with open(dict_folder / 'cond_dict.pkl','rb') as COND_DICT:
        cond_dict = pickle.load(COND_DICT,encoding="bytes")
        
    return encoded_rxn_dict, cond_dict

def construct_chemical_dict(encoded_rxn_dict, target_dict, cond_dict): 
    storage = {}
    chem_dict_start = {}
    chem_dict_tar = {}
    chem_dict_inter = {}
    chem_le={}
    chem_le_rev={}
    _id = 0
    for key,value in encoded_rxn_dict.items():
        for c, stoi in value.items():
            if c=='cond' or c=='score' or c=='penalty':
                continue
            if c not in chem_le:
                chem_le[c]=_id
                chem_le_rev[_id]=c
                _id+=1
            if c in chem_dict_start:
                chem_dict_start[c][1].append(key)
                continue
            elif c in chem_dict_tar:
                chem_dict_tar[c].append(key)
                continue
            elif c in chem_dict_inter:
                if stoi==1:
                    chem_dict_inter[c][0].append(key)
                else:
                    chem_dict_inter[c][1].append(key)
                continue
            elif c in target_dict:
                chem_dict_tar[c]=[key]
            else:
                m = Chem.MolFromSmiles(c)
                nc = 0
                nn = 0
                no = 0
                for a in m.GetAtoms():
                    if a.GetSymbol() == 'C':
                        nc +=1
                    if a.GetSymbol() == 'N':
                        nn +=1
                    if a.GetSymbol() == 'O':
                        no +=1
                if nc<=10 and nn<=3 and no<=5:
                    chem_dict_start[c]=[[],[key]]
                else:
                    if stoi ==1:
                        chem_dict_inter[c] = [[key],[]]
                    else:
                        chem_dict_inter[c] = [[],[key]]

    storage['encoded_rxn_dict'] = encoded_rxn_dict
    storage['chem_dict_start'] = chem_dict_start
    storage['chem_dict_inter'] = chem_dict_inter
    storage['chem_dict_tar'] = chem_dict_tar
    storage['encoded_cond_dict'], storage['cond_le'], storage['cond_le_rev'] = encode_conditions(cond_dict)    
    storage['target_dict'] = target_dict
    
    return storage

def encode_conditions(cond_dict):
    i=0
    cond_le = {}
    cond_le_rev = {}
    encoded_cond_dict = {}
    for key, value in cond_dict.items():
        cond_le[i]=key
        cond_le_rev[key]=i
        encoded_cond_dict[i]=value
        i+=1
    
    return encoded_cond_dict, cond_le, cond_le_rev

def encode_chem(chem_dict, chem_le):
    encoded_dict = {}
    for key, value in chem_dict.items():
        encoded_dict[chem_le[key]]=value
    
    return encoded_dict

def save_dicts(storage, out_folder): 
    if not out_folder.is_dir(): 
        out_folder.mkdir()
    for name, value in storage.items(): 
        out_path = out_folder / f'{name}.pkl'
        with open(out_path, 'wb') as f: 
            pickle.dump(value, f)

if __name__ == '__main__': 
    args = get_args()
    dict_folder = Path(args.dict_folder)
    out_folder = dict_folder / "dictionaries"
    rewards, target_list, target_dict = load_targets_and_rewards(dict_folder)
    encoded_rxn_dict, cond_dict = load_dicts(dict_folder)
    storage = construct_chemical_dict(encoded_rxn_dict, target_dict, cond_dict)
    storage['rewards'] = rewards
    save_dicts(storage, out_folder)
    print('done exporting dictionaries')
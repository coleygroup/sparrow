import os, sys
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located


import argparse 
from pathlib import Path 
from sparrow.json_utils import build_retro_graph_local
import json 
import pandas as pd 



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-dict", action="store", type=str)
    parser.add_argument("--base-dir", action="store", type=str)
    parser.add_argument("--storage", action="store", type=str)
    parser.add_argument("--n", default=0, action="store", type=int)
    parser.add_argument("--sec", default=30, action="store", type=int)
    return parser.parse_args()

def load_targets(filepath): 
    df = pd.read_csv(filepath)
    targets = list(df['SMILES'])
    print(len(targets))
    return targets

def load_storage(filepath): 
    if filepath.exists(): 
        with open(filepath, 'r') as f: 
            storage = json.load(f)
        return storage 
    else: 
        return None


if __name__=='__main__':
    args = get_args()
    kwargs = args.__dict__
    base_dir = Path(kwargs['base_dir'])
    targets = load_targets(base_dir / kwargs['target_dict'])
    storage = load_storage(base_dir / kwargs['storage'])
    tar = targets[kwargs['n']]
    storage = build_retro_graph_local(        
        target_smis=[tar],
        filename=base_dir / kwargs['storage'],        
        time_per_target=kwargs['sec'],
        storage=storage,
    )

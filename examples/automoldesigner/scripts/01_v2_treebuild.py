""" Builds retrosynthesis trees using ASKCOSv2 for all molecules """

import copy
import numpy as np
import requests
import pandas as pd 
from pathlib import Path 
import json 
from tqdm import tqdm 
from rdkit import Chem

query_template = {
    "smiles": "",
    "expand_one_options": {
        "template_max_count": 100,
        "template_max_cum_prob": 0.995,
        "banned_chemicals": [],
        "banned_reactions": [],
        "retro_backend_options": [
            {
                "retro_backend": "template_relevance",
                "retro_model_name": "reaxys",
                "max_num_templates": 100,
                "max_cum_prob": 0.995,
                "attribute_filter": []
            }
        ],
        "use_fast_filter": True,
        "filter_threshold": 0.75,
        "cluster_precursors": False,
        "cluster_setting": {
            "feature": "original",
            "cluster_method": "hdbscan",
            "fp_type": "morgan",
            "fp_length": 512,
            "fp_radius": 1,
            "classification_threshold": 0.2
        },
        "extract_template": False,
        "return_reacting_atoms": False,
        "selectivity_check": False
    },
    "build_tree_options": {
        "expansion_time": 60,
        "max_branching": 25,
        "max_depth": 6,
        "exploration_weight": 1,
        "return_first": False
    },
    "enumerate_paths_options": {
        "path_format": "json",
        "json_format": "nodelink",
        "sorting_metric": "plausibility",
        "validate_paths": True,
        "score_trees": False,
        "cluster_trees": False,
        "cluster_method": "hdbscan",
        "min_samples": 5,
        "min_cluster_size": 5,
        "paths_only": False,
        "max_paths": 500,
    },
}

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def main():
    HOST = "" # insert host address
    PORT = "" # insert port
    output_dir = Path('examples/automoldesigner/trees')
    output_dir.mkdir(exist_ok=True, parents=True)

    smis_data = pd.read_csv('examples/automoldesigner/staph_targets_raw.csv')
    smiles = list(smis_data.SMILES)

    successes = []

    for i, smi in tqdm(enumerate(smiles)):
        output_file = output_dir / f'tree_{i:04d}.json'
        if output_file.exists(): 
            continue

        smi_neutral = Chem.MolToSmiles(neutralize_atoms(Chem.MolFromSmiles(smi)))
        data = copy.deepcopy(query_template)
        data["smiles"] = smi_neutral.strip()

        resp = requests.post(
            url=f"http://{HOST}:{PORT}/api/tree-search/mcts/call-sync-without-token",
            json=data
        ).json()
        
        total_paths = resp["result"]["stats"]["total_paths"]
        print(f"SMILES: {smi_neutral.strip()}, number of paths: {total_paths: .0f}")

        success = total_paths > 0
        successes.append(success)

        reactions = set()
        paths = resp["result"]["paths"]
        for path in paths: 
            rxnsmis = [node['smiles'] for node in path['nodes'] if '>>' in node['smiles']]
            reactions.update(rxnsmis)

        reactions = list(reactions)

        with open(output_file, 'w') as f:
            json.dump(reactions, f, indent='\t')

    print(f"Success rate: {np.mean(successes)}")

if __name__=='__main__': 
    main()


    

 

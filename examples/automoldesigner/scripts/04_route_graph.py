""" Combines reactions, conditions, and scores into a json file input for SPARROW """
import json 
import pandas as pd
from rdkit import Chem

thresh = 0.01 

# load reaction conditions 
with open('examples/automoldesigner/chkpts/conditions.json', 'r') as f: 
    conditions = json.load(f)

# load reaction scores 
with open('examples/automoldesigner/chkpts/scores.json', 'r') as f: 
    scores = json.load(f)

# save json file, exclude reactions with score = 0 
nodes = {
    "Compound Nodes": [], 
    "Reaction Nodes": [{
        "smiles": smi,
        "condition": conditions[smi],
        "score": scores[smi],
    } for smi in scores if scores[smi]>thresh]
}

with open('examples/automoldesigner/chkpts/graph_reactions_only.json', 'w') as f: 
    json.dump(nodes, f, indent='\t')

# generates target csv file
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

raw_data = pd.read_csv("examples/automoldesigner/staph_targets_raw.csv")
neutral_smiles = [Chem.MolToSmiles(neutralize_atoms(Chem.MolFromSmiles(smi))) for smi in raw_data.SMILES]
raw_data["Neutral SMILES"] = neutral_smiles
raw_data.to_csv("examples/automoldesigner/chkpts/neutral_smiles.csv", index=False)

# eliminate compounds that have no routes 
all_rxn_products = set([rxn.split('>>')[1] for rxn, s in scores.items() if s>thresh])
targets_df = raw_data[raw_data["Neutral SMILES"].isin(all_rxn_products)]
targets_df = targets_df[['Neutral SMILES', 'Activity probability', 'Closest known molecule']]
targets_df = targets_df.rename(columns={
    'Neutral SMILES': 'SMILES', 
    'Activity probability': 'Reward', 
})
targets_df.to_csv('examples/automoldesigner/targets.csv', index=False)


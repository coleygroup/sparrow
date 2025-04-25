""" Combines reactions, conditions, and scores into a json file input for SPARROW """
import json 
import pandas as pd
from rdkit import Chem
from sparrow.coster import LookupCoster
from sparrow.route_graph import RouteGraph

thresh = 0.01 

# load reaction conditions 
with open('examples/swanson/chkpts/conditions.json', 'r') as f: 
    conditions = json.load(f)

# load reaction scores 
with open('examples/swanson/chkpts/scores.json', 'r') as f: 
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

with open('examples/swanson/chkpts/graph_reactions_only.json', 'w') as f: 
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

raw_data = pd.read_csv("examples/swanson/targets_raw.csv")
neutral_smiles = [Chem.MolToSmiles(neutralize_atoms(Chem.MolFromSmiles(smi))) for smi in raw_data.SMILES]
raw_data["Neutral SMILES"] = neutral_smiles
raw_data.to_csv("examples/swanson/chkpts/neutral_smiles.csv", index=False)

# eliminate compounds that have no routes 
all_rxn_products = set([rxn.split('>>')[1] for rxn, s in scores.items() if s>thresh])
targets_df = raw_data[raw_data["Neutral SMILES"].isin(all_rxn_products)]
targets_df = targets_df[['Neutral SMILES', 'Reward']]
targets_df = targets_df.rename(columns={
    'Neutral SMILES': 'SMILES', 
})
targets_df.to_csv('examples/swanson/targets.csv', index=False)


graph = RouteGraph(node_filename='examples/swanson/chkpts/graph_reactions_only.json')

coster = LookupCoster(
    lookup_file='examples/automoldesigner/data/enamine_per_g_122024.csv'
)
graph.set_buyable_compounds_and_costs(
    coster=coster,
    save_json_dir='examples/swanson/chkpts',
    save_freq=1e10
)
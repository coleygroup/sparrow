import pandas as pd 
from rdkit.Chem import RDConfig
from rdkit import Chem
from pathlib import Path 
import json 
import os, sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


df = pd.read_csv('examples/garibsingh/targets.csv')
df['SA Score'] = [round(sascorer.calculateScore(Chem.MolFromSmiles(smi)), 4) for smi in df['SMILES']]
df['Combined Score'] = [round(rew-sa , 4) for sa, rew in zip(df['SA Score'], df['Reward'])]
df['#'] = [i for i in range(len(df))]

def cost_from_trees(trees): 
    # creat set of reactions 
    reactions = {}
    bbs = {}
    for tree in trees: 
        for rnode in tree['Reaction Nodes']: 
            reactions[rnode['smiles']] = rnode['score']
        for cnode in tree['Compound Nodes']: 
            if 'terminal' in cnode:
                bbs[cnode['smiles']] = cnode['cost']

    N_rxns = len(reactions)
    average_rxn_score = sum(reactions.values())/N_rxns
    bb_cost = sum(bbs.values())
    return N_rxns, average_rxn_score, bb_cost

tree_dir = Path('examples/garibsingh/baselines/trees')
all_trees = []
for i in range(len(df)):
    with open(tree_dir/f'tree_{i}.json') as f:
        all_trees.append(json.load(f))


top_combined = df.sort_values(by='Combined Score', ascending=False, inplace=False)
top_reward = df.sort_values(by='Reward', ascending=False, inplace=False)
top_sa = df.sort_values(by='SA Score', inplace=False)

strategies = ['Reward Only', 'SA Score Only', 'Combined Score']
dfs = [top_reward, top_sa, top_combined]
entries = []
for strategy, data in zip(strategies, dfs):
    for N in range(1,len(df)+1): 
        compounds = data['#'][:N]
        trees = [all_trees[i] for i in compounds]
        N_rxns, average_rxn_score, bb_cost = cost_from_trees(trees)
        entries.append({
            'strategy': strategy,
            'N': N,
            'Reward': round(sum(data['Reward'][:N]), 4),
            'Cost': bb_cost,
            'N_rxns': N_rxns,
            'Average reaction score': average_rxn_score,
        })

results = pd.DataFrame.from_dict(entries)
results.to_csv('examples/garibsingh/baselines/baseline_curves.csv', index=False)
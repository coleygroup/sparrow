import pandas as pd 
from rdkit.Chem import RDConfig
from rdkit import Chem
from pathlib import Path 
import json 
import os, sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


df = pd.read_csv('examples/amd/targets.csv')
df['SA Score'] = [round(sascorer.calculateScore(Chem.MolFromSmiles(smi)), 4) for smi in df['SMILES']]
df['Combined Score'] = [round(rew-sa , 4) for sa, rew in zip(df['SA Score'], df['Reward'])]
df['#'] = [i for i in range(len(df))]

def cost_from_trees(trees, rewards): 
    # creat set of reactions 
    reactions = {}
    bbs = {}
    expected_reward = 0
    for tree, rew in zip(trees, rewards): 
        p_success = 1
        for rnode in tree['Reaction Nodes']: 
            reactions[rnode['smiles']] = rnode['score']
            p_success = p_success*rnode['score']
        for cnode in tree['Compound Nodes']: 
            if 'terminal' in cnode:
                bbs[cnode['smiles']] = cnode['cost']
        expected_reward += p_success*rew

    N_rxns = len(reactions)
    average_rxn_score = sum(reactions.values())/N_rxns
    bb_cost = sum(bbs.values())
    return N_rxns, average_rxn_score, bb_cost, expected_reward

tree_dir = Path('examples/amd/baselines/trees')
all_trees = []
for i in range(len(df)):
    tree_file = tree_dir/f'tree_{i}.json'
    if tree_file.exists():
        with open(tree_file) as f:
            all_trees.append(json.load(f))
    else: 
        all_trees.append({})

no_route = [i for i in range(len(all_trees)) if all_trees[i]=={}]
df = df[~df['#'].isin(no_route)]

top_combined = df.sort_values(by='Combined Score', ascending=False, inplace=False)
top_reward = df.sort_values(by='Reward', ascending=False, inplace=False)
top_sa = df.sort_values(by='SA Score', inplace=False)

strategies = ['Reward Only', 'SA Score Only', 'Combined Score']
dfs = [top_reward, top_sa, top_combined]
entries = []
for strategy, data in zip(strategies, dfs):
    for N in range(1,len(data)+1, 10): 
        compounds = data['#'][:N]
        rewards = data['Reward'][:N]
        trees = [all_trees[i] for i in compounds]
        N_rxns, average_rxn_score, bb_cost, expected_reward = cost_from_trees(trees, rewards)
        entries.append({
            'strategy': strategy,
            'N': N,
            'Reward': round(sum(data['Reward'][:N]), 4),
            'Cost': bb_cost,
            'N_rxns': N_rxns,
            'Average reaction score': average_rxn_score,
            'Expected Reward': expected_reward, 
        })

results = pd.DataFrame.from_dict(entries)
results.to_csv('examples/amd/baselines/baseline_curves.csv', index=False)
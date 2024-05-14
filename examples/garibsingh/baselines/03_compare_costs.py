import json 
from pathlib import Path 
import pandas as pd 

def create_buyables_from_sparrow(tree_file='examples/garibsingh/trees_w_info.json'):
    with open(tree_file, 'r') as f:
        tree = json.load(f)
    
    buyables = {}
    for node in tree['Compound Nodes']: 
        if node['buyable']: 
            buyables[node['smiles']] = node['cost_per_g']
    
    return buyables

def create_buyable_from_trees(tree_folder='examples/garibsingh/baselines/trees'): 
    buyables = {}
    for tree_file in Path(tree_folder).glob('*.json'): 
        with open(tree_file, 'r') as f:
            tree = json.load(f)
        
        for node in tree['Compound Nodes']: 
            if 'terminal' in node and node['terminal']: 
                buyables[node['smiles']] = node['cost']
    
    return buyables

buyables_2023 = create_buyables_from_sparrow()
buyables_2024 = create_buyable_from_trees()
common = set(buyables_2024).intersection(buyables_2023)
prices = {smi: [buyables_2023[smi], buyables_2024[smi]] for smi in common}
df = pd.DataFrame.from_dict(prices, orient='index', columns=['2023', '2024'])
df.to_csv('examples/garibsingh/baselines/price_comparison.csv')

print('')


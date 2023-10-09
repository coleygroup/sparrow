import csv 
from pathlib import Path
from tqdm import tqdm
from rdkit import Chem

with open('darpa_example/test_mols.csv','r') as csvfile:
    targets = [row[0] for row in csv.reader(csvfile)]

targets = targets[1:]

found_targets = 0

for target in tqdm(targets, desc='Searching in trees for targets'): 
    tar_found = 0
    clean_target = Chem.MolToSmiles(Chem.MolFromSmiles(target))
    for treefile in Path('darpa_example/trees').glob('*.json'):
        with open(treefile, 'r') as f: 
            if target in f.read():
                tar_found = 1
    if tar_found: 
        found_targets += 1

print(f'{found_targets} targets found in trees')

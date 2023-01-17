import pickle
from pathlib import Path
from tqdm import tqdm
import sys
import json


rxns_path = Path('amd_darpa/second_sulfur_explore_scaffold_round4_reactions_for_grouping.json')
with open(rxns_path, 'rb') as file: 
    rxns = json.load(file)

###label rxns
rxn_dict = {}
count = 0
# evaluator = Evaluator() 
# cond_dict = collections.defaultdict(list)

for tree in tqdm(rxns.values()):
    for reactions in tree.values():
        for rxn in reactions.values(): 
            try:
                stoich = { 
                    **{smi: -1 for smi in rxn['reactants']},
                    **{smi: 1 for smi in rxn['products']}
                }
                rxn_dict[rxn['reaction']] = stoich
                count = count + 1
            except BaseException as e:
                print(e)


with open((rxns_path.parent / 'rxn_dict.pickle'),'wb') as RD:
    pickle.dump(rxn_dict, RD)

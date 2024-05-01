import json 
import pandas as pd

results_file = 'results/button_alectinib/30_1_5/optimal_routes.json'
output_file = 'results/button_alectinib/30_1_5/reaction_info.csv'

def results_to_csv(results_file, output_file): 
    with open(results_file, 'r') as f:
        results = json.load(f)
    
    rxn_smis = []
    scores = []
    conditions = []

    for entry in results['Reactions']: 
        if entry['smiles'].startswith('>>'): 
            continue 

        rxn_smis.append(entry['smiles'])
        scores.append(entry['score'])
        conditions.append('.'.join(entry['conditions']))

    df = pd.DataFrame({
        'Reaction SMILES': rxn_smis,
        'ASKCOS score': scores, 
        'ASKCOS conditions': conditions, 
    })

    df.to_csv(output_file, index=False)


if __name__=='__main__':
    results_to_csv(results_file, output_file)
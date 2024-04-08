import requests
import pandas as pd
import time 
import urllib3
import json 
from rdkit import Chem
from pathlib import Path 

from sparrow.coster import ChemSpaceCoster
from sparrow.condition_recommender import AskcosAPIRecommender
from sparrow.scorer import AskcosAPIScorer

chemspace_api_key = ''
host = ''
urllib3.disable_warnings()

def post_and_get(host, params, sleep_time = 10, timeout = 650):
    req = requests.post(host+'/api/v2/tree-builder/', data=params, verify=False)
    try:
        task_id = req.json()['task_id']
    except KeyError:
        return {} , {}
    results = requests.get(host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return req.json(), results.json()

def check_path_valid(path, coster): 
    compound_nodes = []
    reaction_nodes = []
    valid = True 
    for node in path: 
        smi = node['smiles']
        if '>>' in smi:
            reaction_nodes.append({'smiles': smi})
        else: 
            if node['terminal']: 
                buyable, cost = coster.get_buyable_and_cost(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
                compound_nodes.append({'smiles': smi, 'cost': cost, 'terminal': True})
                if not buyable: 
                    valid = False 
            else: 
                compound_nodes.append({'smiles': smi})

    return valid, reaction_nodes, compound_nodes


# load targets.csv
df = pd.read_csv('examples/garibsingh/targets.csv')

# perform retrosynthesis using askcos
params = {
    'buyable_logic': 'and',
    'max_depth': '5',
    'expansion_time': str(60),
    'max_ppg': str(10000),
    'return_first': 'false', 
    'max_branching': str(20),
    'paths_only': 'true',
    'json_format': 'nodelink'
}

def post_and_get(host, params, sleep_time = 10, timeout = 650):
    req = requests.post(host+'/api/v2/tree-builder/', data=params, verify=False)
    try:
        task_id = req.json()['task_id']
    except KeyError:
        return {} , {}
    results = requests.get(host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return req.json(), results.json()

def check_path_valid(path, coster): 
    compound_nodes = []
    reaction_nodes = []
    valid = True 
    for node in path: 
        smi = node['smiles']
        if '>>' in smi:
            reaction_nodes.append({'smiles': smi})
        else: 
            if node['terminal']: 
                buyable, cost = coster.get_buyable_and_cost(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
                compound_nodes.append({'smiles': smi, 'cost': cost, 'terminal': True})
                if not buyable: 
                    valid = False 
            else: 
                compound_nodes.append({'smiles': smi})

    return valid, reaction_nodes, compound_nodes

savedir = Path('examples/garibsingh/baselines/trees')
savedir.mkdir(exist_ok=True)
coster = ChemSpaceCoster(api_key=chemspace_api_key)
coster.build_status_log()

recommender = AskcosAPIRecommender(host=host)
scorer = AskcosAPIScorer(host=host)

for i, smi in enumerate(df['SMILES']): 
    treefile = savedir/f'tree_{i}.json'
    params['smiles'] = smi
    request, result = post_and_get(host, params)
    paths = result['output']

    valid = False 
    c = 0
    # use first path where all BBs are buyable according to ChemSpace
    while valid == False and c <= len(result['output']):
        path = result['output'][c]['nodes']
        valid, reaction_nodes, compound_nodes = check_path_valid(path, coster)
        c += 1

    if not valid: 
        print(f'no valid route to {smi}!!!!')
        continue 
    
    for rnode in reaction_nodes: 
        rnode['condition'] = recommender.recommend_conditions(rnode['smiles'], n_c=1)
        rnode['score'] = scorer.score_rxn(rxn_smi=rnode['smiles'], condition=rnode['condition'])

    tree = {}
    tree['Reaction Nodes'] = reaction_nodes
    tree['Compound Nodes'] = compound_nodes
    with open(treefile, 'w') as f: 
        json.dump(tree, f, indent='\t')
    
    print(f'{len(tree["Compound Nodes"])} compounds, {len(tree["Reaction Nodes"])} reactions')









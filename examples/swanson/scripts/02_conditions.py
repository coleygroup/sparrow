""" Predicts reaction conditions for all reactions in retrosynthetic trees using ASKCOS v2 context recommender """

import json 
import copy 
import time 
import requests 

from pathlib import Path 
from tqdm import tqdm 

HOST = "" # rename with askcos host
PORT = "" # rename with port 

# functions, templates 
query_template = {
    "reactants": "", 
    "products": "", 
    "return_scores": False, 
    "with_smiles": False, 
    "num_results": 1, 
}

def post(host_post, params):
    req = requests.post(host_post, json=params, verify=False)
    try: 
        task_id = req.json()['task_id']
        return task_id
    except KeyError: 
        return None

def get(host_results, task_id, sleep_time = 10, timeout = 650): 
    if task_id is None: 
        return {}, {}
    results = requests.get(host_results + 'task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(host_results + 'task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return results.json()

# load all reactions 
rxn_dirs = Path('examples/swanson/trees').glob('*.json')
allrxns = set()
for rxnfile in rxn_dirs:
    with open(rxnfile, 'r') as f:
        rxns = json.load(f)
    allrxns.update(rxns)

# load all pre-evaluated conditions 
cond_path = Path('examples/swanson/chkpts/conditions.json')
if cond_path.exists():
    with open(cond_path, 'r') as f:
        conditions = json.load(f)
else: 
    conditions = {}

# assign conditions to not previously evaluated reactions 
rxns_eval = [rxn for rxn in allrxns if rxn not in conditions or conditions[rxn]==[[]]]
batch_size = 10
rxn_chunks = [rxns_eval[i:i+batch_size] for i in range(0, len(rxns_eval), batch_size)] 

i = 0
for chunk in tqdm(rxn_chunks, position=0, leave=True):
    task_ids = {}

    # post for each reaction in the chunk 
    for rxn in chunk:   
        params = copy.deepcopy(query_template)
        params["reactants"], params["products"] = rxn.split('>>')
        task_ids[rxn] = post(host_post=f'http://{HOST}:{PORT}/api/legacy/context/', params=params)
    
    # get and update nodes 
    for rxn in chunk: 
        result = get(
            host_results=f'http://{HOST}:{PORT}/api/legacy/celery/',
            task_id=task_ids[rxn],
            sleep_time=0.1,
            timeout=20,
        )
        try: 
            pred = result['output'][0]
        except: 
            print(f'Error with reaction {rxn}, returning empty conditions')
            pred = {}
        
        rxn_cond = [pred[entry] for entry in ['solvent', 'reagent', 'catalyst'] if entry in pred and pred[entry] != '']
        conditions[rxn] = [rxn_cond]
        
        i += 1 

        if i % 100 == 0: 
            with open(cond_path, 'w') as f: 
                json.dump(conditions, f, indent='\t')


# save conditions 
with open(cond_path, 'w') as f: 
    json.dump(conditions, f, indent='\t')
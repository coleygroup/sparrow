""" Scores likelihood of reaction success using ASKCOSv2 WLDN5 forward predictor model """

import json 
import copy 
import time 
import requests 

from pathlib import Path 
from tqdm import tqdm 

HOSTS = [""] #, ""] # add host address(es)
PORT = "" # add port, assumes same port for all addresses
batch_sizes = [1] #, 8] 
host_batch = [i for h, c in zip(HOSTS, batch_sizes) for i in [h]*c]

# functions, templates, etc 
query_template = {
    "backend": "wldn5",
    "model_name": "pistachio",
    "smiles": [],
    "reagents": "",
    "solvent": ""
}

def post(host_post, params):
    req = requests.post(host_post, json=params, verify=False)
    try: 
        task_id = req.json()
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

# load reactions with conditions 
with open('examples/swanson/chkpts/conditions.json', 'r') as f: 
    conditions = json.load(f)

# load pre-evaluated reactions 
score_path = Path('examples/swanson/chkpts/scores.json')
if score_path.exists():
    with open(score_path, 'r') as f:
        scores = json.load(f)
else: 
    scores = {}

# assign conditions to not previously evaluated reactions 
rxns_eval = [smi for smi, cond in conditions.items() if smi not in scores]
step_size = sum(batch_sizes)
rxn_chunks = [rxns_eval[i:i+step_size] for i in range(0, len(rxns_eval), step_size)] 

i = 0 
for chunk in tqdm(rxn_chunks, position=0, leave=True):
    hosts = {}
    task_ids = {}

    # post for each node in the chunk 
    for host, rxn in zip(host_batch[:len(chunk)], chunk): 
        hosts[rxn] = host  
        params = copy.deepcopy(query_template)
        reacs, prod = rxn.split('>>')
        params["smiles"] = [reacs]
        params["reagents"] = '.'.join(conditions[rxn][0])
        if '[OH2--]' in params["reagents"] or 'F[Br](F)F' in params["reagents"]: # forward predictor doesn't like this 
            params["reagents"] = params["reagents"].replace('[OH2--].', '')
            params["reagents"] = params["reagents"].replace('.F[Br](F)F', '')
        task_ids[rxn] = post(host_post=f'http://{hosts[rxn]}:{PORT}/api/forward/controller/call-async/', params=params)
        
    # get and update nodes 
    for rxn in chunk: 
        _, prod = rxn.split('>>')
        result = get(
            host_results=f'http://{hosts[rxn]}:{PORT}/api/legacy/celery/',
            task_id=task_ids[rxn],
            sleep_time=0.1,
            timeout=30,
        )
        
        try: 
            output = result['output']
            outcomes = {entry["outcome"]: entry["prob"] for entry in output["result"][0]}
        except Exception as e:
            print(f'Error with reaction {rxn} on host {hosts[rxn]}, not assigning score')
            print(f'Error: {e}')
            print(f'Result: {result}')
            continue 
            # print(f'Error with reaction {node["smiles"]}, assigning score of 1e-6')
            # outcomes = {}

        if prod in outcomes: 
            scores[rxn] = outcomes[prod]
        else: 
            scores[rxn] = 0

        i += 1

        if i % 100 == 0: 
            with open(score_path, 'w') as f: 
                json.dump(scores, f, indent='\t')

# save scores 
with open(score_path, 'w') as f: 
    json.dump(scores, f, indent='\t')
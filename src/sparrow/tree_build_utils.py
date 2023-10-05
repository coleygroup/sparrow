""" Utilities to build trees using ASKCOS API, adapted from EASIE (https://github.com/itai-levin/easie) """

import requests
import time
from pprint import pprint
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import json
import urllib3
from joblib import Parallel, delayed
from tqdm import tqdm

urllib3.disable_warnings()

def post_and_get(HOST, params, sleep_time = 10, timeout = 650):
    req = requests.post(HOST+'/api/v2/tree-builder/', data=params, verify=False)
    print ('Request submitted')
    print (req.json())
    try:
        task_id = req.json()['task_id']
        print ('Retrieved task id:', task_id)
    except KeyError:
        return {} , {}
    results = requests.get(HOST + '/api/v2/celery/task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(HOST + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return req.json(), results.json()

        
def get_paths(smiles_ls, host, params, filename):
    results = {smiles:{} for smiles in smiles_ls}
    for smiles in tqdm(smiles_ls):
        params['smiles'] = smiles
        try:
            request, result = post_and_get(host, params)
            results[smiles] = result
        except ConnectionError:
            print ("Connection Error from " + host)
    
        with open(filename,'w') as f: 
            json.dump(results, f, indent='\t')

    return results 


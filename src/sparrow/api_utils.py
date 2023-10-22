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
from pathlib import Path

urllib3.disable_warnings()

def post_and_get(host_post, host_results, params, sleep_time = 10, timeout = 650):
    req = requests.post(host_post, data=params, verify=False)
    try:
        task_id = req.json()['task_id']
    except KeyError:
        return {} , {}
    results = requests.get(host_results + 'task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(host_results + 'task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return req.json(), results.json()


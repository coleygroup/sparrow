from typing import List, Dict
from pathlib import Path 
from abc import ABC, abstractmethod, abstractproperty
from tqdm import tqdm 
import requests
import time
import json
import urllib3
from pathlib import Path
from typing import List, Union
from sparrow.utils.json_utils import storage_from_api_response, save_storage_dict

urllib3.disable_warnings()

class PathFinder(ABC):
    """ 
    Determines whether a molecule represented as a SMILES is 
    buyable. If it is, determines its cost 
    TODO: implement LookupCoster so user can input csv of inventory with cost=0
    or maybe allow ChemSpaceCoster to be updated with inventory list and 
    change those buyables to cost = 0
    """

    @abstractmethod
    def get_save_trees() -> Path: 
        """ Stores the retrosynthesis tree object as a single json file, 
        returns the path to the file """


class AskcosAPIPlanner(PathFinder): 
    """ Uses ASKCOS API to get paths to target smiles """
    def __init__(self, 
                 host: str,
                 output_dir: Path,
                 time_per_target: int = 30, 
                 max_ppg: int = 10000,
                 max_branching: int = 20,
                 timeout: int = 600,
                 ):
        
        self.host = host
        self.time_per_target = time_per_target
        self.max_ppg = max_ppg
        self.max_branching = max_branching
        self.output_dir = output_dir
        self.timeout = timeout

        self.params = {
            'buyable_logic': 'or',
            'max_depth': '10',
            'expansion_time': str(self.time_per_target),
            'max_ppg': str(self.max_ppg),
            'return_first': 'false', 
            'max_branching': str(self.max_branching),
        }

    def get_save_trees(self, targets: List[str]): 
        """ Uses ASKCOS MCTS tree to generate paths to each smiles in 
        targets, combines the trees, and stores the output. returns the path
        of the combined trees """

        tree_path = self.output_dir/'combined_tree.json'
        path_ls = self.get_trees(targets, store_dir=self.output_dir/'askcos_trees')
        storage = self.combine_trees(path_ls)
        save_storage_dict(storage, filename=tree_path)

        print(f'Saving combined retrosynthesis tree to {str(tree_path)}')
        return tree_path
    
    def combine_trees(self, path_ls: List[Path]):

        storage = None
        for p in tqdm(path_ls, desc='Combining ASKCOS outputs'): 
            with open(p,'r') as f: 
                entry = json.load(f)
            for smi, path in entry.items(): 
                if 'output' in path and len(path['output'])>0: 
                    storage = storage_from_api_response({"result": path}, storage)

        return storage 
    
    def post_and_get(self, params, sleep_time = 10, timeout = 650):
        req = requests.post(self.host+'/api/v2/tree-builder/', data=params, verify=False)
        try:
            task_id = req.json()['task_id']
        except KeyError:
            return {} , {}
        results = requests.get(self.host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        clock = 0
        while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
            time.sleep(sleep_time)
            results = requests.get(self.host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
            clock += sleep_time
        return req.json(), results.json()

    def get_trees(self, smiles_ls: list, store_dir: Path) -> List[Path]: 

        Path(store_dir).mkdir(parents=True, exist_ok=True)
        
        for i, smiles in tqdm(enumerate(smiles_ls), total=len(smiles_ls), desc='Performing tree search'):
            results = {smiles:{} }
            self.params['smiles'] = smiles
            try:
                request, result = self.post_and_get(params=self.params, timeout=self.timeout)
                results[smiles] = result
            except ConnectionError:
                print ("Connection Error from " + self.host)
        
            with open(store_dir/f'tree_{i}.json','w') as f: 
                json.dump(results, f, indent='\t')

        return list(store_dir.glob('tree*.json')) 
    

class LookupPlanner(PathFinder):
    """ Loads in a json file that is in a retrosynthesis tree structure """
    def __init__(self, output_dir: Union[str, Path] = None, file_list: list = [], json_dir: Union[str, Path] = None) -> None:
        self.file_list = [Path(file_list) for file in file_list] 
        
        if json_dir is not None: 
            [self.file_list.append(p) for p in json_dir.glob('*.json')]
        
        self.tree_path = output_dir/'combined_tree.json'

    def combine_trees(self):
        storage = None
        for p in tqdm(self.file_list, desc='Combining ASKCOS outputs'): 
            with open(p,'r') as f: 
                entry = json.load(f)
            for smi, path in entry.items(): 
                if 'output' in path and len(path['output'])>0: 
                    storage = storage_from_api_response({"result": path}, storage)

        return storage 
        
    def get_save_trees(self, targets=None) -> Path:
        storage = self.combine_trees()        
        save_storage_dict(storage, filename=self.tree_path)
        return self.tree_path

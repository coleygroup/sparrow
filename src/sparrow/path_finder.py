from typing import List, Dict
from pathlib import Path
from abc import ABC, abstractmethod, abstractproperty
from tqdm import tqdm
import os
import requests
import time
import json
import urllib3
from pathlib import Path
from typing import List, Union
from sparrow.utils.json_utils import storage_from_api_response, save_storage_dict
from sparrow.utils.mcts_params import MCTS_PARAMS
import pdb

urllib3.disable_warnings()


#                  time_per_target: int = 30,
#                  max_ppg: int = 10000,
#                  max_branching: int = 20,
#                  timeout: int = 600,

# self.params = {
#             'buyable_logic': 'or',
#             'max_depth': '10',
#             'expansion_time': str(self.time_per_target),
#             'max_ppg': str(self.max_ppg),
#             'return_first': 'false',
#             'max_branching': str(self.max_branching),
# }

# RXN_QUERY_CONFIG = {
#   # "smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1",
#   "description": "",
#   "tags": "",
#   "expand_one_options": {
#     "template_count": 100,
#     "max_cum_template_prob": 0.995,
#     "forbidden_molecules": [],
#     "known_bad_reactions": [],
#     "retro_backend_options": [
#       {
#         "retro_backend": "template_relevance",
#         "retro_model_name": "reaxys",
#         "max_num_templates": 1000,
#         "max_cum_prob": 0.995,
#         "attribute_filter": []
#       }
#     ],
#     "use_fast_filter": "true",
#     "filter_threshold": 0.75,
#     "retro_rerank_backend": "relevance_heuristic",
#     "cluster_precursors": "false",
#     "cluster_setting": {
#       "feature": "original",
#       "cluster_method": "hdbscan",
#       "fp_type": "morgan",
#       "fp_length": 512,
#       "fp_radius": 1,
#       "classification_threshold": 0.2
#     },
#     "extract_template": "false",
#     "return_reacting_atoms": "true",
#     "selectivity_check": "false"
#   },
#   "build_tree_options": {
#     "expansion_time": 10,
#     "max_branching": 25,
#     "max_depth": 5,
#     "exploration_weight": 1,
#     "return_first": "false",
#     "max_trees": 500,
#     "buyable_logic": "and",
#     "max_ppg_logic": "none",
#     "max_ppg": 0,
#     "max_scscore_logic": "none",
#     "max_scscore": 0,
#     "chemical_property_logic": "none",
#     "max_chemprop_c": 0,
#     "max_chemprop_n": 0,
#     "max_chemprop_o": 0,
#     "max_chemprop_h": 0,
#     "chemical_popularity_logic": "none",
#     "min_chempop_reactants": 5,
#     "min_chempop_products": 5,
#     "custom_buyables": []
#   },
#   "enumerate_paths_options": {
#     "path_format": "json",
#     "json_format": "nodelink",
#     "sorting_metric": "plausibility",
#     "validate_paths": "true",
#     "score_trees": "false",
#     "cluster_trees": "false",
#     "cluster_method": "hdbscan",
#     "min_samples": 5,
#     "min_cluster_size": 5,
#     "paths_only": "false",
#     "max_paths": 200
#   },
#   "run_async": "false",
#   "result_id": "22016c62-271d-4958-9ec7-048b6c5fdce3"
# }

BUILD_TREE_OPTIONS = {
    # "smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1",
    "version": 2,
    "max_depth": 10,
    "max_branching": 20,
    "expansion_time": 10,
    "template_count": 500, # TODO: check it
    "max_cum_prob": 0.995,
    "max_chemicals": 500,
    "max_reactions": 500,
    # "max_iterations": 500,
    "max_templates": 1000,
    "buyable_logic": "or",
    "max_ppg_logic": "none",
    "max_ppg": 10000,

}


# BUILD_TREE_OPTIONS = {
#     # "expansion_time": 10,
#     "max_branching": 20,
#     "max_depth": 10,
#     "exploration_weight": 1,
#     "return_first": "false",
#     "max_trees": 500,
#     "buyable_logic": "or",
#     "max_ppg_logic": "none",
#     "max_ppg": 10000,
#     "max_scscore_logic": "none",
#     "max_scscore": 0,
#     "chemical_property_logic": "none",
#     "max_chemprop_c": 0,
#     "max_chemprop_n": 0,
#     "max_chemprop_o": 0,
#     "max_chemprop_h": 0,
#     "chemical_popularity_logic": "none",
#     "min_chempop_reactants": 5,
#     "min_chempop_products": 5,
#     "custom_buyables": []
#   }


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
                 time_per_target: int = 180,
                 max_ppg: int = 10000,
                 max_branching: int = 20,
                 timeout: int = 3000,
                 ):

        self.host = host
        self.time_per_target = time_per_target
        self.max_ppg = max_ppg
        self.max_branching = max_branching
        self.output_dir = output_dir
        self.timeout = timeout

        self.params = MCTS_PARAMS
        self.params["build_tree_options"]["expansion_time"] = str(self.time_per_target)
        self.params["build_tree_options"]["max_branching"] = str(self.max_branching)
        self.params["build_tree_options"]["max_ppg"] = str(self.max_ppg)
        self.params["build_tree_options"]["return_first"] = "false"
        # self.params = {
        #     'buyable_logic': 'or',
        #     'max_depth': '10',
        #     'expansion_time': str(self.time_per_target),
        #     'max_ppg': str(self.max_ppg),
        #     'return_first': 'false',
        #     'max_branching': str(self.max_branching),
        # }

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
        # req = requests.post(self.host+'/api/v2/tree-builder/', data=params, verify=False)
        req = requests.post(self.host+'/api/tree-search/mcts/call-sync-without-token/', json=params, verify=False)

        # try:
        #     task_id = req.json()['task_id']
        # except KeyError:
        #     return {} , {}
        # results = requests.get(self.host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        # clock = 0
        # while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        #     time.sleep(sleep_time)
        #     results = requests.get(self.host + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        #     clock += sleep_time
        # return req.json(), results.json()
        return req.json()

    def get_trees(self, smiles_ls: list, store_dir: Path) -> List[Path]:

        Path(store_dir).mkdir(parents=True, exist_ok=True)

        for i, smiles in tqdm(enumerate(smiles_ls), total=len(smiles_ls), desc='Performing tree search'):
            results = {smiles:{} }
            self.params['smiles'] = smiles
            print(self.params)
            # pdb.set_trace()

            try:
                # result = self.post_and_get(params=self.params, timeout=self.timeout,sleep_time=10)
                req = requests.post(self.host + '/api/tree-search/mcts/call-sync-without-token', json=self.params, verify=False)
                print(req)

                result = req.json()
                results[smiles] = result
            except ConnectionError:
                print ("Connection Error from " + self.host)

            # store_path = os.path.join(store_dir, f'tree_{i}.json')
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

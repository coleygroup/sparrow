from abc import ABC, abstractmethod
from typing import List, Union, Dict 
from pathlib import Path  
import json 
import requests
import time
from sparrow.utils.api_utils import post_and_get

def clean_context(context):
    """Clean a single context tuple from v1 condition recommender.

    Output is intended to be passed to forward evaluator.
    Modified from ASKCOS
    """
    temperature, solvent, reagent, catalyst = context

    # Remove chemicals without parsable smiles
    solvent = remove_names(solvent)
    reagent = remove_names(reagent)
    catalyst = remove_names(catalyst)

    # Remove any trailing periods (seems unnecessary?)
    solvent = solvent.rstrip(".")
    reagent = reagent.rstrip(".")
    catalyst = catalyst.rstrip(".")

    # Keep first solvent only
    if "." in solvent:
        solvent = solvent.split(".")[0]

    # Return as a list of chemicals only
    return [c for c in [solvent, reagent, catalyst] if c]


def remove_names(smiles):
    """Remove chemicals without parsable SMILES from the input"""
    smiles_list = smiles.split(".")
    smiles_list = [smi for smi in smiles_list if "Reaxys" not in smi]
    return ".".join(smiles_list)


class Recommender(ABC): 
    """ 
    Recommends reaction conditions (i.e. contexts) from a reaction smiles.
    """

    @abstractmethod
    def recommend_conditions(rxn_smi: str, n_c: int) -> List: 
        """ Returns a list of n_c recommended conditions (each a list). """

    def __call__(self, rxn_smi: str, n_c: int = 1) -> List: 
        return self.recommend_conditions(rxn_smi, n_c)

class AskcosLookupRecommender(Recommender): 
    """ Gets conditions from the output of the Askcos API. Optionally use a cleaner to clean the contexts """
    def __init__(self,
                 lookup_file: str, 
                 context_cleaner = None 
                 ): 
        
        self.cleaner = context_cleaner

        if Path(lookup_file).is_dir(): 
            self.data = self.explore_dir(Path(lookup_file))
        else: 
            self.data = self.explore_file(lookup_file)
        
    def explore_dir(self, lookup_dir: Path) -> Dict: 
        data = {}
        for file in lookup_dir.glob('*.json'): 
            data.update(self.explore_file(file))
        
        return data
    
    def explore_file(self, lookup_file: Union[str, Path]) -> Dict:
        with open(lookup_file, 'r') as f: 
            data = json.load(f)

        contexts = {}
        for rxnsmi, entry in data.items():
            contexts[rxnsmi] = []
            for context in entry['output']:
                contexts[rxnsmi].append([
                    context['temperature'], 
                    context['solvent'], 
                    context['reagent'], 
                    context['catalyst']]
                )
        
        if self.cleaner: 
            contexts = self.clean_contexts(contexts)
        
        return contexts
    
    def clean_contexts(self, context_dict): 
        clean_contexts = {
            rxn: [self.cleaner.clean_context(context) for context in value ]
            for rxn, value in context_dict.items()
        }
        
        return clean_contexts
    
    def recommend_conditions(self, rxn_smi, n_c) -> List: 
        contexts = self.data[rxn_smi]
        return contexts[:n_c]
        

class AskcosAPIRecommender(Recommender): 
    """ 
    Uses the ASKCOS v2 API to recommend contexts for reactions. 
    For each reaction, first tries the graph model. If that fails, uses the fingerprint model. 
    """
    
    def __init__(self, host: str, port: str = None, model: str = 'graph'): 

        if port: 
            self.host = f'{host}:{port}'
        else: 
            self.host = host 
        
        self.max_attempts = 3

    def process_result(self, result):
        condition = []
        for entry in result['agents']: 
            if entry.get('role')=='REACTANT': 
                continue 
            condition.append(entry['smi_or_name'])
        
        return condition    
            
    def recommend_conditions(self, rxn_smi: str, n_c: int, attempt=0, model='GRAPH') -> List: 
        """ Returns a list of n_c recommended conditions (each a list). """
        data = {
            "smiles": rxn_smi,
            "reagents": [],
            "n_conditions": 1
        }
        try:
            resp = requests.post(
                url=f"http://{self.host}/api/context-recommender/v2/condition/{model}/call-sync",
                json=data
            ).json()
        except requests.exceptions.JSONDecodeError:
            if attempt < self.max_attempts: 
                return self.recommend_conditions(rxn_smi, n_c, attempt=attempt+1, model='FP')
            else:     
                print(f'ASKCOS error recommending conditions for {rxn_smi}, returning empty conditions')
                return [[]]
        
        if len(resp['result']) == 0: 
            if attempt < self.max_attempts: 
                return self.recommend_conditions(rxn_smi, n_c, attempt=attempt+1, model='FP')
            else: 
                print(f'ASKCOS error recommending conditions for {rxn_smi}, returning empty conditions')
                return [[]]
        
        return [self.process_result(resp["result"][0])]
    

class AskcosV1APIRecommender(Recommender): 
    """ Uses the ASKCOS API to recommend contexts for reactions """
    
    def __init__(self, host: str): 
        self.host = host 

    def recommend_conditions(self, rxn_smi: str, n_c: int) -> List:
        contexts = self.get_conditions(rxn_smi, n_c)
        contexts = [self.clean_context(ctxt) for ctxt in contexts]
        return contexts
    
    def get_conditions(self, rxn_smi: str, n_c: int = 1) -> List: 
        """ Returns a list of n_c recommended conditions (each a list). """
        reactant_smis, product_smis = rxn_smi.split('>>')
        params = {
            'reactants': reactant_smis, 
            'products': product_smis, 
            'num_results': n_c, # default is 10
            'return_scores': 'true' # default is false
        }
        try:
            request, result = post_and_get(
                host_post=self.host+'/api/context/', 
                host_results=self.host+'/api/v2/celery/',
                params=params,
                sleep_time=0.1,
            )
        except ConnectionError:
            print ("Connection Error from " + self.host)
        
        if 'output' in result: 
            contexts = [
                [
                context['temperature'],
                context['solvent'],
                context['reagent'],
                context['catalyst']
                ]
                for context in result['output']
            ]
            return contexts
        else: 
            print(f'Context not recommended successfully for reaction {rxn_smi}')
            print(result)
            return [[]]
    
    def clean_context(self, context): 
        context = clean_context(context)
        return context 
    
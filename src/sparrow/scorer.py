""" A Scorer gives one or more scores to a reaction """
from abc import ABC, abstractmethod, abstractproperty
from typing import Union, List
from pathlib import Path 
import requests
from rdkit import Chem

from sparrow.api_utils import post_and_get


def canonicalize(smiles: str) -> str: 
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False)


class Scorer(ABC): 
    """ 
    A Scorer scores a reaction. The score should indicate
    the likelihood that a reaction will be successful. These scores are maximized 
    in the route optimization workflow. 

    """
    @abstractmethod
    def score_rxn(rxn_smi: str, condition: List = None) -> float:
        """ Calculates a reaction score """

    def __call__(self, rxn_smi, condition: List = None) -> float: 
        score = self.score_rxn(rxn_smi, condition)
        return score
    

class AskcosScorer(Scorer): 
    """ 
    Scores reactions using local ASKCOS 
    """
    def __init__(self, askcos_path = None):

        if askcos_path is not None: 
            import sys
            sys.path.append(askcos_path) 
        
        from askcos.synthetic.evaluation.evaluator import Evaluator

        self.requires_contexts = True 
        self.scorer = Evaluator()
    
    def score_rxn(self, rxn_smi: str, condition: List = None) -> float: 
        """ Scores a single reaction with a single set of conditions ('condition') """
        reactant_smi, product_smi = rxn_smi.split('>>')

        evaluation_results = self.scorer.evaluate(
            reactant_smiles=reactant_smi, 
            target=product_smi, 
            contexts=condition,
        )

        score = evaluation_results[0]['target']['prob']

        return score 


class AskcosAPIScorer(Scorer): 
    """ 
    Scores reactions by called ASKCOS API 
    """
    def __init__(self, host: str): 
        self.requires_contexts = True
        self.host = host 
    
    def score_rxn(self, rxn_smi: str, condition: List = None) -> float:
        if type(condition[0]) is list: 
            condition = condition[0]
        reactants, product = rxn_smi.split('>>')

        params = {
            'reactants': reactants, 
            'reagents': '.'.join(condition),
            'num_results': 10, 
        }
        try:
            _, result = post_and_get(
                host_post=self.host+'/api/forward/', 
                host_results=self.host+'/api/v2/celery/',
                params=params,
                sleep_time=1,
                timeout=60,
            )
        except ConnectionError:
            print ("Connection Error from " + self.host)

        outcomes = {
            prod['smiles']: prod['prob']
            for prod in result['output']
        }
        if canonicalize(product) in outcomes: 
            return outcomes[product]
        else: 
            return 0 
        

class LookupScorer(Scorer): 
    """ 
    Scores reactions by looking them up in a stored file 
    TODO: implement this 
    """
    def __init__(self, 
                 lookup_file: Union[str, Path], 
                 descriptor_type: str, 
                 ):

        self.requires_contexts = False 

        return 

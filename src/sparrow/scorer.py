""" A Scorer gives one or more scores to a reaction """
from abc import ABC, abstractmethod, abstractproperty
from typing import Union, List
from pathlib import Path 



class Scorer(ABC): 
    """ 
    A Scorer scores a reaction. The score should indicate
    the likelihood that a reaction will be successful. These scores are maximized 
    in the route optimization workflow. 

    """
    @abstractmethod
    def score_rxn(rxn_smi: str, condition: List = None) -> Union[float, int]:
        """ Calculates a reaction score """

    def __call__(self, rxn_smi, condition: List = None) -> float: 
        self.score_rxn(rxn_smi, condition)
    

class AskcosScorer(Scorer): 
    """ 
    Scores reactions using ASKCOS 
    """
    def __init__(self, 
                 context_recommender = None, 
                 ): 
        
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

    

class LookupScorer(Scorer): 
    """ Scores reactions by looking them up in a stored file """
    def __init__(self, 
                 lookup_file: Union[str, Path], 
                 descriptor_type: str, 
                 ):

        self.requires_contexts = False 

        return 

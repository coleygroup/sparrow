from abc import ABC, abstractmethod, abstractproperty
from typing import List, Union, Dict 
from pathlib import Path  
import json 

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
        

class AskcosRecommender(): 
    """ Uses ASKCOS's context recommender to predict conditions """
    
    def __init__(self): 

        from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
        import askcos.global_config as gc
        import askcos.utilities.contexts as context_cleaner
        
        self.recommender = NeuralNetContextRecommender()
        self.recommender.load_nn_model(
                model_path=gc.NEURALNET_CONTEXT_REC['model_path'], 
                info_path=gc.NEURALNET_CONTEXT_REC['info_path'], 
                weights_path=gc.NEURALNET_CONTEXT_REC['weights_path']
            )
        
        self.cleaner = context_cleaner
        
        return 
    
    def recommend_conditions(self, rxn_smi: str, n_c: int = 1) -> List: 
        
        conditions = self.get_conditions(rxn_smi, n_c)
        clean_condition = [self.clean_condition(cond) for cond in conditions]

        return clean_condition

    def get_conditions(self, rxn_smi: str, n_c: int = 1) -> List: 
        conditions = self.recommender.get_n_conditions(
            rxn_smi, 
            n=n_c, 
            return_separate=True, 
            return_scores=False
        )

        return conditions

    def clean_condition(self, condition: List) -> List:
        """ 
        Takes an 'uncleaned' condition that has temperatures and empty entries
        and returns a 'cleaned' context with only chemical 
        """
        condition = self.cleaner.clean_context(condition[:4])
        return condition

    def __call__(self, rxn_smi: str, n_c: int = 1): 
        
        return self.recommend_conditions(rxn_smi, n_c)



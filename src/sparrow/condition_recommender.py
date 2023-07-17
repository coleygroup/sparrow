from abc import ABC, abstractmethod, abstractproperty
from typing import List 

class Recommender(ABC): 
    """ 
    Recommends reaction conditions (i.e. contexts) from a reaction smiles.
    """

    @abstractmethod
    def recommend_conditions(rxn_smi: str, n_c: int) -> List: 
        """ Returns a list of n_c recommended conditions (each a list). """


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
        clean_conditions = [self.clean_condition[cond] for cond in conditions]

        return clean_conditions
    
    def __call__(self, rxn_smi: str, n_c: int = 1) -> List: 

        return self.recommend_conditions(rxn_smi, n_c)
        

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



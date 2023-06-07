import sys 
sys.path.append('/home/jfromer/sparrow') # need to figure out a way to get rid of this, to load askcos 

from sparrow.node import Node
from typing import Iterable, Optional, List

import numpy as np 


class ReactionNode(Node): 
    """ 
    A ReactionNode is a node of a RouteGraph that is a reaction.
    It is defined by its reactants, products, and conditions 
    """
    def __init__(self, 
                 smiles: str,
                 parents: Optional[List[Node]] = None,
                 children: Optional[List[Node]] = None,
                 dummy: Optional[bool] = False,
                 **kwargs) -> None:
        
        super().__init__(smiles, parents, children, **kwargs)
        
        self.conditions = [] #TODO: add conditions 
        self.score = 0
        self.penalty = np.inf
        self.dummy = dummy

        return 
    
    def calc_conditions_and_score(
            self, 
            context_recommender = None,
            evaluator = None, 
        ): 

        import askcos.utilities.contexts as context_cleaner
        from askcos.synthetic.evaluation.evaluator import Evaluator
        from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
        import askcos.global_config as gc

        if context_recommender is None: 
            context_recommender = NeuralNetContextRecommender()
            context_recommender.load_nn_model(
                model_path=gc.NEURALNET_CONTEXT_REC['model_path'], 
                info_path=gc.NEURALNET_CONTEXT_REC['info_path'], 
                weights_path=gc.NEURALNET_CONTEXT_REC['weights_path']
            )

        if evaluator is None: 
            evaluator = Evaluator()

        self.conditions = self.predict_conditions(context_recommender)
        reactant_smi, product_smi = self.smiles.split('>>')
        evaluation_results = evaluator.evaluate(
            reactant_smiles=reactant_smi, 
            target=product_smi, 
            contexts=[self.conditions],
        )

        self.score = evaluation_results[0]['target']['prob']
        if self.score == 0:
            self.penalty = 20 ## Do we want to remove these reactions from the graph? 
        else:
            self.penalty = np.min([1/self.score, 20])

        return self.conditions, self.score

    def predict_conditions(self, context_recommender: NeuralNetContextRecommender):
        contexts = context_recommender.get_n_conditions(
            self.smiles, 
            n=1, 
            return_separate=False, 
            return_scores=False
        )

        conditions = context_cleaner.clean_context(contexts[0])

        return conditions
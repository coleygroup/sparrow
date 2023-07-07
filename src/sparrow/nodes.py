from typing import Dict, Iterable, Optional, List
import numpy as np 


class Node: 
    """
    A Node is one node of a RouteGraph
    """

    def __init__(self, 
                 smiles: str, 
                 parents: Optional[List] = None, 
                 children: Optional[List] = None,
                 **kwargs) -> None:
        self.smiles = smiles
        self.parents = {}
        if parents: 
            self.update(parents=parents)

        self.children = {}
        if children: 
            self.update(children=children)

        self.selected = 0
        self.id = 0

    def update(self, parents=None, children=None): 
        """ Adds children and parents to node """

        if parents:
            for parent in parents: 
                self.parents[parent.smiles] = parent 
        
        if children: 
            for child in children: 
                self.children[child.smiles] = child

        return 

    def to_dict(self): 
        """ Converts node information to dictionary for saving as json """
        node_info = {
            'smiles': self.smiles, 
            'parents': [parent.smiles for parent in self.parents],
            'children': [child.smiles for child in self.children],
            'selected': self.selected,
            'id': self.id, 
        }

        return node_info

class CompoundNode(Node):
    """
    A CompoundNode is a node of a RouteGraph that represents 
    one compound. It is defined by its smiles 
    """

    def __init__(self, 
                smiles: str,
                parents: Optional[List[Node]],
                children: Optional[List[Node]],
                **kwargs) -> None:
    
        super().__init__(smiles, parents, children, **kwargs)
        
        self.buyable = 0 #TODO: update CompoundNode to determine buyability 
        self.cost_per_g = 0 #TODO: update CompoundNode to determine cost if buyable  
        self.reward = 0

        self.is_target = 0
        self.is_intermediate = 0
        
        return

    def set_as_target(self):
        self.is_target = True 

    def set_reward(self, reward):
        self.reward = reward
    
    def set_as_buyable(self):
        self.buyable = 1

    def set_as_intermediate(self):
        self.is_intermediate = 1
    
    def to_dict(self):
        node_info = super().to_dict()
        cmd_info = {
            'buyable': self.buyable, 
            'cost_per_g': self.cost_per_g, 
            'reward': self.reward, 
            'is_target': self.is_target, 
            'is_intermediate': self.is_intermediate
        }
        return {**node_info, **cmd_info}
        
    

class ReactionNode(Node): 
    """ 
    A ReactionNode is a node of a RouteGraph that is a reaction.
    It is defined by its reactants, products, and conditions 
    """
    def __init__(self, 
                 smiles: str,
                 parents: Optional[List[Node]] = None,
                 children: Optional[List[Node]] = None,
                 score: Optional[float] = None, 
                 conditions: Optional[List] = None,
                 dummy: Optional[bool] = False,
                 **kwargs) -> None:
        
        super().__init__(smiles, parents, children, **kwargs)
        
        if self.conditions is not None: # so that later conditions can be added through json route graph file 
            self.conditions = conditions 
            self.condition_set = True 
        else:
            self.conditions = [] #TODO: add conditions 
            self.conditions_set = False
        
        if self.score is not None: 
            self.score = score 
            self.score_set = True 
        else: 
            self.score = 0
            self.score_set = False 

        self.penalty = np.inf

        self.dummy = dummy # if it is a dummy reaction (no reactants, produces a starting material)

        return 
    
    def update_score(self, score): 
        
        self.score = score
        self.score_set = True 
    
    def update_conditions(self, conditions): 
        """ 
        Currently assuming only one set of conditions will apply to this reaction
        TODO: extend to multiple condition options 
        """
        self.conditions = conditions
        self.condition_set = True  

    def to_dict(self):
        node_info = super().to_dict()
        rxn_info = {
            'conditions': self.conditions, 
            'score': self.score, 
            'penalty': self.penalty,
            'dummy': self.dummy, 
        }
        return {**node_info, **rxn_info}
        
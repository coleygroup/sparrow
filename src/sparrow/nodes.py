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
                 id: Optional[str] = None,
                 selected: Optional[int] = 0,
                ) -> None:
        self.smiles = smiles
        self.parents = {}
        if parents: 
            self.update(parents=parents)

        self.children = {}
        if children: 
            self.update(children=children)

        self.selected = selected
        self.id = id

    def update(self, 
               parents = None, 
               children = None, 
            ) -> None: 
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
            'parents': [parent.smiles for parent in self.parents.values()],
            'children': [child.smiles for child in self.children.values()],
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
                parents: Optional[List[Node]] = None,
                children: Optional[List[Node]] = None,
                buyable: Optional[bool] = None, 
                cost_per_g: Optional[float] = None, 
                reward: Optional[float] = None,
                is_target: Optional[int] = 0,  
                is_intermediate: Optional[int] = 0,
                **kwargs) -> None:
    
        super().__init__(smiles, parents, children, **kwargs)
        
        self.buyable = buyable #TODO: update CompoundNode to determine buyability 
        self.cost_per_g = cost_per_g #TODO: update CompoundNode to determine cost if buyable  
        self.reward = reward

        self.is_target = is_target
        self.is_intermediate = is_intermediate
        
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
    
    def update(self,
               parents: Optional[List[Node]] = None, 
               children: Optional[List[Node]] = None, 
               buyable: Optional[bool] = None, 
               cost_per_g: Optional[float] = None, 
               reward: Optional[float] = None,
               is_target: Optional[int] = None,  
               is_intermediate: Optional[int] = None,
               ) -> None:
        
        super().update(parents, children)

        if buyable: 
            self.buyable = buyable
        
        if cost_per_g: 
            self.cost_per_g = cost_per_g

        if reward: 
            self.reward = reward
        
        if is_target: 
            self.is_target = is_target
        
        if is_intermediate: 
            self.is_intermediate = is_intermediate

        return 

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
                 condition: Optional[List] = None,
                 dummy: Optional[bool] = False,
                 **kwargs) -> None:
        
        super().__init__(smiles, parents, children, **kwargs)
        
        if condition is not None: # so that later conditions can be added through json route graph file 
            self.condition = condition
            self.condition_set = True 
        else:
            self.condition = [] #TODO: add conditions 
            self.condition_set = False
        
        if score is not None: 
            self.score = score 
            self.score_set = True 
            if score is not 0: 
                self.penalty = 1/score
            else: 
                self.penalty = None
        else: 
            self.score = 0
            self.score_set = False 
            self.penalty = None

        self.dummy = dummy # if it is a dummy reaction (no reactants, produces a starting material)

        return 
    
    def update_score(self, score): 
        
        self.score = score
        self.score_set = True 
    
    def update_condition(self, condition): 
        """ 
        Currently assuming only one set of conditions will apply to this reaction
        TODO: extend to multiple condition options 
        """
        self.condition = condition
        self.condition_set = True  

    def to_dict(self):
        node_info = super().to_dict()
        rxn_info = {
            'condition': self.condition, 
            'condition_set': self.condition_set,
            'score': self.score, 
            'score_set': self.score_set,
            'dummy': self.dummy, 
        }
        return {**node_info, **rxn_info} 

    def update(self,
               parents: Optional[List[Node]] = None, 
               children: Optional[List[Node]] = None, 
               condition = None, 
               condition_set = None, 
               score = None, 
               score_set = None, 
               dummy = None,
               **kwargs) -> None: 
        
        super().update(parents, children)

        if condition: 
            self.update_condition(condition)
        
        if condition_set: 
            self.condition_set = condition_set
        
        if score: 
            self.update_score(score)
        
        if score_set: 
            self.score_set = score_set
        
        if dummy: 
            self.dummy = dummy 
        
        return 
        
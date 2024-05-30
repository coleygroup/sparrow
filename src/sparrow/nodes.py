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

        self.blocked = False 

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
                parents: List[Node] = None,
                children: List[Node] = None,
                buyable: bool = None, 
                cost_per_g: float = None, 
                cost_function: list = None,
                reward: float = None,
                is_target: bool = False,  
                is_intermediate: bool = False,
                cost_set: int = False,
                **kwargs) -> None:
    
        super().__init__(smiles, parents, children, **kwargs)
        
        self.cost_set = cost_set 
        self.reward = None 
        self.buyable = None 
        self.cost_per_g = None 
        self.is_target = None 
        self.is_intermediate = None 
        self.cost_set = False 
        self.cost_function = None

        self.update(
            buyable=buyable, 
            cost_per_g=cost_per_g,
            reward=reward,
            is_target=is_target,
            is_intermediate=is_intermediate,
            cost_function=cost_function
        )
        
        return

    def set_as_target(self):
        self.is_target = True 

    def set_reward(self, reward):
        self.reward = reward
    
    def set_as_buyable(self):
        self.buyable = True

    def set_as_intermediate(self):
        self.is_intermediate = True
    
    def to_dict(self):
        node_info = super().to_dict()
        cmd_info = {
            'buyable': self.buyable, 
            'cost_per_g': self.cost_per_g, 
            'reward': self.reward, 
            'is_target': self.is_target, 
            'is_intermediate': self.is_intermediate,
            'cost_set': self.cost_set,
            'cost_function': self.cost_function,
        }
        return {**node_info, **cmd_info}
    
    def update_cost(self, cost_per_g = None):
        if cost_per_g is not None: 
            self.cost_per_g = cost_per_g
            self.cost_set = True
        else: 
            self.cost_per_g = None
        return self.cost_per_g
    
    def update_cost_function(self, amounts_in_g, prices_in_usd):
        self.cost_function = [amounts_in_g, prices_in_usd]
    
    def update(self,
               parents: List[Node] = None, 
               children: List[Node] = None, 
               buyable: bool = None, 
               cost_per_g: float = None, 
               reward: float = None,
               is_target: bool = None,  
               is_intermediate: bool = None,
               cost_function: list = None,
               **kwargs,
               ) -> None:
        
        super().update(parents, children)

        if buyable is not None: 
            self.buyable = buyable
            if buyable is False:
                self.cost_per_g = None
        
        if cost_per_g is not None: 
            self.update_cost(cost_per_g)

        if reward is not None: 
            self.reward = reward
        
        if is_target is not None: 
            if isinstance(is_target, int): 
                self.is_target = {0: False, 1: True}[is_target]
            else: 
                self.is_target = is_target
        
        if is_intermediate is not None: 
            if isinstance(is_intermediate, int): 
                self.is_intermediate = {0: False, 1: True}[is_intermediate]
            else: 
                self.is_intermediate = is_intermediate
        
        if cost_function is not None: 
            self.update_cost_function(cost_function[0], cost_function[1])

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
                 max_penalty: Optional[float] = 20, 
                 penalty: Optional[float] = None,
                 **kwargs) -> None:
        
        super().__init__(smiles, parents, children, **kwargs)
        
        self.max_penalty = max_penalty
        self.score = None
        self.score_set = False
        self.penalty = None
        self.condition = None
        self.condition_set = False
        self.dummy = 0

        self.update(condition=condition, score=score, dummy=dummy, penalty=penalty)

        return 
    
    def update_score(self, score): 
        
        self.score = score
        self.score_set = True 
        if score == 0: 
            self.penalty = self.max_penalty
        elif self.dummy is True: 
            self.penalty = 0
        else: 
            self.penalty = min(self.max_penalty, 1/score)
    
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
               penalty = None, 
               **kwargs) -> None: 
        
        super().update(parents, children)

        if condition is not None: 
            self.update_condition(condition)
        
        if condition_set is not None: 
            self.condition_set = condition_set
        
        if score is not None: 
            self.update_score(score)
        
        if score_set is not None: 
            self.score_set = score_set
        
        if dummy is not None: 
            self.dummy = dummy 
        
        if penalty is not None: 
            self.penalty = penalty
        
        return 
    
    def get_condition(self, n_c: int) -> List[str]:
        if self.condition: 
            return self.condition[:n_c]
        else: 
            return [[] for _ in range(n_c)]
from sparrow.node import Node 
from typing import Iterable, Optional, List 


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

    
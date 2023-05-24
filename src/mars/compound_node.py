from src.mars.node import Node 
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

        return

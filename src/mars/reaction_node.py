from src.mars.node import Node
from typing import Iterable, Optional, List


class ReactionNode(Node): 
    """ 
    A ReactionNode is a node of a RouteGraph that is a reaction.
    It is defined by its reactants, products, and conditions 
    """
    def __init__(self, 
                 smiles: str,
                 parents: Optional[List[Node]],
                 children: Optional[List[Node]],
                 **kwargs) -> None:
        
        super().__init__(smiles, parents, children, **kwargs)

        return 

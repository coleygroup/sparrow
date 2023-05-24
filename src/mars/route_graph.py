from src.mars.reaction_node import ReactionNode
from src.mars.compound_node import CompoundNode

from askcos.retrosynthetic.mcts.tree_builder import MCTS
import askcos.utilities.contexts as context_cleaner
from askcos.synthetic.evaluation.evaluator import Evaluator
from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
import askcos.global_config as gc

from src.mars.node import Node
from typing import Iterable, Dict, Union, Optional, List
from pathlib import Path
from itertools import chain 

import pickle 


class RouteGraph: 
    """
    A RouteGraph is a directed graph consisting of reactions (ReactionNodes)
    and compounds (CompoundNodes).
    """
    def __init__(self, 
                 askcos_MCTS_tree: Optional[MCTS] = None, 
                 paths: List = None, 
                ) -> None:
        
        self.compound_nodes = {}
        self.reaction_nodes = {}
        
        if askcos_MCTS_tree or paths: 
            self.build_from_MCTS(MCTS_tree=askcos_MCTS_tree, paths=paths)

        
        return
    
    def add_reaction_node(self, 
                          smiles: str, 
                          parents: Optional[Iterable[str]] = None, 
                          children: Optional[Iterable[str]] = None,
                          **kwargs) -> ReactionNode:
        """ Adds reaction node to graph from reaction smiles. If parents/chidlren provided, first add them as 
        compound nodes  """ 

        if parents: 
            parent_nodes = [self.add_compound_node(parent) for parent in parents]
        else: 
            parent_nodes = None

        if children: 
            child_nodes = [self.add_compound_node(child) for child in children]
        else: 
            child_nodes = None
        
            
        if smiles in self.reaction_nodes.keys(): 
            self.reaction_nodes[smiles].update(
                parents=parent_nodes, 
                children=child_nodes, 
                **kwargs)
        else:
            self.reaction_nodes[smiles] = ReactionNode(
                smiles, 
                parents=parent_nodes, 
                children=child_nodes, 
                **kwargs)
        
        # go back to parent nodes and add this new node as a child 
        if parents: 
            for parent in parents: 
                self.compound_nodes[parent].update(children=[self.reaction_nodes[smiles]])

        # go back to child nodes and add this new node as a parent 
        if children: 
            for child in children: 
                self.compound_nodes[parent].update(parents=[self.reaction_nodes[smiles]])
        
        return self.reaction_nodes[smiles]
    
    def add_compound_node(self, 
                          smiles: str, 
                          parents: Optional[Iterable[str]] = None, 
                          children: Optional[Iterable[str]] = None,
                          **kwargs) -> CompoundNode: 
        """ Adds compound node to the graph from smiles. If parents/children are provided, first 
        adds those as reaction nodes. """
        
        if parents: 
            parent_nodes = [self.add_reaction_node(parent) for parent in parents]
        else: 
            parent_nodes = None

        if children: 
            child_nodes = [self.add_reaction_node(child) for child in children]
        else: 
            child_nodes = None
        
            
        if smiles in self.compound_nodes.keys(): 
            self.compound_nodes[smiles].update(
                parents=(parent_nodes or None), 
                children=(child_nodes or None), 
                **kwargs
            )
        else:
            self.compound_nodes[smiles] = CompoundNode(
                smiles, 
                parents=parent_nodes, 
                children=child_nodes, 
                **kwargs
            )
        
        # go back to parent nodes and add this new node as a child 
        if parents: 
            for parent in parents: 
                self.reaction_nodes[parent].update(children=[self.compound_nodes[smiles]])

        # go back to child nodes and add this new node as a parent 
        if children: 
            for child in children: 
                self.reaction_nodes[parent].update(parents=[self.compound_nodes[smiles]])

        return self.compound_nodes[smiles]
    
    def nodes(self) -> Dict[str, Node]:
        """ Returns a  dictionary [smi, ReactionNode] of all nodes in the graph """
        return {**self.compound_nodes, **self.reaction_nodes}
    
    def save(self, filename: Union[str, Path] = None) -> None: 
        """ Saves route graph using as pickle file """
        if filename is None: 
            filename = Path.cwd() / "route_graph.pickle"
        with open(filename, 'wb') as file:
            pickle.dump(self, file) 
        
        return 
    
    def add_path(self, path): 

        if len(path['children']) > 0: 
            for p in path['children']:
                self.add_path(p)
        
        parents = [p['smiles'] for p in path['children']] or None
        # ^ confusing because children to ASKCOS MCTS are parents here 

        if ">>" in path['smiles']: 
            self.add_reaction_node(smiles=path['smiles'], parents=parents)
        else: 
            self.add_compound_node(smiles=path['smiles'], parents=parents)

        return 

    def build_from_MCTS(self, MCTS_tree: MCTS = None, paths: List = None): 
        """ 
        Builds RouteGraph from Monte Carlo Tree Search output from 
        ASKCOS or paths output from MCTS.get_buyable_paths
        """
        if paths is None: 
            if MCTS_tree is None: 
                print("Must provide at least one of MCTS_tree or paths")
                return 
            
            paths = MCTS_tree.enumerate_paths()

        for path in paths: 
            self.add_path(path)

        return 

    def calc_reaction_scores(self, 
            context_recommender: NeuralNetContextRecommender = None,
            evaluator: Evaluator = None, 
        ):

        for rxn_node in self.reaction_nodes.values(): 
            rxn_node.calc_conditions_and_score(
                context_recommender,
                evaluator,
            )
        
        return


    
def load_route_graph(filename: Union[str, Path]) -> RouteGraph:
    """ Loads route graph from pickle file """
    with open(filename, 'rb') as file:
        route_graph = pickle.load(file)
    
    return route_graph

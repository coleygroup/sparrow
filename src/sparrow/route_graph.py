from sparrow import ReactionNode, CompoundNode, Node
from sparrow.scorer import Scorer, AskcosScorer

from typing import Iterable, Dict, Union, Optional, List
from pathlib import Path
from tqdm import tqdm
import sys, os 

import pickle 


class RouteGraph: 
    """
    A RouteGraph is a directed graph consisting of reactions (ReactionNodes)
    and compounds (CompoundNodes).
    """
    def __init__(self, 
                 askcos_MCTS_tree = None, 
                 paths: List = None, 
                ) -> None:
        
        self.compound_nodes = {}
        self.reaction_nodes = {}
        
        if askcos_MCTS_tree or paths:

            self.add_from_MCTS(MCTS_tree=askcos_MCTS_tree, paths=paths)

        self.ids = {}
        
        return
    
    def add_reaction_node(self, 
                          smiles: str, 
                          parents: Optional[Iterable[str]] = None, 
                          children: Optional[Iterable[str]] = None,
                          dummy: Optional[bool] = False,
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
                dummy=dummy,
                **kwargs)
        
        # go back to parent nodes and add this new node as a child 
        if parents: 
            for parent in parents: 
                self.compound_nodes[parent].update(children=[self.reaction_nodes[smiles]])

        # go back to child nodes and add this new node as a parent 
        if children: 
            for child in children: 
                self.compound_nodes[child].update(parents=[self.reaction_nodes[smiles]])
        
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
        """ TODO: insert description """
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

    def add_from_MCTS(self, MCTS_tree = None, paths: List = None): 
        """ 
        Builds RouteGraph from Monte Carlo Tree Search output from 
        ASKCOS or paths output from MCTS.get_buyable_paths
        """
        from askcos.retrosynthetic.mcts.tree_builder import MCTS
        import askcos.utilities.contexts as context_cleaner
        from askcos.synthetic.evaluation.evaluator import Evaluator
        from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
        import askcos.global_config as gc

        if paths is None: 
            if MCTS_tree is None: 
                print("Must provide at least one of MCTS_tree or paths")
                return 
            
            paths = MCTS_tree.enumerate_paths()

        for path in paths: 
            self.add_path(path)

        return 
    
    def set_targets(self, targets: Iterable[str]):
        """ defines CompoundNode.is_target = 1 for specified targets, where 
        targets is a list of smiles strings of target molecules """
        for target in targets: 
            if target not in self.compound_nodes: 
                print(f"Target {target} not tagged in route graph because it is not already a node")
            else: 
                self.compound_nodes[target].set_as_target()

        return 
    
    def set_target_rewards(self, target_dict: Dict[str, Union[int, float]]):
        """ sets CompoundNode.reward to specified value for each key in target_dict """
        for target, reward in target_dict.items(): 
            if target not in self.compound_nodes: 
                print(f"Reward not set for target {target} because it is not already a node")
            else:
                self.compound_nodes[target].set_reward(reward)
        
        return 
    
    def set_buyable_compounds(self):
        """ sets CompoundNode.buyable for starting materials listed in ****insert**** 
        (for now using arbitrary definition based on number of carbons, nitrogens, oxygens) """
        for smiles, node in self.compound_nodes.items():
            nc = smiles.count("c") + smiles.count("C")
            nn = smiles.count("N") + smiles.count("n")
            no = smiles.count("O") + smiles.count("o")
            if nc<=12 and nn<=3 and no<=5:
                node.set_as_buyable()

        return 

    def set_intermediates(self):
        """ Run after targets and starting materials are set """
        for node in self.compound_nodes.values():
            if (node.is_target==0) and (node.buyable==0):
                node.set_as_intermediate()
        
        return 

    def set_compound_types(self, target_dict):
        """ TODO: insert description """

        self.set_buyable_compounds()
        self.set_targets(target_dict.keys())
        self.set_target_rewards(target_dict)
        self.set_intermediates()

        return 
    
    def intermediate_nodes(self): 
        """ Returns a list of CompoundNodes that have is_intermediate==1 """
        inter_nodes = [node for node in self.compound_nodes.values() if node.is_intermediate==1 ]
        return inter_nodes
    
    def buyable_nodes(self):
        """ Returns a list of CompoundNodes that have buyable==1 """
        buyable_nodes = [node for node in self.compound_nodes.values() if node.buyable==1]
        return buyable_nodes
    
    def id_nodes(self):
        """ Sets IDs for all nodes """
        
        self.rxn_ids = {}
        for node, i in zip(self.reaction_nodes.values(), range(len(self.reaction_nodes))):
            self.rxn_ids[f"R{i}"] = node
            node.id = f"R{i}"

        self.compound_ids = {}
        for node, i in zip(self.compound_nodes.values(), range(len(self.compound_nodes))):
            self.compound_ids[f"C{i}"] = node
            node.id = f"C{i}"
        
        self.ids = {**self.rxn_ids, **self.compound_ids}

        return self.ids

    def reaction_nodes_only(self) -> List[ReactionNode]
        """ Returns a list of ReactionNodes in this graph """
        return self.reaction_nodes.values()



    
def load_route_graph(filename: Union[str, Path]) -> RouteGraph:
    """ Loads route graph from pickle file """
    with open(filename, 'rb') as file:
        route_graph = pickle.load(file)
    
    return route_graph

from sparrow import ReactionNode, CompoundNode, Node
from sparrow.scorer import Scorer, AskcosScorer
from sparrow.coster import Coster, NaiveCoster
from typing import Iterable, Dict, Union, Optional, List
from pathlib import Path
from tqdm import tqdm
import sys, os 
import json 

import pickle 


class RouteGraph: 
    """
    A RouteGraph is a directed graph consisting of reactions (ReactionNodes)
    and compounds (CompoundNodes).
    """
    def __init__(self, 
                 node_filename: str = None, 
                ) -> None:
        """
        node_filename: filename of a json with information about reaction and compound nodes 
        """
        
        self.compound_nodes = {}
        self.reaction_nodes = {}
        
        if node_filename:

            self.add_from_json(node_filename)

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
                self.reaction_nodes[child].update(parents=[self.compound_nodes[smiles]])

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
    
    def set_targets(self, targets: Iterable[str]):
        """ defines CompoundNode.is_target = 1 for specified targets, where 
        targets is a list of smiles strings of target molecules """
        invalid_targets = []
        for target in targets: 
            if target not in self.compound_nodes: 
                print(f"Target {target} not tagged in route graph because it is not already a node")
                print(f"Target will not be considered!!!!")
                invalid_targets.append(target)
            else: 
                self.compound_nodes[target].set_as_target()

        return invalid_targets
    
    def set_target_rewards(self, target_dict: Dict[str, Union[int, float]]):
        """ sets CompoundNode.reward to specified value for each key in target_dict """
        for target, reward in target_dict.items(): 
            if target not in self.compound_nodes: 
                print(f"Reward not set for target {target} because it is not already a node")
            else:
                self.compound_nodes[target].set_reward(reward)
        
        return 
    
    def set_buyable_compounds_and_costs(self, coster: Coster = None):
        """ sets CompoundNode.buyable for starting materials listed in ****insert**** 
        (for now using arbitrary definition based on number of carbons, nitrogens, oxygens) """
        if coster is None: 
            coster = NaiveCoster()
        
        for node in tqdm(self.compound_nodes_only(), 'Searching ChemSpace'): 
            buyable, cost = coster.get_buyable_and_cost(node.smiles)
            node.update(
                buyable=buyable, 
                cost_per_g=cost,
            )

        return 

    def set_intermediates(self):
        """ Run after targets and starting materials are set """
        for node in self.compound_nodes.values():
            if (node.is_target==0) and (node.buyable==0):
                node.set_as_intermediate()
        
        return 

    def set_compound_types(self, target_dict, coster=None):
        """ TODO: insert description """

        self.set_buyable_compounds_and_costs(coster)
        invalid_targets = self.set_targets(target_dict.keys())
        for invalid_tar in invalid_targets: 
            target_dict.pop(invalid_tar)
        self.set_target_rewards(target_dict)
        self.set_intermediates()

        return target_dict
    
    def intermediate_nodes(self) -> List[CompoundNode]: 
        """ Returns a list of CompoundNodes that have is_intermediate==1 """
        inter_nodes = [node for node in self.compound_nodes.values() if node.is_intermediate==1 ]
        return inter_nodes
    
    def buyable_nodes(self) -> List[CompoundNode]:
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

    def reaction_nodes_only(self) -> List[ReactionNode]:
        """ Returns a list of ReactionNodes in this graph """
        return list(self.reaction_nodes.values())

    def compound_nodes_only(self) -> List[CompoundNode]:
        """ Returns a list of ReactionNodes in this graph """
        return list(self.compound_nodes.values())

    def node_from_id(self, id) -> Node: 
        return self.ids[id]
    
    def node_from_smiles(self, smiles) -> Node: 
        return self.nodes()[smiles]
    
    def id_from_smiles(self, smiles) -> str: 
        return self.nodes()[smiles].id

    def smiles_from_id(self, id) -> str: 
        return self.node_from_id(id).smiles

    def to_json(self, filename) -> None: 
        compound_nodes = [node.to_dict() for node in self.compound_nodes_only()]
        reaction_nodes = [node.to_dict() for node in self.reaction_nodes_only()]
        storage = {
            'Compound Nodes': compound_nodes, 
            'Reaction Nodes': reaction_nodes
        }
        with open(filename, 'w') as f: 
            json.dump(storage, f, indent="\t")

        return 

    def add_from_json(self, filename) -> None:
        """ Adds node information from a json file """
        with open(filename, 'r') as f: 
            storage = json.load(f)
        
        if 'Compound Nodes' in storage.keys(): 
            for node_info in storage['Compound Nodes']: 
                self.add_compound_node(
                    smiles=node_info.pop('smiles'),
                    parents=node_info.pop('parents') if 'parents' in node_info else [],
                    children=node_info.pop('children') if 'children' in node_info else [],
                    **node_info,
                )
        
        if 'Reaction Nodes' in storage.keys(): 
            for node_info in storage['Reaction Nodes']: 
                self.add_reaction_node(
                    smiles=node_info.pop('smiles'),
                    parents=node_info.pop('parents') if 'parents' in node_info else [],
                    children=node_info.pop('children') if 'children' in node_info else [], 
                    **node_info,
                )
        
        return



    
def load_route_graph(filename: Union[str, Path]) -> RouteGraph:
    """ Loads route graph from pickle file """
    with open(filename, 'rb') as file:
        route_graph = pickle.load(file)
    
    return route_graph

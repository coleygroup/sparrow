from sparrow import ReactionNode, CompoundNode, Node
from sparrow.coster import Coster
from typing import Iterable, Dict, Union, Optional, List
from pathlib import Path
from tqdm import tqdm
from datetime import datetime
from scipy.sparse import csr_matrix
import json 
import numpy as np 
import pickle
import warnings  
import networkx as nx 

def find_cycles_nx(adjacency_matrix):
    G = nx.DiGraph(adjacency_matrix)
    cycles = list(nx.simple_cycles(G))
    # filtered_cycles = [cycle for cycle in cycles if len(cycle) <= k]

    return cycles

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
                          dummy: Optional[bool] = None,
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
                dummy=dummy,
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
    
    def remove_dummy_rxns(self) -> None: 
        for node in self.reaction_nodes_only(): 
            if node.smiles.startswith('>>'): 
                self.remove_rxn_node(node.smiles)
        return 
    
    def prune_dummy_rxns(self) -> None: 
        for node in self.dummy_nodes_only(): 
            child = list(node.children.values())[0]
            if not child.buyable: 
                self.remove_rxn_node(node.smiles)
        return         
    
    def remove_rxn_node(self, smi) -> None: 
        self.reaction_nodes.pop(smi)

        for node in self.compound_nodes_only(): 
            if smi in node.parents: 
                node.parents.pop(smi)
            if smi in node.children: 
                node.parents.pop(smi)

        return 
    
    def remove_compound_node(self, smi) -> None: 
        self.compound_nodes.pop(smi)

        for node in self.reaction_nodes_only(): 
            if smi in node.parents: 
                node.parents.pop(smi)
            if smi in node.children: 
                node.parents.pop(smi)

        return         
    
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

        for target in targets: 
            self.node_from_id(target).set_as_target()

        return
    
    def set_target_rewards(self, target_dict: Dict[str, Union[int, float]]):
        """ sets CompoundNode.reward to specified value for each key in target_dict """
        for target_id, reward in target_dict.items(): 
            if self.smiles_from_id(target_id) not in self.compound_nodes: 
                print(f"Reward not set for target {self.smiles_from_id(target_dict)} because it is not already a node")
            else:
                self.node_from_id(target_id).set_reward(reward)
        
        return 
    
    def set_buyable_compounds_and_costs(self, coster: Coster = None, save_json_dir: str = None):
        """ sets CompoundNode.buyable and CompoundNode.cost_per_g for starting materials in ChemSpace """

        if coster is None: 
            return 
        
        coster.build_status_log()
        
        prog_bar = tqdm(total=len(self.compound_nodes_only()), desc= 'Searching for price/buyability', position=0)
        c = 0
        for node in self.compound_nodes_only(): 
            if node.cost_set and node.cost_per_g != 1: 
                prog_bar.update(1)
                continue 

            coster.status_log.set_description_str(f'Searching for {node.smiles}')

            buyable, cost = coster.get_buyable_and_cost(node.smiles)
            node.update(
                buyable=buyable, 
                cost_per_g=cost,
            )
            c += 1
            prog_bar.update(1)

            if c % 500 == 0 and save_json_dir is not None: 
                time = datetime.now().strftime("%H-%M-%S")
                self.to_json(Path(save_json_dir) / f'trees_w_costs_{time}.json')
            
        coster.status_log.set_description_str('Saving retrosynthetic graph with costs')
        self.to_json(Path(save_json_dir) / f'trees_w_costs.json')
        coster.status_log.set_description_str('Done searching for buyability and costs')

        return 

    def set_intermediates(self):
        """ Run after targets and starting materials are set """
        for node in self.compound_nodes.values():
            if (not node.is_target) and (not node.buyable):
                node.set_as_intermediate()
            else: 
                node.update(is_intermediate=False)
        
        return 

    def set_compound_types(self, target_dict, save_dir: str = None, coster=None):
        """ TODO: insert description """

        self.set_targets(target_dict.keys())
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
        for i, node in enumerate(self.reaction_nodes_only()):
            self.rxn_ids[f"R{i}"] = node
            node.id = f"R{i}"

        self.compound_ids = {}
        for i, node in enumerate(self.compound_nodes_only()):
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

    def dummy_nodes_only(self) -> List[ReactionNode]: 
        return [node for node in self.reaction_nodes_only() if node.dummy]
    
    def non_dummy_nodes(self) -> List[ReactionNode]: 
        return [node for node in self.reaction_nodes_only() if not node.dummy]
    
    def child_of_dummy(self, dummy_node_id: str) -> str: 
        """ Returns the ID corresponding to the starting material node that is produced by the dummy reaction """
        dummy_node = self.node_from_id(dummy_node_id)
        child_node = list(dummy_node.children.values())[0]
        return child_node.id 

    def node_from_id(self, id) -> Node: 
        return self.ids[id]
    
    def node_from_smiles(self, smiles) -> Node: 
        return self.nodes()[smiles]
    
    def id_from_smiles(self, smiles) -> str: 
        return self.nodes()[smiles].id

    def smiles_from_id(self, id) -> str: 
        return self.node_from_id(id).smiles

    def unique_reagents(self) -> List[str]: 
        conditions = [node.get_condition(1)[0] for node in self.reaction_nodes_only()]
        reagents = set([i for cond in conditions for i in cond])
        return reagents    

    def compute_adjacency_matrix(self) -> csr_matrix: 
        id_to_ind = {node.id: ind for ind, node in enumerate(self.nodes().values())}

        rows = []
        cols = []

        for idd, ind in id_to_ind.items():
            child_ids = [child.id for child in self.node_from_id(idd).children.values()]
            for id in child_ids: 
                rows.append(id_to_ind[id])
                cols.append(ind)
        
        data = (np.ones(len(rows)), (np.array(rows), np.array(cols)) )
         
        A = csr_matrix( data, shape=(len(id_to_ind), len(id_to_ind)) )        

        return A, id_to_ind

    def dfs_find_cycles_nx(self) -> list: 
        A, id_to_ind = self.compute_adjacency_matrix()
        cycles = find_cycles_nx(A)
        ind_to_id = {ind: idd for idd, ind in id_to_ind.items()}
        cycless = [[ind_to_id[ind] for ind in cyc if ind_to_id[ind].startswith('R') ] for cyc in cycles]
        return cycless 

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

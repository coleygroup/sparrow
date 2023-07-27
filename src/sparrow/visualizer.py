""" Visualizes selected routes """
from sparrow.route_graph import RouteGraph
from sparrow import ReactionNode
from pulp import LpVariable
from typing import List, Dict
import networkx as nx
import json 
import pandas as pd 
import matplotlib.pyplot as plt


color_dict = {
    "Starting": "green",
    "Intermediate": "blue",
    "Target": "red",
    "None": "gray"
}

class Visualizer: 
    """ Description of Visualizer class """
    def __init__(self,
                 route_graph: RouteGraph,
                 nonzero_vars: List[LpVariable],
                 ) -> None:
        self.route_graph = route_graph
        self.vars = nonzero_vars
        self.var_names = [var.name for var in nonzero_vars]

        self.reaction_nodes = set()
        self.intermediate_nodes = set()
        self.starting_nodes = set()
        self.target_nodes = set()

        self.extract_vars(nonzero_vars)

        self.target_verts, self.inter_verts, self.start_verts, self.rxn_verts = self.create_vertices()
        self.vertices = {**self.target_verts, **self.inter_verts, **self.start_verts, **self.rxn_verts}
        self.edges = self.create_edges()

        self.vis_graph = self.create_vis_graph()

    def reactions_from_vars(self, filename): 
        rxnfortar = [var.name for var in self.vars if var.name.startswith('rxnfortarget')]
        layers = self.multipartite_from_digraph()
        target_ids = list(self.target_verts.keys())
        storage = {}
        for tar_id in target_ids:
            rxns = [rxn.split('_')[1] for rxn in rxnfortar if rxn.endswith(tar_id)]
            tar = self.route_graph.smiles_from_id(tar_id)
            storage[tar] = [self.route_graph.smiles_from_id(rxn_id) for rxn_id in rxns]
            # for rxn_id in rxns: 
            #     c = 0
            #     if self.route_graph.node_from_id(rxn_id).smiles.startswith('>>'):
            #         storage[tar][f'starting_material_{c}'] = self.route_graph.smiles_from_id(rxn_id)
            #         c+=1
            #     else: 
            #         storage[tar][layers[rxn_id]] = [self.route_graph.smiles_from_id(rxn_id), self.route_graph.node_from_id(rxn_id).get_condition(1)[0]]

        with open(filename,'w') as f: 
            json.dump(storage, f, indent="\t")



    def extract_vars(self, vars: List[LpVariable]) -> None:

        for var in vars: 
            
            if 'rxnflow' in var.name: 
                id = var.name.split('_')[1]
                node = self.route_graph.node_from_id(id)
                self.reaction_nodes.add(node)
                if not node.dummy: 
                    self.add_compounds_from_rxn(node)
            
            if var.name.startswith('target'):
                id = var.name.split('_')[1]
                node = self.route_graph.ids[id]
                self.target_nodes.add(node)
                
        return 

    def add_compounds_from_rxn(self, rxn_node: ReactionNode) -> None:

        for child in rxn_node.children.values(): 
            if child.is_target: 
                continue 

            self.intermediate_nodes.add(child)

        for parent in rxn_node.parents.values(): 
            # if parent.is_target: 
            #     continue 

            # figure out if the parent is a dummy starting reaction 
            dummy_parent_selected = False
            for par in parent.parents.values(): 
                if f"rxnflow_{par.id}" in self.var_names and self.route_graph.ids[par.id].dummy:
                    dummy_parent_selected = True 

            if dummy_parent_selected: 
                self.starting_nodes.add(parent)

            else:
                self.intermediate_nodes.add(parent)


        return 
    
    def create_vertices(self) -> Dict[str, Dict]: 

        target_verts = { node.id:{
                    "ID": node.id,
                    "Node Type": "Compound",
                    "Compound Class": "Target", 
                }
            for node in self.target_nodes}
        
        inter_verts = { node.id:{
                    "Node Type": "Compound",
                    "Compound Class": "Intermediate", 
                }
            for node in self.intermediate_nodes}
        
        start_verts = { node.id:{
                    "Node Type": "Compound",
                    "Compound Class": "Starting", 
                }
            for node in self.starting_nodes}
        
        rxn_verts = { node.id:{
                    "Node Type": "Reaction",
                    "Compound Class": "None", 
                }
            for node in self.reaction_nodes if not node.dummy}
        
        return target_verts, inter_verts, start_verts, rxn_verts
    
    def create_edges(self) -> List:
        
        edges = []
        for rxn_node in self.reaction_nodes: 
            if not rxn_node.dummy: 
                for parent in rxn_node.parents.values(): 
                    edges.append([parent.id, rxn_node.id])
                for child in rxn_node.children.values(): 
                    edges.append([rxn_node.id, child.id])

        return edges 
    
    def create_vis_graph(self): 

        vis_graph = nx.DiGraph(directed=True)
        vis_graph.add_edges_from(self.edges)

        return vis_graph
    
    def multipartite_from_digraph(self): 
        layers = {}
        for id in self.vis_graph.nodes: 
            if id in self.target_verts: 
                lay = 0
            else: 
                lay = self.max_distance_to_target(id)
            layers[id] = lay
        
        nx.set_node_attributes(self.vis_graph, values = layers, name='layer')

        return layers

    def max_distance_to_target(self, node_id): 
        dists = []
        for target in self.target_verts: 
            try: 
                long_path = max(nx.all_simple_paths(self.vis_graph, node_id, target), key=lambda x: len(x))
            except: 
                continue 
            dists.append(len(long_path)-1)
        
        return max(dists)
    
    def plot_graph(self, path: str = 'debug/graph_fig.png'): 
        
        self.multipartite_from_digraph()

        # pos = nx.spring_layout(self.vis_graph)
        pos = nx.multipartite_layout(self.vis_graph, subset_key = "layer")   

        options = {
            'nodelist': list(self.vis_graph),
            'node_color': [color_dict[self.vertices[id]["Compound Class"]] for id in list(self.vis_graph)],
            'node_size': 600,
            'width': 3,
            'arrowstyle': '-|>',
            'arrowsize': 12,
        }

        plt.figure(1, figsize=(30,30)) 


        # targets 
        nx.draw_networkx_nodes(self.vis_graph, pos, 
                               nodelist=self.target_verts, 
                               node_size=1000, 
                               node_color=color_dict["Target"],
                               label='Target')
        
        
        # intermediates
        nx.draw_networkx_nodes(self.vis_graph, pos,
                               nodelist=self.inter_verts, 
                               node_size=1000, 
                               node_color=color_dict["Intermediate"],
                               label='Intermediate')
        
        # starting
        nx.draw_networkx_nodes(self.vis_graph, pos,
                               nodelist=self.start_verts, 
                               node_size=1000, 
                               node_color=color_dict["Starting"],
                               label='Starting Material')
        
        # reactions
        nx.draw_networkx_nodes(self.vis_graph, pos,
                               nodelist=self.rxn_verts, 
                               node_size=100, 
                               node_color=color_dict["None"],
                               label='Reaction')
        
        # edges 
        nx.draw_networkx_edges(self.vis_graph, pos, width=3, arrowstyle='-|>', arrowsize=12, node_size=1000)
        # add nodelist and nodesize to this so that arrows are not hidden 

        nx.draw_networkx_labels(self.vis_graph, pos)

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        lgnd = plt.legend(by_label.values(), by_label.keys())

        for handle in lgnd.legendHandles: 
            handle._sizes = [30]

        plt.savefig(path)
   
        prefix = path.split('.png')[0]

        ids = list(self.vis_graph.nodes)
        smis = [self.route_graph.ids[id].smiles for id in ids]
        df = pd.DataFrame({"ID": list(ids), "SMILES": smis})
        df = df.set_index('ID').sort_values('ID')
        df.to_csv(f'{prefix}_key.csv')

        return 







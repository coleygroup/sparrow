""" Visualizes selected routes """
from sparrow.route_graph import RouteGraph
from sparrow.node import Node
from sparrow.reaction_node import ReactionNode
from pulp import LpVariable
from typing import List, Dict
import pandas as pd 
import networkx as nx
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

        self.plot_graph()

        print('done')



    def extract_vars(self, vars: List[LpVariable]) -> None:

        for var in vars: 
            
            if 'rxnflow' in var.name: 
                id = var.name.split('_')[1]
                node = self.route_graph.ids[id]
                self.reaction_nodes.add(node)
                if not node.dummy: 
                    self.add_compounds_from_rxn(node)
            
            if var.name.startswith('target'):
                id = var.name.split('_')[1]
                node = self.route_graph.ids[id]
                self.target_nodes.add(node)
                
        return 

    def add_compounds_from_rxn(self, rxn_node: ReactionNode) -> None:
        compounds = [*(rxn_node.parents.values()), *(rxn_node.children.values())]

        for child in rxn_node.children.values(): 
            if child.is_target: 
                continue 

            self.intermediate_nodes.add(child)

        for parent in rxn_node.parents.values(): 
            if parent.is_target: 
                continue 

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
    
    def plot_graph(self, path: str = 'debug/graph_fig.png'): 

        pos = nx.nx_pydot.graphviz_layout(self.vis_graph)
                              
        options = {
            'nodelist': list(self.vis_graph),
            'node_color': [color_dict[self.vertices[id]["Compound Class"]] for id in list(self.vis_graph)],
            'node_size': 600,
            'width': 3,
            'arrowstyle': '-|>',
            'arrowsize': 12,
        }

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

        lgnd = plt.legend()
        for handle in lgnd.legendHandles: 
            handle._sizes = [30]

        plt.savefig(path)


        
        return 







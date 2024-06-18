import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import json 
import pandas as pd 
import numpy as np 
from pathlib import Path 
import networkx as nx 
from sparrow.route_graph import RouteGraph
from tqdm import tqdm
import random 
import colorsys

theme_colors = ['#405E90','#D45127','#AA2165','#818084']

cluster_colors = {
    'None': '#0E713E',
    'fps': '#44AA99',
    'objs': '#CC6677',
    'both': '#882255'
}

it_colors = ['#3491C1', '#7D2AC1', '#B9305C', '#DC5501', '#DE9A00', '#377501', '#B4B5B4']

def set_style():
    """set_style"""
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["font.family"] = "sans-serif"
    sns.set(context="paper", style="ticks") 
    mpl.rcParams["text.color"] = "black"
    mpl.rcParams["axes.labelcolor"] = "black"
    mpl.rcParams["axes.edgecolor"] = "black"
    mpl.rcParams["axes.labelcolor"] = "black"
    mpl.rcParams["xtick.color"] = "black"
    mpl.rcParams["ytick.color"] = "black"
    mpl.rcParams["xtick.major.size"] = 2.5
    mpl.rcParams["ytick.major.size"] = 2.5

    mpl.rcParams["xtick.major.width"] = 0.45
    mpl.rcParams["ytick.major.width"] = 0.45

    mpl.rcParams["axes.edgecolor"] = "black"
    mpl.rcParams["axes.linewidth"] = 0.45
    mpl.rcParams["font.size"] = 8
    mpl.rcParams["axes.labelsize"] = 8
    mpl.rcParams["axes.titlesize"] = 8
    mpl.rcParams["figure.titlesize"] = 8
    mpl.rcParams["figure.titlesize"] = 8
    mpl.rcParams["legend.fontsize"] = 7
    mpl.rcParams["legend.title_fontsize"] = 8
    mpl.rcParams["xtick.labelsize"] = 7
    mpl.rcParams["ytick.labelsize"] = 7
    mpl.rcParams['figure.dpi'] = 300

    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['legend.fancybox'] = False
    mpl.rcParams['legend.facecolor'] = "none"

    mpl.rcParams['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth

def set_size(w, h, ax=None):
    """w, h: width, height in inches
    Resize the axis to have exactly these dimensions
    """
    if not ax:
        ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w) / (r - l)
    figh = float(h) / (t - b)
    ax.figure.set_size_inches(figw, figh)

def calc_expected_reward(route_file: Path):

    with open(route_file, 'r') as f: 
        routes = json.load(f)

    expected_reward = 0
    for entry in routes.values(): 
        rew = entry['Reward']
        rxn_scores = [rentry['score'] for rentry in entry['Reactions'] if 'score' in rentry]
        expected_reward += rew*np.prod(rxn_scores)

    return expected_reward

def df_from_dir_void(results_dir):
    result_dir = Path(results_dir)
    summaries = []
    for dir in sorted(list(result_dir.glob('lam*'))): 
        with open(dir/'summary.json', 'r') as f:
            summ = json.load(f)
        
        l1, l2, l3 = summ.pop('Weights')
        summ['Utility Weight'] = l1
        summ['Starting Material Weight'] = l2
        summ['Reaction Weight'] = l3
        summ['Expected Reward'] = calc_expected_reward(dir / 'routes.json')
        summaries.append(summ)

    df = pd.DataFrame.from_dict(summaries)
    return df

def df_from_dir(results_dir):
    result_dir = Path(results_dir)
    summaries = []
    for dir in sorted(list(result_dir.glob('*'))): 
        with open(dir / 'solution' / 'summary.json', 'r') as f:
            summ = json.load(f)
        
        # l1, l2, l3 = summ.pop('Weights')
        # summ['Utility Weight'] = l1
        # summ['Starting Material Weight'] = l2
        # summ['Reaction Weight'] = l3
        summ['Expected Reward'] = calc_expected_reward(dir / 'solution' / 'routes.json')
        summaries.append(summ)

    df = pd.DataFrame.from_dict(summaries)
    return df

def graph_vis(tree_path: str, routes_path: str, cleaned_tar_path: str, percent_plot=0.1): 
    with open(routes_path,'r') as f: 
        routes = json.load(f)
    
    tar_df = pd.read_csv(cleaned_tar_path)
    
    target_smis = list(tar_df['SMILES'])
    
    route_graph = RouteGraph(node_filename=tree_path)
    route_graph.id_nodes()
    _, id_to_ind = route_graph.compute_adjacency_matrix()

    sel_rxns = set([item['smiles'] for entry in routes.values() for item in entry['Reactions']])
    sel_tars = list(routes.keys())
    sel_cpds = set([*sel_tars, *[item for entry in routes.values() for item in entry['Compounds']] ])

    G = nx.DiGraph(directed=True)
    sel_ids = []
    tar_ids = []
    cpd_ids = []
    for cpd_node in route_graph.compound_nodes_only(): 
        G.add_node(id_to_ind[cpd_node.id])
        cpd_ids.append(id_to_ind[cpd_node.id])
        if cpd_node.smiles in sel_cpds: 
            sel_ids.append(id_to_ind[cpd_node.id])
        if cpd_node.smiles in sel_tars: 
            tar_ids.append(id_to_ind[cpd_node.id])
    
    rxn_ids = []
    start_ids = []
    for rxn_node in route_graph.reaction_nodes_only(): 
        G.add_node(id_to_ind[rxn_node.id])
        rxn_ids.append(id_to_ind[rxn_node.id])
        if rxn_node.smiles in sel_rxns: 
            sel_ids.append(id_to_ind[rxn_node.id])
            if rxn_node.smiles.startswith('>>'):
                start_ids.append(id_to_ind[rxn_node.id])

    sel_edges = []
    for rxn_node in route_graph.reaction_nodes_only(): 
        if random.uniform(0,1) > percent_plot and rxn_node.smiles not in sel_rxns: 
            continue 
        for parent in rxn_node.parents.values(): 
            G.add_edge(id_to_ind[parent.id], id_to_ind[rxn_node.id])
            if rxn_node.smiles in sel_rxns: 
                sel_edges.append((id_to_ind[parent.id], id_to_ind[rxn_node.id]))
        for child in rxn_node.children.values(): 
            G.add_edge(id_to_ind[rxn_node.id], id_to_ind[child.id])
            if rxn_node.smiles in sel_rxns: 
                sel_edges.append((id_to_ind[rxn_node.id], id_to_ind[child.id]))

    layers = {}
    target_verts = [id_to_ind[route_graph.id_from_smiles(smi)] for smi in target_smis]

    for node_ind in tqdm(range(len(id_to_ind))): 
        if node_ind in target_verts: 
            layers[node_ind] = 0
        else: 
            try: 
                lay = calc_layer(G, node_ind, target_verts)
                layers[node_ind] = lay
            except: 
                # layers[node_ind] = lay
                G.remove_node(node_ind)
    
    nx.set_node_attributes(G, values = layers, name='layer')

    fig1, ax1 = plt.subplots()
    # pos = nx.multipartite_layout(G, subset_key = "layer")
    pos = nx.spring_layout(G, )
    nodelist2 = [node for node in G.nodes() if node in sel_ids]
    nodecolor = []
    for n in nodelist2:
        if n in tar_ids: 
            nodecolor.append(theme_colors[2])
        elif n in cpd_ids:
            nodecolor.append(theme_colors[3])
        elif n in start_ids: 
            nodecolor.append(theme_colors[0])
        else: 
            nodecolor.append(theme_colors[1])
    # nodecolor = [theme_colors[0] if n in cpd_ids elif n in else theme_colors[2] for n in nodelist2]
    edgelist2 = [edge for edge in G.edges() if edge in sel_edges]
    ns = 20
    nx.draw_networkx_edges(G, pos, edgelist=G.edges(), edge_color='lightgray', node_size=ns, width=0.5, arrowsize=5, ax=ax1)
    nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_color='lightgray', node_size=ns, ax=ax1, edgecolors='lightgray', linewidths=0.1)

    fig2, ax2 = plt.subplots()
    nx.draw_networkx_edges(G, pos, edgelist=edgelist2, edge_color='dimgrey', node_size=ns, width=0.5, arrowsize=5, ax=ax2)
    nx.draw_networkx_nodes(G, pos, nodelist=nodelist2, node_color=nodecolor, node_size=ns, ax=ax2, edgecolors='k', linewidths=0.1)

    
    return fig1, ax1, fig2, ax2
    
def calc_layer(G: nx.DiGraph, node_ind: int, target_verts: list): 
    ps = list(nx.all_simple_paths(G, node_ind, target_verts, cutoff=10))    
    return len(min(ps, key=lambda x: len(x)))

def df_to_latex(df: pd.DataFrame, path: str): 
    df = df.rename(columns={"Utility Weight": "$\lambda_1$", "Starting Material Weight": "$\lambda_2$", "Reaction Weight": "$\lambda_3$"})
    df = df.rename(columns={"Total reward": "Cumulative utility", "Cost starting materials": "Cost of starting materials", \
                            "Number reaction steps": "Number of reaction steps"})
    cols = ['$\lambda_1$','$\lambda_2$','$\lambda_3$',"Cumulative utility", "Cost of starting materials", "Number of reaction steps", "Average reaction score" ]
    df = df[cols]
    df = df.sort_values(by=['$\lambda_1$','$\lambda_2$', '$\lambda_3$'])
    df['$\lambda_1$'] = [f'{a:0.2f}' for a in df['$\lambda_1$']]
    df['$\lambda_2$'] = [f'{a:0.2f}' for a in df['$\lambda_2$']]
    df['$\lambda_3$'] = [f'{a:0.2f}' for a in df['$\lambda_3$']]

    df = df.round({
        # '$\lambda_1$':1,'$\lambda_2$':1,'$\lambda_3$':3, 
        "Cumulative utility":1, "Cost of starting materials": 1,"Number of reaction steps":0,"Average reaction score": 2})
    str_df = df.to_latex(escape=False, index=False, multicolumn_format='c')
    with open(path,'w') as f:
        f.writelines(str_df)
    return df

def df_to_latex_baselines(df: pd.DataFrame, path: str): 
    df = df.round({
        # '$\lambda_1$':1,'$\lambda_2$':1,'$\lambda_3$':3, 
        'Cumulative Reward':1, 'Starting Material Cost': 1,'Number of Reactions':0,"Average Reaction Score": 2})
    str_df = df.to_latex(escape=False, index=False, multicolumn_format='c')
    with open(path,'w') as f:
        f.writelines(str_df)
    return df

def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

def make_color_darker(scale, color: str): 
    rgb = mpl.colors.ColorConverter.to_rgb(color)
    return scale_lightness(rgb, scale)
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import json 
import pandas as pd 
from pathlib import Path 

method_colors = {
    'ei':'#405E90',
    'pi':'#D45127',
    'nds':'#AA2165',
    'random':'#818084',
    'greedy': '#AA2165', 
}


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

def df_from_dir(results_dir):
    result_dir = Path(results_dir)
    summaries = []
    for dir in sorted(list(result_dir.glob('lam*'))): 
        with open(dir/'summary.json', 'r') as f:
            summ = json.load(f)
        
        l1, l2, l3 = summ.pop('Weights')
        summ['Reward Weight'] = l1
        summ['Starting Material Weight'] = l2
        summ['Reaction Weight'] = l3

        summaries.append(summ)

    df = pd.DataFrame.from_dict(summaries)
    return df
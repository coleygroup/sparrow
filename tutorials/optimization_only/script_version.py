import sys, os 
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located
os.environ["CUDA_VISIBLE_DEVICES"]="-1"

import pandas as pd 

target_dict = {
    "COC5=CC(N=CN=C6NC7=CC=CC(OC)=C7)=C6C=C5OC" : 3,
    "COC1=CC(N=CN=C2NC3=CC=C(OCC4=CC=CC=C4)C=C3)=C2C=C1OC" : 4,
    "COC8=CC(N=CN=C9NC%10=CC=CC(Cl)=C%10)=C9C=C8OC" : 2,
    "COC%11=CC(N=CN=C%12NC%13=CC=CC(O)=C%13)=C%12C=C%11OC" : 11,
    "COC%14=CC(N=CN=C%15NC%16=CC=C(NC(C%17=CC=CC=C%17)=O)C=C%16)=C%15C=C%14OC" : 14,
}


from sparrow.tree_build_utils import build_rxn_graph
route_graph = build_rxn_graph(
    target_smis=list(target_dict.keys()),
    n_cpus=10,
    time_per_target=15
)

from sparrow.route_selector import RouteSelector
from sparrow.route_graph import load_route_graph

# route_graph = load_route_graph('tutorials/with_askcos/route_graph_w_scores.pickle')

route_selector = RouteSelector(
    route_graph,
    target_dict,
    calc_reaction_scores=True,
)

# route_selector = RouteSelector(
#     route_graph,
#     target_dict,
#     calc_reaction_scores=False,
# )

route_selector.define_variables()
route_selector.set_objective()
route_selector.set_constraints()
route_selector.optimize()

from sparrow.visualizer import Visualizer
import matplotlib.pyplot as plt

vis = Visualizer(
    route_graph,
    nonzero_vars=route_selector.optimal_variables(),
    )

fig_out_dir = "tutorials/with_askcos/tree_vis.png"
vis.plot_graph(path=fig_out_dir)

plt.show()


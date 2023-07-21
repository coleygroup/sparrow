import sys, os 
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located
sys.path.append('/home/jfromer/sparrow')
os.environ["CUDA_VISIBLE_DEVICES"]="-1" # use askcos on CPUs only 

from sparrow.coster import ChemSpaceCoster, NaiveCoster
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
# from sparrow.tree_build_utils import build_rxn_graph
from sparrow.visualizer import Visualizer
import matplotlib.pyplot as plt
from keys import chemspace_api_key

# smis = ["CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C", "CC(C)(C)OC(=O)N1CCCCCC1C(O)=O", "CNC(C)(C)C#C", "BrCc1cccc(Br)c1", "Fc1cc(CBr)ccc1Br",]

api_key = "dZH7vZYK2JDKWxgMSCKIBQZcKfteL395UuYtCuHoVk1WUcpq1MIeiPn95mBLsXOh"

# coster = ChemSpaceCoster(api_key=api_key)
# costs = coster(smis)
target_dict = {
    "COC5=CC(N=CN=C6NC7=CC=CC(OC)=C7)=C6C=C5OC" : 3,
    "COC1=CC(N=CN=C2NC3=CC=C(OCC4=CC=CC=C4)C=C3)=C2C=C1OC" : 4,
    "COC8=CC(N=CN=C9NC%10=CC=CC(Cl)=C%10)=C9C=C8OC" : 2,
    "COC%11=CC(N=CN=C%12NC%13=CC=CC(O)=C%13)=C%12C=C%11OC" : 11,
    "COC%14=CC(N=CN=C%15NC%16=CC=C(NC(C%17=CC=CC=C%17)=O)C=C%16)=C%15C=C%14OC" : 14,
}

# build_rxn_graph(
#     target_smis=list(target_dict.keys()),
#     n_cpus=10,
#     time_per_target=15,
#     filename='debug/askcos_paths.json'
# )

route_graph = RouteGraph()
route_graph.add_from_json('debug/askcos_paths_w_scores.json')

route_selector = RouteSelector(
    route_graph,
    target_dict,
    condition_recommender=None,
    rxn_scorer=None,
    coster=ChemSpaceCoster(api_key=chemspace_api_key),
    weights=[1,1,1,1], 
) # conditions and scoring will happen during initialization with the current setup 

route_selector.graph.to_json('debug/askcos_paths_v2.json')

route_selector.define_variables()
route_selector.set_objective()
route_selector.set_constraints()
route_selector.optimize(solver=None) # solver='GUROBI' for GUROBI (license needed)

vis = Visualizer(
    route_graph,
    nonzero_vars=route_selector.optimal_variables(),
    )

fig_out_dir = "debug/optimal_routes_ChemSpace.png"
vis.plot_graph(fig_out_dir)

plt.show()
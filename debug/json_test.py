""" To see if json export works for RouteGraph """
import sys, os 
sys.path.append('/home/jfromer/sparrow/askcos-core') # change this for your system to where askcos folder is located
os.environ["CUDA_VISIBLE_DEVICES"]="-1" # use askcos on CPUs only 

from sparrow.tree_build_utils import build_rxn_graph
from sparrow.route_graph import RouteGraph

def make_tree_and_graph(): 

    target_dict = {
        'O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1': 10
    }

    build_rxn_graph(
        target_smis=list(target_dict.keys()),
        n_cpus=10,
        time_per_target=15
    )

    return 

def read_json(filename='debug/json_test.json'):
    graph = RouteGraph()
    graph.add_from_json(filename)
    print('done')

# make_tree_and_graph()
read_json('debug/askcos_paths.json')
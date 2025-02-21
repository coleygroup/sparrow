import json
import pandas as pd
import networkx as nx 
from tqdm import tqdm 
from sparrow.coster import LookupCoster
from sparrow.route_graph import RouteGraph
from sparrow.utils.parallel_utils import chunked_parallel

graph = RouteGraph(node_filename='examples/automoldesigner/chkpts/graph_reactions_only.json')

# prune
cutoff_distance = 16
graph.id_nodes() 

nx_graph, id_to_ind = graph.create_nx_graph()
nx_graph = nx_graph.to_undirected()
graph_inds = [ind for _, ind in id_to_ind.items()] 
graph_nodeids = [nid for nid, _ in id_to_ind.items()] 

def keep_node_id(ind, id_to_ind, nx_graph, target_ids, cutoff_distance): 
    for target_id in target_ids: 
        if nx.has_path(nx_graph, source=ind, target=id_to_ind[target_id]):
            if nx.shortest_path_length(nx_graph, source=ind, target=id_to_ind[target_id]) < cutoff_distance:
                return True
    return False   

target_ids = [graph.compound_nodes[s].id for s in pd.read_csv('examples/automoldesigner/targets.csv').SMILES]
parallel_fn = lambda x: keep_node_id(x, id_to_ind, nx_graph, target_ids, cutoff_distance)
keep_results = chunked_parallel(input_list=graph_inds, chunks=1000, function=parallel_fn, max_cpu=100)

# keep_results = [parallel_fn(i) for i in graph_inds]
remove_ids = [idd for idd, res in zip(graph_nodeids, keep_results) if not res]

print(len(remove_ids))
print(f'Removing {len(remove_ids)/len(id_to_ind)*100:0.2f}% of nodes by pruning with cutoff distance {cutoff_distance}')

for nid in remove_ids: 
    if nid.startswith('C'): 
        graph.remove_compound_node(smi=graph.smiles_from_id(nid), remove_neighbors=False)
    else: 
        graph.remove_rxn_node(smi=graph.smiles_from_id(nid), remove_neighbors=False)

# set buyability 
coster = LookupCoster(
    lookup_file='examples/automoldesigner/data/enamine_per_g_122024.csv'
)
graph.set_buyable_compounds_and_costs(
    coster=coster,
    save_json_dir='examples/automoldesigner/chkpts',
    save_freq=1e10
)
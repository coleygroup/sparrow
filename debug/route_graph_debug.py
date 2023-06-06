from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.visualizer import Visualizer
import pickle 


# with open('debug/simple_tree.pickle', 'rb') as file:
    # paths = pickle.load(file)

# route_graph = RouteGraph(paths=paths)

with open('debug/RouteGraph_withScores.pickle', 'rb') as file: 
    route_graph = pickle.load(file)

target_dict = {
    'O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1': 10
}

route_selector = RouteSelector(
    route_graph,
    target_dict,
    calc_reaction_scores=False
)

route_selector.define_variables()
route_selector.set_objective()
route_selector.set_constraints()

route_selector.optimize()

nonzero_vars = route_selector.routes_from_optimal_variables()

visualizer = Visualizer(route_graph=route_selector.graph, nonzero_vars=nonzero_vars)

print('Optimization problem done')
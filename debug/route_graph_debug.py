from mars.route_graph import RouteGraph
from mars.route_selector import RouteSelector
import pickle 


with open('simple_tree.pickle', 'rb') as file:
    paths = pickle.load(file)

route_graph = RouteGraph(paths=paths)

target_dict = {
    'O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1': 100000
}

route_selector = RouteSelector(
    route_graph,
    target_dict,
)

route_selector.define_variables()
route_selector.set_objective()
route_selector.set_constraints()

route_selector.optimize()

print('Optimization problem done')
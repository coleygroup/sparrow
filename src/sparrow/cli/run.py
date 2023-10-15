from pathlib import Path
import pandas as pd 
import json 
import sys 

from sparrow.path_finder import AskcosAPIPlanner, LookupPlanner
from sparrow.route_graph import RouteGraph
from sparrow.route_selector import RouteSelector
from sparrow.condition_recommender import AskcosRecommender, AskcosAPIRecommender
from sparrow.scorer import AskcosScorer, AskcosAPIScorer
from sparrow.coster import ChemSpaceCoster, NaiveCoster
from sparrow.cli.args import get_args



def get_target_dict(filepath: Path): 
    df = pd.read_csv(filepath)
    target_dict = {
        smiles: reward
        for smiles, reward in zip(df['SMILES'], df['Reward'])
    }
    return target_dict, list(df['SMILES'])

def optimize(selector, params):

    selector.define_variables()
    selector.set_objective()
    selector.set_constraints()
    solver = {'pulp': None, 'gurobi': 'GUROBI'}[params['solver']]
    selector.optimize(solver=solver) # solver='GUROBI' for GUROBI (license needed)

    return selector 

def export_selected_nodes(selector: RouteSelector, rxn_list, starting_list, target_list, base_dir): 
    storage = {'Starting Materials': [], 'Reactions': [], 'Targets': []}
    graph = selector.graph

    for rxn_id in rxn_list:
        node = graph.node_from_id(rxn_id)
        if node.dummy: 
            storage['Reactions'].append({
                'smiles': node.smiles,
                'starting material cost ($/g)': selector.cost_of_dummy(node), 
            })
        else: 
            storage['Reactions'].append({
                'smiles': node.smiles,
                'conditions': node.get_condition(1)[0], 
                'score': node.score,
            })

    for cpd_id in starting_list: 
        node = graph.node_from_id(cpd_id)
        storage['Starting Materials'].append({
            'smiles': node.smiles,
            'cost': node.cost_per_g
        })

    for cpd_id in target_list: 
        node = graph.node_from_id(cpd_id)
        storage['Targets'].append({
            'smiles': node.smiles,
            'reward': node.reward,
        })
    
    with open(base_dir/'optimal_routes.json','w') as f:
        json.dump(storage, f, indent='\t')

    return storage 

def extract_vars(selector: RouteSelector, base_dir): 
    sys.setrecursionlimit(100000)
    nonzero_vars = [
            var for var in selector.problem.variables() if var.varValue > 0.01
        ]
    rxn_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('rxn')]
    mol_ids = [var.name.split('_')[1] for var in nonzero_vars if var.name.startswith('mol')]

    selected_targets = set(mol_ids) & set(selector.targets)
    starting_mats = set([node.id for node in selector.graph.buyable_nodes()])
    selected_starting = set(mol_ids) & starting_mats
    print(f'{len(selected_targets)} targets selected using {len(rxn_ids)} reactions and {len(selected_starting)} starting materials')
    export_selected_nodes(selector, rxn_ids, selected_starting, selected_targets, base_dir)

    storage = {}
    for target in selected_targets: 
        store_dict = {'Compounds':[], 'Reactions':[]}
        smi = selector.graph.smiles_from_id(target)
        storage[smi] = find_mol_parents(store_dict, target, mol_ids, rxn_ids, selector.graph)
        storage[smi]['Reward'] = selector.target_dict[target]

    with open(base_dir/f'routes_{len(selected_targets)}tars.json','w') as f: 
        json.dump(storage, f, indent='\t')

    return storage 

def find_rxn_parents(store_dict, rxn_id, selected_mols, selected_rxns, graph): 
    par_ids = [n.id for n in graph.node_from_id(rxn_id).parents.values()]
    selected_pars = set(par_ids) & set(selected_mols)
    for par in selected_pars: 
        store_dict['Compounds'].append(graph.smiles_from_id(par))
        store_dict = find_mol_parents(store_dict, par, selected_mols, selected_rxns, graph)
    return store_dict

def find_mol_parents(store_dict, mol_id, selected_mols, selected_rxns, graph): 
    par_ids = [n.id for n in graph.node_from_id(mol_id).parents.values()]
    selected_pars = set(par_ids) & set(selected_rxns)
    for par in selected_pars: 
        node = graph.node_from_id(par)
        if node.dummy: 
            store_dict['Reactions'].append({
                'smiles': node.smiles,
                'starting material cost ($/g)': selector.cost_of_dummy(node), 
            })
        else: 
            store_dict['Reactions'].append({
                'smiles': node.smiles,
                'conditions': node.get_condition(1)[0], 
                'score': node.score,
            })
        store_dict = find_rxn_parents(store_dict, par, selected_mols, selected_rxns, graph)
    return store_dict

def get_path_storage(params, targets): 
    if params['path_finder'] == 'lookup': 
        planner = LookupPlanner(params['tree_lookup'])
    elif params['path_finder'] == 'api': 
        planner = AskcosAPIPlanner( 
            host=params['tree_host'],
            base_dir=Path(params['base_dir']),
            time_per_target=params['time_per_target'], 
            max_ppg=params['max_ppg'],
            max_branching=params['max_branching'],
        )

    return planner.get_save_trees(targets)

def build_recommender(params): 
    rec = params['recommender']
    if rec == 'api': 
        return AskcosAPIRecommender(host=params['context_host'])
    elif rec == 'local': 
        return AskcosRecommender(askcos_path=params['askcos_path'])
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Context recommender {rec} not implemented')

def build_scorer(params): 
    rec = params['recommender']
    if rec == 'api': 
        return AskcosAPIScorer(host=params['scorer_host'])
    elif rec == 'local': 
        return AskcosScorer(askcos_path=params['askcos_path'])
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Scorer {rec} not implemented')
    
def build_coster(params): 
    rec = params['coster']
    if rec == 'chemspace': 
        sys.path.append(params['key_path'])
        from keys import chemspace_api_key
        return ChemSpaceCoster(api_key=chemspace_api_key)
    elif rec == 'naive': 
        return NaiveCoster()
    elif rec == 'lookup': 
        raise NotImplementedError
    else:
        raise NotImplementedError(f'Scorer {rec} not implemented')
     
def build_selector(params, target_dict, storage_path):

    graph = RouteGraph(node_filename=storage_path)

    weights = [params['reward_weight'], params['start_cost_weight'], params['reaction_weight']]

    selector = RouteSelector(
        target_dict=target_dict,
        route_graph=graph, 
        condition_recommender=build_recommender(params), 
        output_dir=Path(params['base_dir']),
        rxn_scorer=build_scorer(params),
        coster=build_coster(params),
        weights=weights,
        constrain_all_targets=params['constrain_all']
    )

    filepath = Path(params['base_dir'])/'trees_w_info.json'
    print(f'Saving route graph with contexts, reaction scores, and costs to {filepath}')
    selector.graph.to_json(filepath)
    
    return selector

def run():
    args = get_args()
    params = vars(args)

    print('SPARROW will be run with the following parameters:')    
    for k, v in sorted(params.items()):
        if v is not None: 
            print(f"  {k}: {v}")
    print(flush=True)

    base_dir = Path(params['base_dir'])
    target_dict, targets = get_target_dict(params['target_csv']) 
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path)
    selector = optimize(selector, params)
    storage = extract_vars(selector, base_dir)

if __name__ == '__main__':
    run()
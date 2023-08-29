import requests 

HOST = 'http://18.4.94.12'

params = {
    'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1', # required
    
    # optional with defaults shown
    'max_depth': 4,
    'max_branching': 25,
    'expansion_time': 20,
    'max_ppg': 10,
    'template_count': 100,
    'max_cum_prob': 0.995,
    'chemical_property_logic': 'none',
    'max_chemprop_c': 0,
    'max_chemprop_n': 0,
    'max_chemprop_o': 0,
    'max_chemprop_h': 0,
    'chemical_popularity_logic': 'none',
    'min_chempop_reactants': 5,
    'min_chempop_products': 5,
    'filter_threshold': 0.75,
    
    'return_first': 'false'
}
# resp = requests.get(HOST+'/api/treebuilder/', params=params, verify=False)

# resp = requests.get('http://ASKCOS/askcos/askcos_site/api/tree_builder.py', params=params, verify=False)

params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', #required
    'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1', #required
    
    'num_results': 5, # default is 10
    'return_scores': 'true' # default is false
}

resp = requests.get(HOST+'/api/context/', params=params, verify=False)

print(resp.json())
from utils import create_tree_html, construct_tree_from_graph
from pathlib import Path
from pulp import LpProblem, GUROBI
import pickle


weights = [1,1,1,500]

case = 'amd_darpa'
with open(Path(case)/'rxn_le.pkl', 'rb') as f: 
    rxn_le = pickle.load(f)

with open(Path(case)/'problem.pkl', 'rb') as f: 
    prob_dict = pickle.load(f)

dict_path = Path(case) / 'dictionaries'
stor = {}
for d in list(dict_path.glob('*')): 
	with open(d, 'rb') as f: 
		stor[d.stem] = pickle.load(f)

encoded_rxn_dict = stor['encoded_rxn_dict']
encoded_cond_dict = stor['encoded_cond_dict']
chem_dict_start = stor['chem_dict_start']
chem_dict_inter = stor['chem_dict_inter']
chem_dict_tar = stor['chem_dict_tar']
target_dict = stor['target_dict']

_id = len(encoded_rxn_dict)
dummy_edg_start = len(encoded_rxn_dict)
for key in chem_dict_start.keys():
    encoded_rxn_dict[_id]={key:1,"cond":[['']]*10,'score':[1]*10,'penalty':[0]*10}
    chem_dict_start[key][0].append(_id)
    _id+=1
dummy_edg = range(dummy_edg_start,_id)
real_edg = range(dummy_edg_start)


_, prob = LpProblem.from_dict(prob_dict)
prob.writeLP('StockOptimizationNetworkAnalysis.lp')
prob.solve(GUROBI(timeLimit=86400))


print('\n')
print("***************************************************************")
print("Summarizing optimal solution...")
print("***************************************************************")
print('\n')

selected_rxns_list=[]
selected_rxns_id_list=[]
selected_rxns_times=[]
selected_cond_dict= {}
for v in prob.variables():
#     if 'reactions' in v.name and v.varValue>0.1 and 'bin' not in v.name:
	if 'option' in v.name and v.varValue>0.1:
	#         print(v.name, '=', v.varValue)
		rxn_no = int(v.name.split('_')[0][6:])
		cond_no = int(v.name.split('_')[1])
		if rxn_no in real_edg:
			selected_rxns_id_list.append(rxn_no)
			selected_rxns_list.append(rxn_le[rxn_no])
			selected_rxns_times.append(int(v.varValue))
			selected_cond_dict[rxn_no]=cond_no

used_start = set()
used_inter = set()
used_cond = set()
used_rxns = set()
for rxn_id in selected_rxns_id_list:
	used_rxns.add(rxn_id)
	for chem in encoded_rxn_dict[int(rxn_id)].keys():
		if chem in chem_dict_start:
			used_start.add(chem)
		elif chem in chem_dict_inter or (chem in chem_dict_tar and encoded_rxn_dict[int(rxn_id)][chem]==-1):
			used_inter.add(chem)
	for cond in encoded_rxn_dict[int(rxn_id)]['cond'][selected_cond_dict[int(rxn_id)]]:
		if cond!='':
			used_cond.add(cond)
		
# print(len(selected_cond_dict))           
print("number of starting materials: ", len(used_start))
print("number of intermediates: ",len(used_inter))
print("number of catalysts/solvents/reagents: ",len(used_cond))


####
print('\n')
print("***************************************************************")
print('Writing results files...')
print("***************************************************************")
print('\n')

storage = {'used_inter': used_inter, 'used_start': used_start, 'used_rxns': used_rxns}
with open(Path(case)/'optimal_vars.pkl', 'wb') as f: 
	pickle.dump(storage, f)

suffix = '_'.join(map(str,weights))+'_no_decomp'
trees = []
for target in target_dict:
    tree = construct_tree_from_graph(target, used_inter, used_start, used_rxns)
    trees.append({'tree':tree})
#     print(get_cum_score_from_tree(tree))
create_tree_html(trees, case+'/'+suffix)
print('Done! See results file at {}'.format(case+'/'+suffix+'.html'))




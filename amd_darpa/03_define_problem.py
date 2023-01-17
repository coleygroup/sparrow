import sys 
sys.path.append('/home/jfromer/Molecule_library_synthesis/')
import pickle
import pandas as pd
from rdkit import Chem
import collections
import urllib
import argparse
from tqdm import tqdm
from pathlib import Path
# from pulp import *
from askcos.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI
from gurobipy import *
import csv

# options, change for each case 
case = 'amd_darpa'
# weights = [10,1,0.1,22] [start, CSR, penalty, acuisiton]
weights = [1,1,1,500]
constrain_all_targets = 0

dict_path = Path(case) / 'dictionaries'

stor = {}
for d in list(dict_path.glob('*')): 
	with open(d, 'rb') as f: 
		stor[d.stem] = pickle.load(f)

with open(Path(case)/'rxn_le.pkl', 'rb') as f: 
	rxn_le = pickle.load(f)

encoded_rxn_dict = stor['encoded_rxn_dict']
encoded_cond_dict = stor['encoded_cond_dict']

chem_dict_start = stor['chem_dict_start']
chem_dict_inter = stor['chem_dict_inter']
chem_dict_tar = stor['chem_dict_tar']

target_dict = stor['target_dict']
rewards = stor['rewards']

cond_le = stor['cond_le']
cond_le_rev = stor['cond_le_rev']


###create dummy rxns###
_id = len(encoded_rxn_dict)
dummy_edg_start = len(encoded_rxn_dict)
for key in chem_dict_start.keys():
    encoded_rxn_dict[_id]={key:1,"cond":[['']]*10,'score':[1]*10,'penalty':[0]*10}
    chem_dict_start[key][0].append(_id)
    _id+=1
dummy_edg = range(dummy_edg_start,_id)
real_edg = range(dummy_edg_start)

# print('\n')
# print('number of starting materials: ',len(encoded_start_dict))
# print('number of intermediates: ',len(inter_set))
# print('number of catalysts/solvents/reagents: ',len(encoded_cond_dict))
# print('number of reactions: ', num_rxns)
# print('\n')


####combined planning without decomposition

print('\n')
print("***************************************************************")
print("Performing combined synthesis planning without decomposition...")
print("***************************************************************")
print('\n')

# suffix = '_'.join(map(str,weights))+'_no_decomp'
r = LpVariable.dicts('reactions',encoded_rxn_dict.keys(),lowBound=0,upBound=len(target_dict),cat="Continuous")
rb = LpVariable.dicts('start', encoded_rxn_dict.keys(),cat="Binary")
m = LpVariable.dicts('cond', encoded_cond_dict.keys(), lowBound=0,upBound=1,cat='Continuous')
O = LpVariable.dicts('option', encoded_rxn_dict.keys())
for i in encoded_rxn_dict.keys():
	O[i]=LpVariable.dicts('option'+str(i), range(len(encoded_rxn_dict[i]['cond'])),lowBound=0,upBound=1,cat='Binary')


prob = LpProblem('Network_analysis', LpMinimize)
M=len(target_dict)
#objective
prob += weights[0]*lpSum([rb[i] for i in dummy_edg])\
		+ weights[1]*lpSum([m[k] for k in cond_le])\
		+ weights[2]*lpSum([lpSum([encoded_rxn_dict[i]['penalty'][j]*O[i][j] for j in range(len(encoded_rxn_dict[i]['cond']))]) for i in real_edg])
		

#constraints
# for any target, choose exactly one pathway
# TODO: loop through targets within encoded_rxn_dict instead
for target in tqdm(list(target_dict), desc='Target constraints'):
	# rxn_to_target = []
	# rxn_from_target = []
	# for key, value_dict in encoded_rxn_dict.items(): 
	# 	if target in value_dict and (value_dict[target]==1): 
	# 		rxn_to_target.append(key)
	# 	if target in value_dict and (value_dict[target]==-1): 
	# 		rxn_from_target.append(key)
	rxn_to_target = [key for key,value_dict in encoded_rxn_dict.items() if (target in value_dict) and (value_dict[target]==1)]
	rxn_from_target = [key for key,value_dict in encoded_rxn_dict.items() if (target in value_dict) and (value_dict[target]==-1)]
	prob += lpSum([r[i] for i in rxn_to_target])<=1+lpSum([r[i] for i in rxn_from_target])
	prob += prob.objective - weights[3]*rewards[target]*(lpSum([r[i] for i in rxn_to_target]) - lpSum([r[i] for i in rxn_from_target]))
		

# for all intermediates either some reaction in and some reaction out, or no reaction in or no reaction out

for inter,related_rxns in tqdm(chem_dict_inter.items(), desc='Intermediate constraints'):
	rxn_to_inter = related_rxns[0]
	rxn_from_inter = related_rxns[1]
	prob+= lpSum([r[i] for i in rxn_to_inter])==lpSum([r[j] for j in rxn_from_inter])


for start,related_rxns in tqdm(chem_dict_start.items(), desc='Starting material constraints'):
	rxn_to_start = related_rxns[0]
	rxn_from_start = related_rxns[1]
	prob+= lpSum([r[i] for i in rxn_to_start])==lpSum([r[j] for j in rxn_from_start])

#if one reaction uses the starting material is selected, the starting material is selected
for key in tqdm(encoded_rxn_dict.keys(), desc='Flow constraints'):
	prob += rb[key]*M>=r[key]
	prob += rb[key]<=r[key]

for key in real_edg:
	prob+= lpSum([O[key][j] for j in range(len(encoded_rxn_dict[key]['cond']))])==rb[key]

for i in real_edg:
	for j in range(len(encoded_rxn_dict[i]['cond'])):
		for k in encoded_rxn_dict[i]['cond'][j]:
			if k !="":
				try:
					prob+= m[cond_le_rev[k]]>=O[i][j]
				except: 
					prob+= m[cond_le_rev[k.encode('utf-8')]]>=O[i][j]
					
print('problem is set up, start solving...')

prob_dict = prob.to_dict()
with open(Path(case)/'problem.pkl', 'wb') as f: 
	pickle.dump(prob_dict, f)





import sys 
sys.path.append('/home/jfromer/Molecule_library_synthesis/')
import pickle
import pandas as pd
from rdkit import Chem
import collections
import urllib
import argparse
from tqdm import tqdm
# from pulp import *
from askcos.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from pulp import LpVariable, LpProblem, LpMinimize, lpSum, GUROBI
from gurobipy import *
from rdkit.Chem import QED

# options, change for each case 
case = 'case_4'
# weights = [10,1,0.1,22] [start, condition, penalty, acuisiton]
weights = [1,0.1,0.1,5]
constrain_all_targets = 0

## load reactions 
with open(case+'/reaction_dict.pickle','rb') as RD:
    rxn_dict = pickle.load(RD)

##load the list of targets
target_df = pd.read_csv(case+'/target_list.csv')
target_list = list(target_df.loc[target_df['numberofpaths']>0]['SMILES'])
target_list = [Chem.MolToSmiles(Chem.MolFromSmiles(target),False) for target in target_list]
print(str(len(target_list)) + " targets")
target_dict = set(target_list)

# calculate rewards
rewards = {smi: (QED.default(Chem.MolFromSmiles(smi))) for smi in target_list}


###label rxns
i=0
rxn_le = {}
encoded_rxn_dict = {}
for key,value in rxn_dict.items():
    rxn_le[i]=key
    encoded_rxn_dict[i]=value
    i+=1

###load evaluated reaction dictionary
with open(case+'/encoded_rxn_dict_with_cond.pkl','rb') as RXN_DICT:
	try:
		encoded_rxn_dict = pickle.load(RXN_DICT,encoding="byte")
	except: 
		encoded_rxn_dict = pickle.load(RXN_DICT,encoding="latin1")
	
with open(case+'/cond_dict.pkl','rb') as COND_DICT:
    cond_dict = pickle.load(COND_DICT,encoding="bytes")
num_rxns=len(encoded_rxn_dict)
####calculate penalty

for key,value in encoded_rxn_dict.items():
	try:
		encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
			for key,score in encoded_rxn_dict[key]['score'].items()}
	except:
		encoded_rxn_dict[key]['penalty']= {key:1/score if score>0.05 else 20 \
			for key,score in encoded_rxn_dict[key]['score'.encode('utf-8')].items()}


###construct chemical dictionary
chem_dict_start = {}
chem_dict_tar = {}
chem_dict_inter = {}
chem_le={}
chem_le_rev={}
_id = 0
for key,value in encoded_rxn_dict.items():
    for c, stoi in value.items():
#         print c
        if c=='cond' or c=='score' or c=='penalty':
            continue
        if c not in chem_le:
            chem_le[c]=_id
            chem_le_rev[_id]=c
            _id+=1
        if c in chem_dict_start:
            chem_dict_start[c][1].append(key)
            continue
        elif c in chem_dict_tar:
            chem_dict_tar[c].append(key)
            continue
        elif c in chem_dict_inter:
            if stoi==1:
                chem_dict_inter[c][0].append(key)
            else:
                chem_dict_inter[c][1].append(key)
            continue
        elif c in target_dict:
            chem_dict_tar[c]=[key]
        else:
            m = Chem.MolFromSmiles(c)
            nc = 0
            nn = 0
            no = 0
            for a in m.GetAtoms():
                if a.GetSymbol() == 'C':
                    nc +=1
                if a.GetSymbol() == 'N':
                    nn +=1
                if a.GetSymbol() == 'O':
                    no +=1
            if nc<=10 and nn<=3 and no<=5:
                chem_dict_start[c]=[[],[key]]
            else:
                if stoi ==1:
                    chem_dict_inter[c] = [[key],[]]
                else:
                    chem_dict_inter[c] = [[],[key]]

###create dummy rxns###
_id = len(encoded_rxn_dict)
dummy_edg_start = len(encoded_rxn_dict)
for key in chem_dict_start.keys():
    encoded_rxn_dict[_id]={key:1,"cond":[['']]*10,'score':[1]*10,'penalty':[0]*10}
    chem_dict_start[key][0].append(_id)
    _id+=1
dummy_edg = range(dummy_edg_start,_id)
real_edg = range(dummy_edg_start)

###create a copy of reaction dictionary with only chemical information
encoded_rxn_dict_chem_only = {}
for key,value in encoded_rxn_dict.items():
    encoded_rxn_dict_chem_only[key]={}
    for c,stoi in value.items():
        if c=='cond' or c=='score' or c=='penalty':
            continue
        else:
            encoded_rxn_dict_chem_only[key][c]=stoi

##calculate SCSCores###
scscorer = SCScorePrecursorPrioritizer()
scscorer.load_model(model_tag='1024bool')

encoded_start_dict = {}
start_set= set([])
scscore_dict = {}
for key, value in chem_dict_start.items():
    start_set.add(chem_le[key])
    encoded_start_dict[chem_le[key]]=value
    scscore_dict[chem_le[key]] = scscorer.get_score_from_smiles(key, noprice=True)

###encode the intermediates
inter_set = set([])
encoded_inter_dict = {}
for key, value in chem_dict_inter.items():
    encoded_inter_dict[chem_le[key]]=value
    inter_set.add(chem_le[key])

####encode condition dictionary######
i=0
cond_le = {}
cond_le_rev = {}
encoded_cond_dict = {}
for key, value in cond_dict.items():
    cond_le[i]=key
    cond_le_rev[key]=i
    encoded_cond_dict[i]=value
    i+=1

###encode target dictionary MAYBE NOT NECESSARY
tar_set = set([])
encoded_tar_dict = {}
for key, value in chem_dict_tar.items():
    tar_set.add(chem_le[key])
    encoded_tar_dict[chem_le[key]]=value

print('\n')
print('number of starting materials: ',len(encoded_start_dict))
print('number of intermediates: ',len(inter_set))
print('number of catalysts/solvents/reagents: ',len(encoded_cond_dict))
print('number of reactions: ', num_rxns)
print('\n')

####combined planning without decomposition

print('\n')
print("***************************************************************")
print("Performing combined synthesis planning without decomposition...")
print("***************************************************************")
print('\n')

suffix = '_'.join(map(str,weights))+'_no_decomp'
r = LpVariable.dicts('reactions',encoded_rxn_dict.keys(),lowBound=0,upBound=len(target_dict),cat="Continuous")
rb = LpVariable.dicts('start', encoded_rxn_dict.keys(),cat="Binary")
m = LpVariable.dicts('cond', encoded_cond_dict.keys(), lowBound=0,upBound=1,cat='Continuous')
O = LpVariable.dicts('option', encoded_rxn_dict.keys())
for i in encoded_rxn_dict.keys():
	O[i]=LpVariable.dicts('option'+str(i), range(len(encoded_rxn_dict[i]['cond'])),lowBound=0,upBound=1,cat='Binary')


prob = LpProblem('Network analysis', LpMinimize)
M=len(target_dict)
#objective
prob += weights[0]*lpSum([rb[i] for i in dummy_edg])\
		+ weights[1]*lpSum([m[k] for k in cond_le])\
		+ weights[2]*lpSum([lpSum([encoded_rxn_dict[i]['penalty'][j]*O[i][j] for j in range(len(encoded_rxn_dict[i]['cond']))]) for i in real_edg])
		

#constraints
# for any target, choose exactly one pathway
for target in tqdm(list(target_dict)):
	rxn_to_target = [key for key,value_dict in encoded_rxn_dict.items() if (target in value_dict) and (value_dict[target]==1)]
	rxn_from_target = [key for key,value_dict in encoded_rxn_dict.items() if (target in value_dict) and (value_dict[target]==-1)]
	if constrain_all_targets:  
		prob += lpSum([r[i] for i in rxn_to_target])==1+lpSum([r[i] for i in rxn_from_target])
	else: 
		prob += lpSum([r[i] for i in rxn_to_target])<=1+lpSum([r[i] for i in rxn_from_target])
		prob += prob.objective - weights[3]*rewards[target]*(lpSum([r[i] for i in rxn_to_target]) - lpSum([r[i] for i in rxn_from_target]))

# for all intermediates either some reaction in and some reaction out, or no reaction in or no reaction out

for inter,related_rxns in tqdm(chem_dict_inter.items()):
	rxn_to_inter = related_rxns[0]
	rxn_from_inter = related_rxns[1]
	prob+= lpSum([r[i] for i in rxn_to_inter])==lpSum([r[j] for j in rxn_from_inter])


for start,related_rxns in tqdm(chem_dict_start.items()):
	rxn_to_start = related_rxns[0]
	rxn_from_start = related_rxns[1]
	prob+= lpSum([r[i] for i in rxn_to_start])==lpSum([r[j] for j in rxn_from_start])

#if one reaction uses the starting material is selected, the starting material is selected
for key in tqdm(encoded_rxn_dict.keys()):
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


prob.writeLP('stock_optimization_network_analysis.lp')
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

def construct_tree_from_graph(target, used_inter, used_start, used_rxns):
    def find_buyable_path(node, seen_chemicals):
        if node in used_start:
            return {'is_chemical':True, 'smiles':node, 'children':[]}
        if node==target:
            seen_chemicals.add(node)
            node = {'is_chemical':True, 'smiles':target, 'children':[] }
            rxns = []
            rxns_to_target = [key for key,value_dict in encoded_rxn_dict.items() if (target in value_dict) and (value_dict[target]==1)]
            for rxn in rxns_to_target:
                if rxn in used_rxns:
                    rxns.append(rxn)
            for rxn in rxns:
                
                node['children'].append(find_buyable_path(rxn, seen_chemicals))
            if not any(node['children']):
                return
            return node
        if node in used_rxns:
#             if node in seen_rxns:
#                 return
#             seen_rxns.add(node)
            rxn = node
            node = {'is_reaction':True, 'smiles':rxn_le[rxn],'forward_score':encoded_rxn_dict[rxn]['score'][selected_cond_dict[rxn]],'context':encoded_rxn_dict[rxn]['cond'][selected_cond_dict[rxn]], 'children':[]}
            chemicals = [key for key,value in encoded_rxn_dict[rxn].items() if value ==-1]
            for chemical in chemicals:
                node['children'].append(find_buyable_path(chemical, seen_chemicals))
            if not any(node['children']):
                return
            return node
        if node in used_inter:
            if node in seen_chemicals:
                return
            seen_chemicals.add(node)
            inter = node
            node = {'is_chemical':True, 'smiles':inter, 'children':[]}
            rxns = []
            try:
                rxns_to_inter = chem_dict_inter[inter][0]
            except:
                rxns_to_inter = [key for key,value_dict in encoded_rxn_dict.items() if (inter in value_dict) and (value_dict[inter]==1)]
            for rxn in rxns_to_inter:
                if isinstance(used_rxns, dict):
                    if rxn in used_rxns and used_rxns[rxn]>0:
                        rxns.append(rxn)
                        used_rxns[rxn]-=1
                        break
                else:
                    if rxn in used_rxns:
                        rxns.append(rxn)
            for rxn in rxns:
                node['children'].append(find_buyable_path(rxn, seen_chemicals))
            if not any(node['children']):
                return
            return node
    return find_buyable_path(target, set())

def construct_tree_for_d3_visualization(tree,depth,new_tree = {}):
    if 'is_chemical' in tree.keys():
#         new_tree['smiles']='http://askcos.mit.edu/draw/smiles/'+str(urllib.quote(tree['smiles'],safe= ''))
        new_tree['smiles']='http://askcos.mit.edu/draw/smiles/'+str(tree['smiles']).replace('#','%23')
        new_tree['rc_type'] ='chemical'
        try:
            new_tree['freq'] = chemical_freq_dict[tree['smiles']]
        except:
            pass
        new_tree['smiles']=str(new_tree['smiles'])
    else:
        new_tree['smiles']='http://askcos.mit.edu/draw/smiles/'
        new_tree['names']=''
        new_tree['score']='%.3f' % tree['forward_score']
        new_tree['_id']=int(depth)
        for i in range(5):
            for c in tree['context'][i].split('.'):   
                if 'Reaxys' not in c:
                    new_tree['smiles']+='.'+str(urllib.parse.quote(c,safe= '')).replace('#','%23')
                elif 'Reaxys Name' in c:
                    new_tree['names']+=str(c)[11:]+'.'
                else:
                    new_tree['names']+=str(c)
        if new_tree['smiles'] == "http://askcos.mit.edu/draw/smiles/.....":
            new_tree['smiles'] = ''
        new_tree['names']=str(new_tree['names'])
        new_tree['rc_type']='reaction'
    new_tree['children'] = []
    if tree['children']:
        for child in tree['children']:
            new_tree['children'].append({})
            construct_tree_for_d3_visualization(child,depth+0.5, new_tree['children'][-1])
    return new_tree

def get_cum_score_from_tree(tree):
    cum_score=1

    if tree['children']:
        if 'is_reaction' in tree:
            cum_score = tree['forward_score']
        for child in tree['children']:
            if 'is_chemical' in tree:
                cum_score = get_cum_score_from_tree(child)
            else:
                cum_score *= get_cum_score_from_tree(child)
    return cum_score

def create_tree_html(trees,file_name):
	height = 200*len(trees)
	try:
		outfile = file_name

	except Exception as e:
		print(e)
		print('Need to specify file name to write results to')

	trees_for_visualization = {'name': 'dummy_root','children':[]}
	for i,tree in enumerate(trees):
		if tree['tree']: 
			trees_for_visualization['children'].append(construct_tree_for_d3_visualization(tree['tree'],depth = 0.5, new_tree={}))
			trees_for_visualization['children'][-1]['_id'] = ('T%d' % i)
			trees_for_visualization['children'][-1]['score'] = '%.3f'%get_cum_score_from_tree(tree['tree'])
	fid_out = open(outfile+'.html','w',encoding='utf-8')
	fid_out.write('<!DOCTYPE html>\n')
	fid_out.write('  <head>\n')
	fid_out.write('    <meta charset="utf-8">\n')
	fid_out.write('    <title>{}</title>\n'.format(outfile))
	fid_out.write('    <style>\n')
	fid_out.write('	.node circle {\n')
	fid_out.write('	  fill: #fff;\n')
	fid_out.write('	  stroke: steelblue;\n')
	fid_out.write('	  stroke-width: 3px;\n')
	fid_out.write('	}\n')
	fid_out.write('	.node rect {\n')
	fid_out.write('		fill: #fff;\n')
	fid_out.write('		stroke: steelblue;\n')
	fid_out.write('		stroke_width: 3px;\n')
	fid_out.write('	}\n')
	fid_out.write('	.node text { font: 12px sans-serif; }\n')
	fid_out.write('	.link {\n')
	fid_out.write('	  fill: none;\n')
	fid_out.write('	  stroke: #ccc;\n')
	fid_out.write('	  stroke-width: 2px;\n')
	fid_out.write('	}\n')
	fid_out.write('    </style>\n')
	fid_out.write('  </head>\n')
	fid_out.write('  <body>\n')
	fid_out.write('<!-- load the d3.js library -->	\n')
	fid_out.write('<script src="http://d3js.org/d3.v3.min.js"></script>\n')
	fid_out.write('<script>\n')
	fid_out.write('var treeData = [\n')
	fid_out.write('{}\n'.format(trees_for_visualization))
	fid_out.write('];\n')
	fid_out.write('var margin = {top: 20, right: 120, bottom: 20, left: 0},\n')
	fid_out.write('	width = 3000 - margin.right - margin.left,\n')
	fid_out.write('	height = {} - margin.top - margin.bottom;\n'.format(height))
	fid_out.write('var i = 0;\n')
	fid_out.write('var tree = d3.layout.tree()\n')
	fid_out.write('	.size([height, width]);\n')
	fid_out.write('var diagonal = d3.svg.diagonal()\n')
	fid_out.write('	.projection(function(d) { return [d.y, d.x]; });\n')
	fid_out.write('var svg = d3.select("body").append("svg")\n')
	fid_out.write('	.attr("width", width + margin.right + margin.left)\n')
	fid_out.write('	.attr("height", height + margin.top + margin.bottom)\n')
	fid_out.write('  .append("g")\n')
	fid_out.write('	.attr("transform", \n')
	fid_out.write('	      "translate(" + margin.left + "," + margin.top + ")");\n')
	fid_out.write('root = treeData[0];\n')
	fid_out.write('update(root);\n')
	fid_out.write('function update(source) {\n')
	fid_out.write('  // Compute the new tree layout.\n')
	fid_out.write('  var nodes = tree.nodes(root).reverse(),\n')
	fid_out.write('	  links = tree.links(nodes);\n')
	fid_out.write('  // Normalize for fixed-depth.\n')
	fid_out.write('  nodes.forEach(function(d) { d.y = d.depth * 120; });\n')
	fid_out.write('  // Declare the nodes…\n')
	fid_out.write('  var node = svg.selectAll("g.node")\n')
	fid_out.write('	  .data(nodes, function(d) { return d.id || (d.id = ++i); });\n')
	fid_out.write('  // Enter the nodes.\n')
	fid_out.write('  var nodeEnter = node.enter().append("g")\n')
	fid_out.write('	  .attr("class", "node")\n')
	fid_out.write('	  .attr("transform", function(d) { \n')
	fid_out.write('		  return "translate(" + d.y + "," + d.x + ")"; });\n')
	fid_out.write('  nodeEnter.append("image")\n')
	fid_out.write('      .attr("xlink:href", function(d) { return d.smiles; })\n')
	fid_out.write('      .attr("x", "-80px")\n')
	fid_out.write('      .attr("y", "-35px")\n')
	fid_out.write('      .attr("width", "80px")\n')
	fid_out.write('      .attr("height", "80px");\n')
	fid_out.write('  nodeEnter.append("path")\n')
	fid_out.write('  	  .style("stroke", "black")\n')
	fid_out.write('  	  .style("fill", function(d) { if (d.freq==1) { return "white"; }\n')
	fid_out.write('  	  								else if (d.freq==2) { return "yellow";}\n')
	fid_out.write('  	  								else if (d.freq==3) { return "orange"; }\n')
	fid_out.write('  	  								else if (d.freq>=4) { return "red"; }\n')
	fid_out.write('  	  								else {return "white";}\n')
	fid_out.write('  	  								})\n')
	fid_out.write('  	  .attr("d", d3.svg.symbol()\n')
	fid_out.write('  	  				.size(20)\n')
	fid_out.write('  	  				.type(function(d) {if\n')
	fid_out.write('  	  					(d.rc_type == "chemical") {return "circle";} else if\n')
	fid_out.write('  	  					(d.rc_type == "reaction") {return "cross";}\n')
	fid_out.write('  	  				}));\n')
	fid_out.write('  nodeEnter.append("text")\n')
	fid_out.write('	  .attr("x", 0)\n')
	fid_out.write('	  .attr("y", 35)\n')
	fid_out.write('	  .attr("text-anchor", function(d) { \n')
	fid_out.write('		  return d.children || d._children ? "end" : "start"; })\n')
	fid_out.write('	  .text(function(d) { return d.names; })\n')
	fid_out.write('	  .style("fill-opacity", 1);\n')
	fid_out.write('  nodeEnter.append("text")\n')
	fid_out.write('	  .attr("x", -60)\n')
	fid_out.write('	  .attr("y", -30)\n')
	fid_out.write('	  .attr("text-anchor", function(d) { \n')
	fid_out.write('		  return d.children || d._children ? "end" : "start"; })\n')
	fid_out.write('	  .text(function(d) { return d._id; })\n')
	fid_out.write('	  .style("fill-opacity", 1);\n')
	fid_out.write('  nodeEnter.append("text")\n')
	fid_out.write('	  .attr("x", 0)\n')
	fid_out.write('	  .attr("y", -30)\n')
	fid_out.write('	  .attr("text-anchor", function(d) { \n')
	fid_out.write('		  return d.children || d._children ? "end" : "start"; })\n')
	fid_out.write('	  .text(function(d) { return d.score; })\n')
	fid_out.write('	  .style("fill-opacity", 1);\n')
	fid_out.write('  // Declare the links…\n')
	fid_out.write('  var link = svg.selectAll("path.link")\n')
	fid_out.write('	  .data(links, function(d) { return d.target.id; });\n')
	fid_out.write('  // Enter the links.\n')
	fid_out.write('  link.enter().insert("path", "g")\n')
	fid_out.write('	  .attr("class", "link")\n')
	fid_out.write('	  .style("stroke", function(d) { return d.target.level; })\n')
	fid_out.write('	  .attr("d", diagonal);\n')
	fid_out.write('  // remove the first level, leaving the targets as the first level\n')
	fid_out.write('  node.each(function(d){\n')
	fid_out.write('	if (d.name == "dummy_root")\n')
	fid_out.write('    	d3.select(this).remove();});\n')
	fid_out.write('	link.each(function(d){\n')
	fid_out.write('	if (d.source.name == "dummy_root") \n')
	fid_out.write('	    d3.select(this).remove();});\n')
	fid_out.write('}\n')
	fid_out.write('</script>\n')
	fid_out.write('  </body>\n')
	fid_out.write('</html>\n')

	fid_out.close()

trees = []
for target in target_dict:
    tree = construct_tree_from_graph(target, used_inter, used_start, used_rxns)
    trees.append({'tree':tree})
#     print(get_cum_score_from_tree(tree))
create_tree_html(trees, case+'/'+suffix)

save_info = [{'trees': trees, 'used_inter':used_inter, 'used_start':used_start, 'used_rxns': used_rxns}]
with open(case+'/'+suffix + '.pkl', 'wb') as handle:
    pickle.dump(save_info, handle, protocol=pickle.HIGHEST_PROTOCOL)


print('Done! See results file at {}'.format(case+'/'+suffix+'.html'))
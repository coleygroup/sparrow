def construct_tree_from_graph(
    target, 
    used_inter, 
    used_start, 
    used_rxns, 
    encoded_rxn_dict, 
    selected_cond_dict, 
    rxn_le, 
    chem_dict_inter
):
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
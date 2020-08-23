# from makeit.retrosynthetic.tree_builder import TreeBuilder
from makeit.retrosynthetic.mcts.tree_builder import MCTS
from makeit.synthetic.evaluation.tree_evaluator import TreeEvaluator
from makeit.utilities.io.logger import MyLogger
import makeit.global_config as gc
import pandas as pd
from makeit.utilities.io import name_parser
from rdkit import Chem
import pickle
import pprint
import sys
# from construct_path_chemical_matrix import get_chemical_table_from_routes

MyLogger.initialize_logFile()
results_folder = sys.argv[1]+'/'

print('start to load targets')
# dataset = pd.read_csv('common_substructure.csv')
# dataset = pd.read_csv('Top30DrugsByPrescription.csv')
# dataset.head()

# with open('truncated_mol_lib.pickle','r') as MOLLIB:
# with open('WHO_EM/WHO_truncated.pickle','r') as MOLLIB:
#     truncated_mol_lib = pickle.load(MOLLIB)
smiles_df=pd.read_csv('targets.csv')
smiles_list = list(smiles_df['smiles'])

n_paths = []

truncated_mol_lib = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
print(len(truncated_mol_lib))# for name in dataset['SMILES']:

print("{} target molecules to expand".format(len(truncated_mol_lib)))
celery = False
NCPUS = 4
# treeBuilder = TreeBuilder(celery=celery, mincount=25, mincount_chiral=10)
Tree = MCTS(nproc=NCPUS, mincount=gc.RETRO_TRANSFORMS_CHIRAL['mincount'], 
        mincount_chiral=gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'],
        celery=celery)

all_routes = []
all_trees = []
all_routes_contexts = []
index_list_all_routes = []
is_first_target = True
status = 0

status_file = open(results_folder+'status for molecules.dat','w')
status_file.write('SMILES\tnumberofpaths\tnotes\n')
reaction_dict = {}
for mol in truncated_mol_lib:
# for mol in [Chem.MolFromSmiles('OC(C1CCCCN1)c2cc(nc3c2cccc3C(F)(F)F)C(F)(F)F')]:
    if is_first_target:
        soft_reset = False
        is_first_target = False
    else:
        soft_reset = True
# for mol in [Chem.MolFromSmiles('O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1')]:
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)   
        print('expanding target {}'.format(smiles)) 
        status, paths = Tree.get_buyable_paths(smiles,
                                            nproc=NCPUS,
                                            max_depth=10,
                                            expansion_time=60,
                                            max_trees=10000, 
                                            max_cum_template_prob=0.995,
                                            max_branching=25, apply_fast_filter=True, filter_threshold=0.75,
                                            template_count=1000,
                                            max_natom_dict={'C':10,'N':3,'O':5,'logic': 'or'},
                                            max_ppg=0,
                                            # min_chemical_history_dict={'as_reactant':1, 'as_product':1,'logic':'and'},
                                            soft_reset=soft_reset,
                                            soft_stop=True)
        print('done for target {}'.format(smiles))

        for chemical in Tree.Chemicals:
            chemical_obj = Tree.Chemicals[chemical]
            for _id, result in chemical_obj.template_idx_results.items():
                if result.waiting:
                    continue
                for rsmi,R in result.reactions.items():
                    if (not R.valid) or R.price == -1:
                        continue
                    rxn = '.'.join(sorted(R.reactant_smiles)) + '>>' + chemical
                    if rxn not in reaction_dict:
                        reaction_dict[rxn]={}
                        for rct in rsmi.split('.'):
                            reaction_dict[rxn][rct]=-1
                        reaction_dict[rxn][chemical]=1
        print('number of reactions explored: {}'.format(len(reaction_dict)))
        feasible_trees = []
        n_paths.append(len(paths))
    except Exception as e:
        status_file.write('{}\t{}\t{}\t\n'.format(smiles,0,e))
        n_paths.append(0)
        pass    # index_list_all_routes.extend(index_list_one_target)
    print(len(all_routes))

with open(results_folder+'reaction_dict.pickle','wb') as RD:
    pickle.dump(reaction_dict, RD)

target_list = pd.DataFrame({'SMILES':smiles_list,
                             'numberofpaths':n_paths})
target_list.to_csv(results_folder+'target_list.csv')

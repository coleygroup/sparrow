""" Runs nonlinear formulation of SPARROW over a grid of SM budgets and maximum reactions. Currently just for button_alectinib case """
import subprocess
from pathlib import Path 
import numpy as np 

def opt_rxn_constr(min_rxns=10, max_rxns=60, n_rxn=6):
    """ Optimizes and saves results for Case 1 in preprint (Garibsingh et al DOI:10.1073/pnas.2104093118) """
    conf_path = Path('examples/garibsingh/config_opt.ini')

    cmds = []

    n_rxns = np.r_[0:70:5] # [*np.r_[0:100:10], *np.r_[100:200:20]]

    out_dir = Path(f'results/garib_rxnconstr_nonlinear')
    for max_rxn in n_rxns: 
        out_folder = out_dir / f'rxn_{max_rxn}'
        cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --formulation expected_reward --prune-distance 16 --max-rxns {max_rxn}'
        cmds.append(cmd)

    return cmds


def opt_budget_constr():
    """ Optimizes and saves results for Case 1 in preprint (Garibsingh et al DOI:10.1073/pnas.2104093118) """
    conf_path = Path('examples/garibsingh/config_opt.ini')

    cmds = []

    budgets = [*np.r_[10:160:10], *[175, 200, 225, 250]]
    
    out_dir = Path(f'results/garib_budgetconstr_nonlinear')
    for budget in budgets: 
        out_folder = out_dir / f'budget_{budget}'
        cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --formulation expected_reward --prune-distance 16 --starting-material-budget {budget}'
        cmds.append(cmd)

    return cmds

def opt_budget_constr_linear(): 
    cmds = []
    budgets = [*np.r_[10:160:10], *[175, 200, 225, 250]]
    out_dir = Path(f'results/garib_budgetconstr_linear')
    conf_path = Path('examples/garibsingh/config_opt.ini')
    
    for budget in budgets: 
        out_folder = out_dir / f'budget_{budget}'
        cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --formulation linear --start-cost-weight 0 --starting-material-budget {budget} --reward-weight 40'
        cmds.append(cmd)
    
    return cmds 

def opt_rxn_constr_linear(): 
    cmds = []
    n_rxns = np.r_[0:70:5] # [*np.r_[0:70:10], *np.r_[80:200:20]]
    out_dir = Path(f'results/garib_rxnconstr_linear')
    conf_path = Path('examples/garibsingh/config_opt.ini')
    
    for max_rxn in n_rxns: 
        out_folder = out_dir / f'rxn_{max_rxn}'
        cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --formulation linear --start-cost-weight 0 --max-rxns {max_rxn} --reward-weight 50'
        cmds.append(cmd)
    
    return cmds 

def run_commands(cmds):
    print(f'Running SPARROW {len(cmds)} times')
    for cmd in cmds: 
        print(cmd)
        subprocess.call(cmd, shell=True)

def print_commands(cmds): 
    cmdss = [f'{cmd}\n' for cmd in cmds]
    with open('scripts/cmd_list.txt', 'w') as f: 
        f.writelines(cmdss)

if __name__=='__main__':
    cmds_linear = opt_rxn_constr_linear()
    cmds_nonlinear = opt_rxn_constr()
    # cmds_budget = opt_budget_constr_linear()
    # run_commands([*cmds_rxn, *cmds_budget])
    # print_commands([*cmds_rxn, *cmds_budget])
    print_commands(cmds_nonlinear)
    run_commands(cmds_linear)

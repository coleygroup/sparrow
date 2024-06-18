""" Runs nonlinear formulation of SPARROW over a grid of SM budgets and maximum reactions. Currently just for button_alectinib case """
import subprocess
from pathlib import Path 
import numpy as np 

def opt_grid(min_rxns=10, max_rxns=210, n_rxn=11, min_budget=10, max_budget=210, n_budget=11):
    """ Optimizes and saves results for Case 1 in preprint (Garibsingh et al DOI:10.1073/pnas.2104093118) """
    conf_path = Path('examples/amd/config_opt.ini')

    cmds = []

    n_rxns = np.linspace(min_rxns, max_rxns, n_rxn, dtype = int)
    budgets = np.linspace(min_budget, max_budget, n_budget)
    
    out_dir = Path(f'results/amd_grid_nonlinear')
    for max_rxn in n_rxns: 
        for budget in budgets: 
            out_folder = out_dir / f'rxn_{max_rxn}_budget_{budget}'
            summary_file = out_folder / 'solution' /'summary.json'
            if not summary_file.exists(): 
                cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --formulation expected_reward --prune-distance 16 --max-rxns {max_rxn} --starting-material-budget {budget}'
                
                cmds.append(cmd)

    return cmds



def run_commands(cmds):
    print(f'Running SPARROW {len(cmds)} times')
    for cmd in cmds: 
        print(cmd)
        subprocess.call(cmd, shell=True)

if __name__=='__main__':
    cmds = opt_grid()
    run_commands(cmds)

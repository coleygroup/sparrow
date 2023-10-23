""" Applies SPARROW to the three case studies shown in preprint, creates results used for all plots and routes """
import subprocess
from pathlib import Path 
import numpy as np 

def opt_garib(vary_l1=True, vary_l2=True, vary_l3=True):
    """ Optimizes and saves results for Case 1 in preprint (Garibsingh et al DOI:10.1073/pnas.2104093118) """
    conf_path = Path('examples/garibsingh/config_opt.ini')

    cmds = []

    # vary lambda_1
    out_dir = Path(f'results/garibsingh_rewvary')
    reward_lam = [2, 5, 8, 12.5, 15, 20, 30, 50, 60, 70, 80]
    if vary_l1: 
        for i, lam in enumerate(reward_lam):
            tags = [f'lam-{i}']
            out_folder = out_dir / '_'.join(tags)
            cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --reward-weight {lam} --start-cost-weight 1 --reaction-weight 1'
            
            cmds.append(cmd)


    # vary lambda_2 
    out_dir = Path(f'results/garibsingh_costvary')
    rxn_cost = [0.2, 0.35, 0.7, 1.3, 2, 2.5, 3, 3.5, 6, 12, 15]
    if vary_l2: 
        for i, lam in enumerate(rxn_cost):
            tags = [f'lam-{i}']
            out_folder = out_dir / '_'.join(tags)
            cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --reward-weight 20 --start-cost-weight {lam} --reaction-weight 1'
            
            cmds.append(cmd)


    # vary lambda_3 
    out_dir = Path(f'results/garibsingh_rxnvary')
    rxn_lam = [0.1, 0.25, 0.35, 0.5, 1, 1.5, 1.7, 2, 3, 4, 5]
    if vary_l3: 
        for i, lam in enumerate(rxn_lam):
            tags = [f'lam-{i}']
            out_folder = out_dir / '_'.join(tags)
            cmd = f'sparrow --config {conf_path} --output-dir {out_folder} --reward-weight 20 --start-cost-weight 1 --reaction-weight {lam}'
            
            cmds.append(cmd)

    return cmds


def opt_amd():
    """ Optimizes Case 2 in preprint (Koscher et al DOI:10.26434/chemrxiv-2023-r7b01)
    Retrosynthesis trees provided by Koscher et al """
    conf_path = Path(f'examples/amd/config_opt.ini')
    lam = [None, 1, 1]
    lam_1s = np.linspace(1, 20, 20).round(2)
    out_dir = Path(f'results/amd/lam_{lam[0]}_{lam[1]}_{lam[2]}')
    cmds = []

    for l1 in lam_1s:
        out_dir = Path(f'results/amd/lam_{l1}_{lam[1]}_{lam[2]}')
        cmd = f'sparrow --config {conf_path} --output-dir {out_dir} --reward-weight {l1} --start-cost-weight {lam[1]} --reaction-weight {lam[2]}'
        cmds.append(cmd)

    return cmds

def opt_button():
    """ Optimizes and saves results for Case 3 in preprint (Button et al DOI:10.1038/s42256-019-0067-7) """
    conf_path = Path(f'examples/button_alectinib/config_opt.ini')
    lam = [None, 1, 5]
    lam_1s = [30] # np.linspace(1, 20, 20).round(2)
    cmds = []

    for l1 in lam_1s:
        out_dir = Path(f'results/button_alectinib/lam_{l1}_{lam[1]}_{lam[2]}')
        cmd = f'sparrow --config {conf_path} --output-dir {out_dir} --reward-weight {l1} --start-cost-weight {lam[1]} --reaction-weight {lam[2]}'
        cmds.append(cmd)

    return  

def run_commands(cmds):
    for cmd in cmds: 
        print(cmd)
        subprocess.call(cmd, shell=True)

if __name__=='__main__':
    cmds = opt_garib()
    run_commands(cmds)

    cmds = opt_amd()
    run_commands(cmds)

    cmds = opt_button
    run_commands(cmds)

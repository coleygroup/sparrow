import subprocess
from pathlib import Path 
import numpy as np 
import time

vary_l1 = False
vary_l2 = False 
vary_l3 = True
conf_path = Path('examples/garibsingh/config_opt.ini')


# reward_lam = np.linspace(4, 50, 24)
reward_lam = [2, 5, 8, 12.5, 15, 20, 30, 50, 60, 70, 80]
cmds = []

# vary lambda_1
out_dir = Path(f'results/garibsingh_rewvary')
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

print(f'Running {len(cmds)} molpal runs:')
for cmd in cmds: 
    print(cmd)
    time.sleep(0.5)
    subprocess.call(cmd, shell=True)
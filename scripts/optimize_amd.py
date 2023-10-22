import subprocess
from pathlib import Path 
import numpy as np 

conf_path = Path(f'examples/amd/config_opt.ini')
lam = [None, 1, 1]
lam_1s = np.linspace(1, 20, 20).round(2)
out_dir = Path(f'results/amd/lam_{lam[0]}_{lam[1]}_{lam[2]}')
cmds = []

for l1 in lam_1s:
    out_dir = Path(f'results/amd/lam_{l1}_{lam[1]}_{lam[2]}')
    cmd = f'sparrow --config {conf_path} --output-dir {out_dir} --reward-weight {l1} --start-cost-weight {lam[1]} --reaction-weight {lam[2]}'
    cmds.append(cmd)
        
for cmd in cmds: 
    print(cmd)
    subprocess.call(cmd, shell=True)
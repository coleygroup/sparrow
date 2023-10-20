import subprocess
from pathlib import Path 
import numpy as np 

conf_path = Path(f'examples/button_alectinib/config_opt.ini')
lam = [None, 0.1, 0.1]
lam_1s = [40] # np.linspace(1, 20, 20).round(2)
cmds = []

for l1 in lam_1s:
    out_dir = Path(f'results/button_alectinib/lam_{l1}_{lam[1]}_{lam[2]}')
    cmd = f'sparrow --config {conf_path} --output-dir {out_dir} --reward-weight {l1} --start-cost-weight {lam[1]} --reaction-weight {lam[2]}'
    cmds.append(cmd)
        
for cmd in cmds: 
    print(cmd)
    subprocess.call(cmd, shell=True)
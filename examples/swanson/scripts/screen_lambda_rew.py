import subprocess 
import numpy as np 

lambda_rews = np.r_[0.01:1.01:0.01]
target_csv = 'examples/swanson/targets.csv'
graph = 'examples/swanson/data/trees_w_costs.json'
output_base = 'results/swanson_vary_rew'
max_rxns = 100

cmd_template = f"sparrow --target-csv {target_csv} --graph {graph} --dont-buy-targets"

for lam in lambda_rews: 
    cmd = f'{cmd_template} --max-rxns {max_rxns} --reward-weight {lam} --reaction-weight {1-lam} --start-cost-weight 0 --output-dir {output_base}/lamrew_{lam:0.2f}'
    subprocess.call(cmd, shell=True)

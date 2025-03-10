import subprocess 
import numpy as np 

lambda_rews = np.r_[0.01:1.01:0.01]
target_csv = 'examples/automoldesigner/targets.csv'
graph = 'examples/automoldesigner/chkpts/trees_w_costs.json'
output_base = 'results/automoldesigner_vary_rew'
max_rxns = 500

cmd_template = f"sparrow --target-csv {target_csv} --graph {graph} --dont-buy-targets"

for lam in lambda_rews: 
    cmd = f'{cmd_template} --max-rxns {max_rxns} --reward-weight {lam} --reaction-weight {1-lam} --start-cost-weight 0 --output-dir {output_base}/lamrew_{lam:0.2f}'
    subprocess.call(cmd, shell=True)

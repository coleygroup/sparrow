import subprocess 

base_cmd = 'python examples/li_covid_example/add_tree.py --base-dir examples/li_covid_example --target-dict li_rewards_copy.csv --storage tree_lazy.json --sec 60 '

for i in range(1): 
    cmd = f'{base_cmd} --n {i}'
    print(cmd)
    subprocess.call(cmd, shell=True)
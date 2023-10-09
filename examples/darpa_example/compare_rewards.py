import json 
from pathlib import Path 
import pandas as pd 
import seaborn as sns


base_dir = Path('examples/darpa_example')

with open(base_dir/'routes_34tars.json','r') as f: 
    routes_34 = json.load(f)

rewards_34 = [entry['Reward'] for entry in routes_34.values()]

df = pd.read_csv(base_dir/'cleaned_tar_dict.csv')
rewards_all = list(df['Reward'])

classes = ['Dataset']*len(rewards_all) + ['SPARROW']*len(rewards_34)
rew = rewards_all + rewards_34
df = pd.DataFrame({'Set': classes, 'Rewards': rew})

ax = sns.histplot(data=df, x="Rewards", hue="Set")
fig = ax.get_figure()
fig.savefig(base_dir/"rew.png") 
import json
import yaml
import pandas as pd
from pathlib import Path

 
properties_path = Path('amd_darpa/second_sulfur_explore_scaffold_round4_compiled_chemprop_predictions.yaml')

def obj(logP, photodeg, uv_vis, logP_unc=0, photodeg_unc=0, uv_vis_unc=0): 
    wfs = [1/20, 1/6, 1/600]
    scalar = -wfs[0]*logP + -wfs[1]*photodeg + wfs[2]*uv_vis 
    return scalar


with open(properties_path, 'r') as file:
    properties = yaml.safe_load(file)

objs = []
for smi, props in properties.items(): 
    print(smi)
    objs.append(obj(
        logP = props['logp']['predicted_value'],
        photodeg = props['photodeg']['predicted_value'],
        uv_vis = props['uv_vis']['predicted_value'],
    ))

df = pd.DataFrame(data={'smiles': list(properties.keys()), 'utility': objs})
df = df.sort_values(by=['utility'], ascending=False)
df.to_csv((properties_path.parent / 'target_list.csv'), index=False)
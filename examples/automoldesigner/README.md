# AutoMolDesigner case study for demonstrating SPARROW
Code and data to reproduce the results of applying SPARROW to downselect from putative _Staphylococcus aureus_ proposed by AutoMolDesigner. 

Reference: Shen, T.; Guo, J.; Han, Z.; Zhang, G.; Liu, Q.; Si, X.; Wang, D.; Wu, S.; Xia, J. AutoMolDesigner for Antibiotic Discovery: An AI-Based Open-Source Software for Automated Design of Small-Molecule Antibiotics. _J. Chem. Inf. Model._ 2024, 64 (3), 575–583. https://doi.org/10.1021/acs.jcim.3c01562.

## Data 
The [enamine_per_g_122024.csv](data/enamine_per_g_122024.csv) file contains the SMILES for compounds in Enamine's Building Block database as of December 2024. [staph_targets_raw.csv](data/staph_targets_raw.csv) was obtained directly from the AutoMolDesigner paper's Supplementary Information. 

## Obtaining Results 
The general workflow for this case study follows these steps, using the specified file in the [scripts folder](scripts): 
1. Perform retrosynthesis with ASKCOS for all candidates ([01_v2_treebuild.py](scripts/01_v2_treebuild.py))
2. Obtain conditions for all reactions in the predicted routes ([02_conditions.py](scripts/02_conditions.py))
3. Obtain likelihood of reaction success scores for all reactions in the predicted routes, incorporating the previously predicted conditions ([03_scores.py](scripts/03_scores.py))
4. Generate a json file with all reactions that have a likelihood score greater than 0.01 ([04_route_graph.py](scripts/04_route_graph.py))
5. Eliminate any node _not_ within a distance of 16 from at least one candidate and assign buyability according to Enamine's building block database ([05_costs_and_prune.py](scripts/05_costs_and_prune.py))
6. Perform the SPARROW runs defined in [cmd_list.txt](scripts/cmd_list.txt) ([run_cmds.batch](scripts/run_cmds.batch) or [run_cmds.py](scripts/run_cmds.py))

All data from steps 1-5 is included in the [trees_w_costs.json.gz](data/trees_w_costs.json.gz) file. Simply decompress this file and begin running SPARROW, after following the installation instructions in [main README](../../README.md). 

To obtain results that depend on reaction classes, use NameRxn to classify all reactions and store results in a csv file following [this template](../templates/reaction_classes.csv). The commands in [cmd_list.txt](scripts/cmd_list.txt) assume this file is named ``reaction_classes.csv`` and is in the [data](data) directory. 

## Analyzing results
The [figures](figures/figures.ipynb) notebook contains code to generate figures and analyze results, assuming that the output directories in [cmd_list.txt](scripts/cmd_list.txt) are unchanged and that all runs are called from the root directory of the SPARROW repository. 


# SPARROW (Synthesis Planning And Rewards-based Route Optimization Workflow)

A workflow to simultaneously select molecules and their synthetic routes for lead optimization and design-make-test cycles. This optimization approach aims to minimize synthetic cost while selecting the molecules that are most likely to fulfill design constraints. More details about SPARROW can be found in [our preprint](https://arxiv.org/abs/2311.02187).   

## Overview 
This repository performs the following steps: 
1. Performs a tree search to identify synthetic routes for a set of provided candidate molecules using ASKCOS 
2. Combines synthetic routes for all molecules into one synthesis network (defined as a **RouteGraph** object).
3. Calculates confidence scores for all reactions in the network using the ASKCOS forward predictor. These scores indicate the confidence of the model in the success of the reaction and serve as a proxy to estimate the likelihood that the reaction will succeed. 
4. Defines optimization problem variables using PuLP
5. Sets relevant constraints on the optimization variables and sets the objective function. 
6. Solves optimization problem.  
7. Outputs the resulting selected routes and target molecules. 

## Table of Contents 
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running SPARROW](#running-sparrow)
- [Reproducing Results](#reproducing-results)
- [Future Goals](#future-goals)

## Requirements 
To use ASKCOS to perform retrosynthesis searches, propose conditions, and score reactions, an API address to a deployed version of ASKCOS is required. Directions to deploy ASKCOS can be found [here](https://github.com/ASKCOS/ASKCOS). A ChemSpace API key is required use ChemSpace to assign compound buyability and cost. Refer to [these directions](https://api.chem-space.com/docs/#:~:text=Get%20API%20Key&text=The%20API%20key%20is%20unique,%40chem%2Dspace.com.) to attain the key. The key should be entered into the [keys.py](keys.py) file, and the path to this file is required to run SPARROW with ChemSpace. 

Dependencies: 
* Python >= 3.7
* Remaining requirements in [requirements.txt](requirements.txt)

SPARROW has been tested on Python version 3.7 and 3.8 and on Linux and Windows machines. 

## Installation

On a standard desktop computer, installation of SPARROW should typically take less than one hour. To begin installation, create a conda environment for SPARROW using [mamba](https://mamba.readthedocs.io/en/latest/installation.html) and install additional requirements through pip. If you are using Gurobi, obtain a Gurobi license and follow the instructions below to use with SPARROW. 

```
mamba env create -f environment.yml
conda activate sparrow
pip install -r requirements.txt
```
Finally, install this package. 
```
python setup.py develop 
```

#### Installing Mamba
If on Linux with x86 cores: 
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
Otherwise, use the correct link from [Mambaforge installation page](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh). 

#### Installing and using Gurobi 
If you intend to use Gurobi to solve the optimization, first ensure that you have [obtained a Gurobi license](https://www.gurobi.com/solutions/licensing/). Then, in the relevant conda environment, install your Gurobi license using the following: 
```
grbgetkey <your key>
```
You can check the status and expiration date of your license using `gurobi_cl --license`. 

## Running SPARROW
The general command to run SPARROW is:
`sparrow --target-csv <path/to/target_csv> --path-finder {api, lookup} --recommender {api, lookup} --coster {naive, chemspace} [additional arguments]`

Alternatively, you may run the following command: 
`sparrow --config <path/to/config>`
which directs SPARROW to a configuration file with all arguments. Some example config files are in the [examples folder](examples).

If more than one run is performed on the same set of compounds, and only the weighting factors are changed, the route graph information from the first run can be easily incorporated into following runs. By providing the path to the `trees_w_info.json` file from the initial run as the argument for `--graph`, all potential paths, compound costs, and reaction scores are incorporated. In these cases, entries for `--path-finder`, `--recommender`, `--scorer`, and `--coster` are not required. 

#### Settings 
All arguments and descriptions can be viewed by running `sparrow --help`. Below we include the arguments we expect to be most relevant to users. 
 - `--config`: the filepath of the configuration file
 - `--formulation {linear, expected_reward}` (default: `linear`): which formulation to use. `linear` is the formulation optimized in our preprint, while `expected_reward` directly optimizes expected reward (see preprint, Figure 2)
 - `--target-csv`: the filepath of the target csv file
 - `--output-dir`: where to save checkpoints files and paths generated by SPARROW
 - `--graph`: path to route graph json file. If provided, no route planning is performed
 - `--reward-weight` <sup>L</sup>: weighting factor for reward objective 
 - `--start-cost-weight` <sup>L</sup>: weighting factor for starting material cost objective
 - `--reaction-weight` <sup>L</sup>: weighting factor for reaction objective
 - `--max-rxns` <sup>E</sup>: maximum number of reaction steps to select
 - `--starting-material-budget` <sup>E</sup>: maximum starting material cost 
 - `--cluster {custom, similarity}` (default: `None`): How to define clusters if desired. If `custom`, they must be defined in the `target-csv` file. If `similarity`, they are automatically defined using Tanimoto similarity-based clusters with a maximum cutoff of `--cluster-cutoff` (default: `0.7`).
 - `--N-per-cluster`: Constrains that all clusters are represented by at least `N` compounds 
 - `--diversity-weight` (default: `0`) <sup>L</sup>: A weighting factor that encourages selection of many similarity-based clusters. 
 - `--path-finder {lookup,api}`: type of tree builder to use
 - `--tree-lookup-dir`: path of lookup json file with combined retrosynthesis tree
 - `--time-per-target`: expansion time in seconds for each target
 - `--max-ppg`: maximum price per gram in dollars for starting materials for ASKCOS MCTS tree search
 - `--max-branching`: maximum branch factor for ASKCOS MCTS tree search
 - `--tree-host`: host address for tree builder, if using ASKCOS API path finder
 - `--recommender {lookup,api}`: type of context recommender to use
 - `--context-host`: host address for context recommender, if using API recommender
 - `--context-lookup`: path of lookup csv file for lookup context recommender
 - `--scorer {lookup,api}`: type of scorer to use
 - `--scorer-host`: host address for reaction scorer, if using API recommender
 - `--scorer-lookup`: path of reaction scorer csv file for lookup reaction scorer (not implemented yet)
 - `--coster {lookup, naive, chemspace}`: type of compound coster to use
 - `--key-path` path that includes the file keys.py with chemspace api key
 - `--inventory`: path to CSV file with SMILES strings and costs, if using lookup coster
 - `--dont-buy-targets`: ensures that solution does not propose directly buying any target
 - `--solver {pulp, gurobi}` (default: `pulp`) <sup>L</sup>: the solver to use. The expected reward formulation requires and defaults to Gurobi. 

<sup>L</sup> only used in the linear formulation 
<sup>E</sup> only used in the expected reward formulation

**A note about required arguments:** The only required argument in SPARROW in `--target-csv`. However, providing this alone will not be sufficient to run SPARROW. In addition to candidates and rewards, SPARROW's optimization requires a set of potential reactions and scores for each reaction. If a provided `--graph` argument corresponds to a file that includes both potential reactions as a retrosynthesis tree _and_ reaction scores, that is sufficient to run SPARROW. However, if the file only contains a retrosynthesis tree, without reaction scores, SPARROW will require a `--scorer` argument. Likewise, if no `--graph` is provided, a valid entry for `--path-finder` (and any corresponding arguments) are required. We are currently working on expanding the documentation for SPARROW and improving its usability.

#### Outputs
SPARROW outputs a set of results files, and checkpoints if relevant, to the specified output directory (`--output-dir`). All parameters used for that run are output in `params.ini`. A summary of SPARROW's output, including the number of candidates, reactions, and starting materials selected, is included in `summary.json`. The selected routes are provided in two separate formats, both json files. `routes.json` provides the synthetic route to each selected candidate, and `solution_list_format.json` individually lists the selected candidates, selected reactions, and selected buyable materials. A pruned target SMILES/rewards csv file will also be output by SPARROW _if_ at least one target cannot be found in the retrosynthetic graph.

#### Run time 
The run time associated with SPARROW depends on the information provided. If a route graph is provided, SPARROW will typically run in seconds (< ~100 candidates) or minutes (> ~100 candidates). When retrosynthesis, condition recommendation, reaction scoring, and compound buyability must all be performed, the total run time for SPARROW will scale linearly with the number of candidates. For the [example of 300 candidate molecules](examples/button_alectinib/), the total runtime of the SPARROW workflow was approximately 13 hours. Approximately 5 hours of retrosynthesis planning, 4 hours of searching for buyability and cost, and 4 hours of condition recommendation and scoring contributed to this computation cost.

##  Optimization Problem Formulation 
The formulation of the optimization problem can be found at our [preprint](https://arxiv.org/abs/2311.02187). 

## Reproducing Results
The results shown in our [preprint](https://arxiv.org/abs/2311.02187) can be reproduced using the [optimize_preprint](scripts/optimize_preprint.py) script. This uses SPARROW to select routes from previously generated retrosynthesis trees with previously computed conditions and reaction scores. For each case study, this information is stored in the relevant [examples folder](examples) as a `trees_w_info.json` file. Each example contains a `config_opt.ini` that builds a route graph from the existing `trees_w_info.json` file and (if applicable) a `config.ini` that builds a route graph from scratch using ASKCOS. In order to use the sample `config.ini` files, you must enter an IP address corresponding to an ASKCOS instance where indicated and enter a ChemSpace API key in [keys.py](keys.py). 

## Future Goals
1. Incorporate reaction conditions into objective function 
2. Update tree visualization
3. Modify optimization function formulation to better capture information gain and cost 




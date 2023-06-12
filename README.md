# SPARROW (Synthesis Planning And Rewards-based Route Optimization Workflow)

A workflow to simultaneously select molecules and their synthetic routes for lead optimization and design-make-test cycles. This optimization approach aims to minimize synthetic cost while selecting the molecules that are most likely to fulfill design constraints.  

## Overview 
This repository performs the following steps: 
1. Performs a tree search to identify synthetic routes for a set of provided candidate molecules using ASKCOS 
2. Combines synthetic routes for all molecules into one synthesis network (defined as a **RouteGraph** object).
3. Calculates confidence scores for all reactions in the network using the ASKCOS forward predictor. These scores indicate the confidence of the model in the success of the reaction and serve as a proxy to estimate the likelihood that the reaction will succeed. 
4. Defines optimization problem variables using PuLP
5. Sets relevant constraints on the optimization variables and sets the objective function. 
6. Solves optimization problem with Gurobi.  
7. Visualizes the resulting selected routes and target molecules. 

## Install 
Create conda environment using [mamba](https://mamba.readthedocs.io/en/latest/installation.html) and install additional requirements through pip. 

```
mamba env create -f environment.yml
conda activate sparrow
pip install -r requirements.txt
```

Install ASKCOS to the local repository. SPARROW can still be used without ASKCOS if a **RouteGraph** with previously calculated scores is available.  
```
git clone *link to askcos-core*
```

You also need to add askcos-data to the askcos path to use the models and ensure that the large files are also pulled. (git lfs must be installed.)
```
cd askcos-core/askcos/
git clone *link to askcos-data*
cd askcos-data/
git lfs pull
```
Then change the name of askcos-data folder to "data". 

Finally, setup this package: 
```
python setup.py develop
```
Note: currently trying to automatically install askcos when running setup.py, but it isn't working. For now, all scripts should start by adding the askcos-core folder to path so that askcos can be imported.  

### Installing Mamba
If on Linux with x86 cores: 
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
Otherwise, use correct link from [Mambaforge installation page](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh)

### ASKCOS Location 
A publically available version is available [here](https://github.com/ASKCOS). Members of the [MLPDS](https://mlpds.mit.edu/) Consortium have access to the most recent version of ASKCOS on GitLab. 

## Requirements 
Currently requires ASKCOS to be downloaded in the directory. However, a quick tutorial on the optimization part is provided that does not require ASKCOS. A set of routes and reaction scores has already been calculated and saved. 


## Running optimization workflow 
TODO: describe how to complete a run 

##  Optimization Problem Formulation 
TODO: describe optimization problem, constraints, and objective function 

## TO DO list
1. Search for starting materials through ChemSpace API instead of naive definition 
2. Get starting material cost from ChemSpace and include in objective function 
3. Add reaction conditions and incorporate into constraints/objective function 
4. Update tree visualization to show compounds and conditions (not just IDs)
5. Add CLI functionality 

## Troubleshooting 

**Recursion error during pickling** 
For route_graphs with many (> ~1000) reactions, you may run into a an recursion error when trying to save the object. As stated [here](http://docs.python.org/library/pickle.html#what-can-be-pickled-and-unpickled), you fix this error by raising the recursion limit (default is 1000): 
```
import sys 
sys.setrecursionlimit(2000)
```

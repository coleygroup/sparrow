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
Currently requires ASKCOS to be downloaded in the directory. (TODO: allow ASKCOS to not be in directory if RouteGraph already defined with scores.)

## Running optimization workflow 
TODO: describe how to complete a run 

##  Optimization Problem Formulation 
TODO: describe optimization problem, constraints, and objective function 

## TO DO list
1. Complete tutorials, both that require ASKCOS and that do not require ASKCOS 
2. Search for starting materials through ChemSpace API instead of naive definition 
3. Get starting material cost from ChemSpace and include in objective function 
4. Add reaction conditions and incorporate into constraints/objective function 
5. Update tree visualization to show compounds and conditions (not just IDs)
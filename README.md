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
Create conda environment and install this module. 

```
conda env create -f environment.yml (Won't work yet)
conda activate sparrow
pip install -r requirements.txt
python setup.py develop
```

## Requirements 
Currently requires ASKCOS to be downloaded in the directory. (TODO: allow ASKCOS to not be in directory if RouteGraph already defined with scores.)

## Running optimization workflow 
TODO: describe how to complete a run 

##  Optimization Problem Formulation 
TODO: describe optimization problem, constraints, and objective function 


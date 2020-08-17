# Molecule_library_synthesis
This is the code repository for optimization of synthesis plans for multiple molecules.
The code will require external dependency of the ASKCOS software, which need to be downloaded at
git clone https://github.com/Coughy1991/ASKCOS.

After downloading the repository, the folder ASKCOS as well as the subfolder ASKCOS/askcos should be added to PYTHONPATH

Another important dependency for running the code is the GUROBI optimizer. The instruction for installation and obtaining a free academic license can be found at https://www.gurobi.com/

The pipeline for running the code include the following steps, using case study 1 as an example:
1. run KERAS_BACKEND=theano python case_1/tree_builder.py to construct the reaction network for multiple molecule targets
2. run KERAS_BACKEND=theano python predict_conditions_and_evaluate.py case_1 to evaluate the reactions
3. run jupyter notebook and run through the blocks of code to optimize pathway selection


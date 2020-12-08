# Alchemical_FE

The alchemical_fe folder contains the following subfolders:
    1. Tools – containing scripts to help set up TI simulations
    2. Examples – containing test cases for relative binding free energy (rbfe) and relative solvation free energy (rsfe) calculations.
    3. Documentation – containing documentation specific to AMBER20_DD_BOOST 
    4. Tutorials – containing tutorials for using AMBER20_DD_BOOST


The afe_setup_clean.sh script, located in the Tools directory, is designed to help setup the directory structure and input files that are necessary for running alchemical free energy simulations with AMBER20_DD_BOOST. 
afe_setup_clean.sh accepts a simplified input file named ‘input’ that is described in detail in the tutorial located in the Tutorials directory. 
Briefly, for a given system, such as a specific protein target, or a collection of small molecules, 
	a) a list of desired transformations can be provided, 
	b) key simulation settings can be specified, and 
	c) initial configuration files (MD equilibrated parameter (parm) and coordinate (rst) files associated with the system and specified transformations) must be provided. 

afe_setup_clean.sh can then be used to generate a hierarchy of directories containing relevant parameter, coordinate, and AMBER input files and job submission scripts. 
afe_setup_clean.sh can subsequently be used to run and check all equilibration and production free energy simulations of that particular system

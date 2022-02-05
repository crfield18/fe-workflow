# Alchemical_FE

The alchemical_fe folder contains the following subfolders:
    1. bin – containing workflow tools to help set up and analyze TI simulations
    2. Examples – containing test cases for relative binding free energy (rbfe) and relative solvation free energy (rsfe) calculations.
    3. Documentation – containing documentation specific to AMBER20_DD_BOOST 
    4. Tutorials – containing tutorials for using AMBER20_DD_BOOST


Workflow Tools are a set of scripts that are designed to facilitate the setup, execution, and analysis of alchemical
free energy (AFE) simulations using AMBER DD Boost. Currently, the Workflow Tools can be used to perform relative binding
free energy (RBFE), relative solvation free energy (RSFE), and absolute solvation free energy (RSFE) calculations. The scripts
use a simplified input file, which is described in detail later in the user-guide, that provides top-level control on various important
aspects of the intended AFE simulations.

Briefly, for a given system, such as a specific protein target, or a collection of small molecules,
a) a list of desired transformations can be provided,
b) key simulation settings can be specified, and
c) initial configuration files (MD equilibrated parameter (parm) and coordinate (rst) files associated with the system and
specified transformations) must be provided.

The Workflow Tools can then be used to generate a hierarchy of directories containing relevant parameter, coordinate, and AMBER 
input files and job submission scripts. The Workflow Tools can also be used to analyze the free energy simulations using FE-ToolKit.


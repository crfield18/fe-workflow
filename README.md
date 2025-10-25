# FE-Workflow

This is a forked copy of the original [GitLab repo](https://gitlab.com/RutgersLBSR/FE-Workflow) by [@pgbarletta](https://github.com/pgbarletta)!

This free energy workflow package provides a collection of programs that simplify the setup and analysis of alchemical free energy simulations such as relative binding free energy (RBFE), relative solvation free energy (RSFE), and absolute solvation free energy (ASFE) to be run with Amber24.

Some introductory information and examples are available in the links below.

 - [Setting up alchemical free energy simulations](https://ambertutorials-rutgerslbsr-c744272d5a9c1169e0dc9e19b8d800019105.gitlab.io/workshop/05_afeSims/03_asfe.html)
 - [Analyzing alchemical free energy simulations](https://ambertutorials-rutgerslbsr-c744272d5a9c1169e0dc9e19b8d800019105.gitlab.io/workshop/05_afeSims/04_rbfe.html)


## Requirements:

 - AmberTools.
 - [FE-ToolKit](https://gitlab.com/RutgersLBSR/fe-toolkit).

## Installation

AMBER should be installed and `AMBERHOME` defined. AmberTools' `cpptraj`, `parmed`, `edgembar`, and `fetkutils` are used by `FE-Workflow`.
`bin/` subdirectories within `FE-Workflow` should be in the `$PATH`.

The `setup_fe` script represents the main executable of `FE-Workflow`, and can be created by running the script `makesetup_fe.sh` located in the `FE-Workflow` repository.
`makesetup_fe.sh` will also write a `setup_directives` file with 3 paths: `MDEngine` (your `$AMBERHOME`) `ToolKit` (`$AMBERHOME`/bin) and `Workflow` (the `FE-Workflow` path), and the `FE-Workflow.bashrc` file that should be sourced before executing `setup_fe`.


In the default setup, `setup_fe` is meant to be kept in the `FE-Workflow/bin` directory and, since `FE-Workflow/bin` is added to the `$PATH` variable, `setup_fe` should be available as a command line program.
In most cases, this should be the most convenient way of using `setup_fe`. However, if needed, the location of `setup_fe` can be changed by changing the `$PATH` variable in `setup_fe` accordingly.

### Environment settings on Linux

If you have both `ambertools` and a standalone version of `fe-toolkit` installed (from source or pip), you must ensure the amber is not overriding your standalone installation.
This conflict can prevent FE-Workflow from finding the correct fe-toolkit version.
To fix this, you can add the following to your `.bashrc`:

#### if you installed `fe-toolkit` from source

```
export BACKUP_PATH=$PATH
export BACKUP_PYTHONPATH=${PYTHONPATH}
source <amber installation path>/amber.sh
export PATH=$BACKUP_PATH:$PATH
export PYTHONPATH=${BACKUP_PYTHONPATH}:$PYTHONPATH
export PATH=$PATH:<FE-Workflow path>/bin
```

This will effectively put the amber binaries and python libraries at the end of each path.

#### if you installed `fe-toolkit` through pip

```
source <amber installation path>/amber.sh
export PATH=<fe-toolkit installation path>/local/bin:$PATH
export PYTHONPATH=<fe-toolkit installation path>/local/lib/python3.11/site-packages:${PYTHONPATH}
export PATH=$PATH:<FE-Workflow path>/bin
```

## References:

1. Tai-Sung Lee, Hsu-Chun Tsai, Abir Ganguly, Timothy J. Giese, and Darrin M. York. Robust, Efficient
and Automated Methods for Accurate Prediction of Protein-Ligand Binding Affinities in AMBER Drug
Discovery Boost, volume 1397 of ACS Symposium Series. November 2021. ISBN 9780841298064. doi:
10.1021/bk-2021-1397.ch007
2. FE-ToolKit: A Versatile Software Suite for Analysis of High-Dimensional Free Energy Surfaces and Alchemical Free Energy Networks Timothy J. Giese, Ryan Snyder, Zeke Piskulich, German P. Barletta, Shi Zhang, Erika McCarthy, Solen Ekesan, and Darrin M. York J. Chem. Inf. Model (2025). doi: 10.1021/acs.jcim.5c00554
3. Zoe Cournia, Christophe Chipot, Benoît Roux, Darrin M. York, and Woody Sherman. Free Energy
Methods in Drug Discovery—Introduction, volume 1397 of ACS Symposium Series. November 2021.
ISBN 9780841298064. doi: 10.1021/bk-2021-1397.ch001.
4. Timothy J. Giese and Darrin M. York. Variational Method for Networkwide Analysis of Relative
Ligand Binding Free Energies with Loop Closure and Experimental Constraints. J. Chem. Theory
Comput., 17(3):1326–1336, 2021
5. Tai-Sung Lee, Bryce K. Allen, Timothy J. Giese, Zhenyu Guo, Pengfei Li, Charles Lin, T. Dwight McGee
Jr., David A. Pearlman, Brian K. Radak, Yujun Tao, Hsu-Chun Tsai, Huafeng Xu, Woody Sherman,
and Darrin M. York. Alchemical Binding Free Energy Calculations in AMBER20: Advances and Best
Practices for Drug Discovery. J. Chem. Inf. Model., 60:5595–5623, 2020
6. Tai-Sung Lee, Zhixiong Lin, Bryce K. Allen, Charles Lin, Brian K. Radak, Yujun Tao, Hsu-Chun
Tsai, Woody Sherman, and Darrin M. York. Improved Alchemical Free Energy Calculations with
Optimized Smoothstep Softcore Potentials. J. Chem. Theory Comput., 16:5512–5525, September 2020. ISSN 1549-9626. doi: 10.1021/acs.jctc.0c00237
7. Hsu-Chun Tsai, Yujun Tao, Tai-Sung Lee, Kenneth M. Merz, and Darrin M. York. Validation of Free
Energy Methods in AMBER. J. Chem. Inf. Model., 60:5296–5300, November 2020. ISSN 1549-960X.
doi: 10.1021/acs.jcim.0c00285

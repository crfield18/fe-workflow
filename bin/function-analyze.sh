#!/bin/bash
function write_edgembar {

        ticalc=$1; shift
        cat << EOF > DiscoverEdges.py
#!/usr/bin/env python3

import edgembar
import os
from pathlib import Path

#
# The output directory (where the edge xml input files are to be written
#
odir = Path("analysis")

#
# The format string describing the directory structure.
# The {edge} {env} {stage} {trial} {traj} {ene} placeholders are used
# to extract substrings from the path; only the {edge} {traj} and {ene}
# are absolutely required.  If the {env} placeholder is missing, then
# 'target' evironment is assumed.
#
# Full example:
#    s = r"dats/{trial}/free_energy/{edge}_ambest/{env}/{stage}/efep_{traj}_{ene}.dat"
# Minimal example:
#    s = r"dats/{edge}/efep_{traj}_{ene}.dat"

EOF

        if [ "${ticalc}" == "rbfe" ]; then
        cat << EOF >> DiscoverEdges.py
s = r"data/{edge}/{env}/{trial}/efep_{traj}_{ene}.dat"
exclusions=None
edges = edgembar.DiscoverEdges(s,exclude_trials=exclusions,
                               target="com",
                               reference="aq" )
EOF
        elif [ "${ticalc}" == "rsfe" ]; then
        cat << EOF >> DiscoverEdges.py
s = r"data/{edge}/{env}/{trial}/efep_{traj}_{ene}.dat"
exclusions=None
edges = edgembar.DiscoverEdges(s,exclude_trials=exclusions,
                               target="com",
                               reference="aq" )
EOF

        elif [ "${ticalc}" == "asfe" ]; then
        cat << EOF >> DiscoverEdges.py
s = r"data/{edge}/{env}/{trial}/efep_{traj}_{ene}.dat"
exclusions=None
edges = edgembar.DiscoverEdges(s,exclude_trials=exclusions,
                               target="aq",
                               reference="vac")
EOF
        fi
        cat << EOF >> DiscoverEdges.py

#
# In some instances, one may have computed a stage with lambda values
# going in reverse order relative to the thermodynamic path that leads
# from the reactants to the products. We can reverse the order of the
# files to effectively negate the free energy of each state (essentially
# treating the lambda 0 state as the lambda 1 state).
#
#for edge in edges:
#    for trial in edge.GetAllTrials():
#        if trial.stage.name == "STAGE":
#            trial.reverse()

if not odir.is_dir():
    os.makedirs(odir)

for edge in edges:
    fname = odir / (edge.name + ".xml")
    edge.WriteXml( fname )
                                                                                      
EOF
chmod a+x DiscoverEdges.py

}

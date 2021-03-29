list=(1h1q_com_dry)
rbuf=15
        cat <<EOF > tleap.in
# load the AMBER force fields
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip4pew
loadamberparams frcmod.ff14SB
loadamberparams frcmod.tip4pew
loadAmberParams frcmod.ionsjc_tip4pew
loadoff tip4pewbox.off


# load force field parameters for ligand
loadamberparams ${list[0]}.frcmod
loadoff ${list[0]}.lib

# load the coordinates and create the complex
complex = loadPdb ${list[0]}.pdb

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

#savepdb complex ${list[0]}-solv.pdb
saveamberparm complex ${list[0]}-solv.parm7 ${list[0]}-solv.rst7

quit
EOF
        tleap -s -f tleap.in > output

#complex = loadPdbUsingSeq ${list[0]}.pdb { $(cat ${list[0]}.seq) }

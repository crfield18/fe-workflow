#!/bin/bash

if [ ! -e "methane.frcmod" ]; then
   parmchk2 -i methane.mol2 -o methane.frcmod -f mol2
fi

if [ -e "methane.frcmod" ]; then
if [ "$(grep -l ATTN methane.frcmod)" != "" ]; then
   echo "The frcmod file "methane.frcmod" contains unset parameters"
   echo "You need to supply a properly parameterized frcmod file to timutate"
   echo "or you'll need to modify ticopy.sh"
   exit 1  
fi
fi

if [ ! -e "methane.lib" ]; then
   cat <<EOF > ticopy.lib.inp
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff
source leaprc.water.tip4pew
methane = loadmol2 methane.mol2
saveoff methane methane.lib
quit
EOF
   tleap -s -f ticopy.lib.inp
   rm leap.log
   rm ticopy.lib.inp
fi

cat << 'EOF' > ticopy.sh.cmds
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff
source leaprc.water.tip4pew
loadOff ticopy.lib
loadAmberParams ticopy.frcmod
loadOff methane.lib
loadAmberParams methane.frcmod

mol0 = loadPdbUsingSeq ticopy.mol0.pdb { L1 }
mol1 = loadPdbUsingSeq ticopy.mol1.pdb { methane }

solvent = loadPdb ticopy.solvent.pdb



x = combine { mol0 mol1 solvent }

setbox x centers

saveAmberParm x ticopy.parm7 ticopy.rst7
quit
EOF

tleap -s  -f ticopy.sh.cmds | grep -v "+---" | grep -v "+Currently" > ticopy.sh.out
cat ticopy.sh.out

if [ ! -e ticopy.parm7 -o ! -e ticopy.rst7 ]; then echo "Failed to make parameter file"; exit 1; fi

# Set the box dimensions
ChBox -X 32.220537000000 -Y 31.784920000000 -Z 31.349282000000 -al 90.000000000000 -bt 90.000000000000 -gm 90.000000000000 -c ticopy.rst7 -o ticopy.rst7.tmp; mv ticopy.rst7.tmp ticopy.rst7



################################


parmutils-tigen.py --uniti -p ticopy.parm7 -c ticopy.rst7 \
--molmask1="@1-6" --timask1="@1-6" --scmask1="@6" --molmask2="@7-11" --timask2="@7-11" --scmask2="" \
--nlambda=21 \
--nlambda-softcore=-1 \
--nmropt=0 \
--cut=10.00 \
--nstlim=1000000 \
--numexchg=0 \
--ntpr=20000 \
--ntwx=20000 \
--cpu-partition="main" \
--cpus-per-node=24 \
--gpu-partition="gpu" \
--gpus-per-node=4 \
--qos="" \
--account="" \
--exclude="" \
--constraint="" \
--min-nodes=-1 \
--max-nodes=-1 \
--days=2 \
--launch="srun --mpi=pmi2" \
--steps=1 \
--steps-per-slurm=1 \


#POST


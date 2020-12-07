#!/bin/bash

if [ ! -e "ejm55.frcmod" ]; then
   parmchk2 -i ejm55.mol2 -o ejm55.frcmod -f mol2
fi

if [ -e "ejm55.frcmod" ]; then
if [ "$(grep -l ATTN ejm55.frcmod)" != "" ]; then
   echo "The frcmod file "ejm55.frcmod" contains unset parameters"
   echo "You need to supply a properly parameterized frcmod file to timutate"
   echo "or you'll need to modify ticopy.sh"
   exit 1  
fi
fi

if [ ! -e "ejm55.lib" ]; then
   cat <<EOF > ticopy.lib.inp
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff
source leaprc.water.tip4pew
ejm55 = loadmol2 ejm55.mol2
saveoff ejm55 ejm55.lib
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
loadOff ejm55.lib
loadAmberParams ejm55.frcmod

mol0 = loadPdbUsingSeq ticopy.mol0.pdb { L1 }
mol1 = loadPdbUsingSeq ticopy.mol1.pdb { ejm55 }

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
ChBox -X 40.552674800000 -Y 40.552674800000 -Z 40.552674800000 -al 109.471219000000 -bt 109.471219000000 -gm 109.471219000000 -c ticopy.rst7 -o ticopy.rst7.tmp; mv ticopy.rst7.tmp ticopy.rst7

# Reset the ifbox flag
sed -i 's|0       0       1|0       0       2|' ticopy.parm7




################################


parmutils-tigen.py --uniti -p ticopy.parm7 -c ticopy.rst7 \
--molmask1="@1-35" --timask1="@1-35" --scmask1="@22,32" --molmask2="@36-68" --timask2="@36-68" --scmask2="" \
--nlambda=12 \
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


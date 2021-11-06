function calcnonwaterinparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set Nonwat = residues inmask !:WAT,Na+,Cl-
show Nonwat

quit
EOF
        cpptraj < cpp.in > log
        local nonwat=$(grep "Nonwat =" log|tail -1| awk -F "'" '{print $2}')
        echo $nonwat
        rm -rf cpp.in
}

function calcwaterinparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set Nwat = residues inmask :WAT
show Nwat

quit
EOF
        cpptraj < cpp.in > log
        local nwat=$(grep "Nwat =" log|tail -1| awk -F "'" '{print $2}')
        echo $nwat
        rm -rf cpp.in
}

function calcionsinparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set Nions = residues inmask :Na+
show Nions

quit
EOF
        cpptraj < cpp.in > log
        local nions=$(grep "Nions =" log|tail -1| awk -F "'" '{print $2}')
        echo $nions
        rm -rf cpp.in
}

function calcsodiuminparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set Nions = residues inmask :Na+,K+
show Nions

quit
EOF
        cpptraj < cpp.in > log
        local nions=$(grep "Nions =" log|tail -1| awk -F "'" '{print $2}')
        echo $nions
        rm -rf cpp.in
}

function calcchlorideinparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set Nions = residues inmask :Cl-
show Nions

quit
EOF
        cpptraj < cpp.in > log
        local nions=$(grep "Nions =" log|tail -1| awk -F "'" '{print $2}')
        echo $nions
        rm -rf cpp.in
}


function calctotalresinparm
{
        local parm="$1"
        cat << EOF > cpp.in
parm ${parm}
set allres = residues inmask :*
show allres

quit
EOF
        cpptraj < cpp.in > log
        local allres=$(grep "allres =" log|tail -1| awk -F "'" '{print $2}')
        echo $allres
        rm -rf cpp.in
}


function write_tleap_merged {
        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
	local lig1=$1; shift
        local numnonstd=$1; shift
	local lig2=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
	local addions=$1; shift
	local boxbuild=$1; shift
	local s=$1; shift
        
        truncate -s0 tleap.in

        # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> tleap.in
        fi

        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> tleap.in
	elif [ "${lff}" == "gaff" ]; then
		printf "source leaprc.gaff\n" >> tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> tleap.in
                printf "loadoff tip4pewbox.off\n" >> tleap.in
                boxkey="TIP4PEWBOX"
	elif [ "${wm}" == "tip3p" ]; then
		printf "source leaprc.water.tip3p\n" >> tleap.in
		printf "loadamberparams frcmod.tip3p\n" >> tleap.in
		boxkey="TIP3PBOX"
        fi

        # check and load non-standard residue parameter files
        i=0
	numnonstd=$(($numnonstd+0))
        while [ "$i" -lt "${numnonstd}" ]; do
                printf "loadamberparams ${lig1}_${i}.frcmod\n" >> tleap.in
                printf "loadoff ${lig1}_${i}.lib\n" >> tleap.in
                i=$(($i+1))
        done

	# load ligand 2 parameter files
	printf "loadamberparams ${lig2}_0.frcmod\n" >> tleap.in
	printf "loadoff ${lig2}_0.lib\n" >> tleap.in

        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
        if [ "${load}" == "pdb" ]; then
                printf "x = loadPdb ${inpfile}.pdb\n" >> tleap.in
        elif [ "${load}" == "pdbseq" ]; then
                printf "x = loadPdbUsingSeq ${inpfile}.pdb { $(cat ${inpfile}.seq) }\n" >> tleap.in
        else
                printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in
        fi

        # add S-S cysteine linkkages if present
        if [ -f ${inpfile}_sslinks ] && [ "$(cat ${inpfile}_sslinks | wc -l)" -gt 0 ]; then
                while read line; do
                        IFS=' ' read -ra args <<< $line
                        printf "bond x.${args[0]}.SG x.${args[1]}.SG\n" >> tleap.in
                done < ${inpfile}_sslinks
        fi

        # build box and neutralize with Na+ Cl-
	if [ "${boxbuild}" == 0 ] && [ "${s}" == "com" ]; then
                printf "setbox x vdw \n" >> tleap.in
	else
                printf "${boxcmd} x ${boxkey} ${rbuf}\n" >> tleap.in
                printf "addions x Na+ 0\n" >> tleap.in
                printf "addions x Cl- 0\n" >> tleap.in

                # add additional Na+ Cl- if needed
                if [ "${addions}" != "0" ]; then
                        printf "addions x Na+ ${addions}\n" >> tleap.in
                        printf "addions x Cl- ${addions}\n" >> tleap.in
                fi
        fi

        # save parm and rst with -solv suffix
        printf "saveamberparm x merged.parm7 merged.rst7\n\n" >> tleap.in
        printf "quit\n" >> tleap.in
}



function write_tleap_merged_head {
        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local lig1=$1; shift
        local numnonstd=$1; shift
        local lig2=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local nna=$1; shift
        local nch=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift
	
        truncate -s0 tleap.in

        # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> tleap.in
        fi
	

        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> tleap.in
                printf "loadoff tip4pewbox.off\n" >> tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> tleap.in
                boxkey="TIP3PBOX"
        fi

        # check and load non-standard residue parameter files
        i=0
	numnonstd=$(($numnonstd+0))
        while [ "$i" -lt "${numnonstd}" ]; do
                printf "loadamberparams ${lig1}_${i}.frcmod\n" >> tleap.in
                printf "loadoff ${lig1}_${i}.lib\n" >> tleap.in
                i=$(($i+1))
        done

        # load ligand 2 parameter files
        printf "loadamberparams ${lig2}_0.frcmod\n" >> tleap.in
        printf "loadoff ${lig2}_0.lib\n" >> tleap.in

        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
        if [ "${load}" == "pdb" ]; then
                printf "x = loadPdb ${inpfile}.pdb\n" >> tleap.in
        elif [ "${load}" == "pdbseq" ]; then
                printf "x = loadPdbUsingSeq ${inpfile}.pdb { $(cat ${inpfile}.seq) }\n" >> tleap.in
        else
                printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in
        fi

        # add S-S cysteine linkkages if present
        if [ -f ${inpfile}_sslinks ] && [ "$(cat ${inpfile}_sslinks | wc -l)" -gt 0 ]; then
                while read line; do
                        IFS=' ' read -ra args <<< $line
                        printf "bond x.${args[0]}.SG x.${args[1]}.SG\n" >> tleap.in
                done < ${inpfile}_sslinks
        fi

        # build box and neutralize with Na+ Cl-
        if [ "${boxbuild}" == 0 ] && [ "${s}" == "com" ]; then
                printf "setbox x vdw \n" >> tleap.in
        else
                printf "${boxcmd} x ${boxkey} ${rbuf}\n" >> tleap.in

                # add additional Na+ Cl- if needed
                if [ "${nna}" != "0" ]; then
                        printf "addions x Na+ ${nna}\n" >> tleap.in
		fi
                if [ "${nch}" != "0" ]; then
                        printf "addions x Cl- ${nch}\n" >> tleap.in
                fi
        fi

}


function RunTLEAP
{

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
	local lig1=$1; shift
        local numnonstd=$1; shift
	local lig2=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local ns=$1; shift
        local nc=$1; shift
	local boxbuild=$1; shift
	local s=$1; shift

        if ! [ -z ${1+x} ]; then local fdel=${1}; fi
        if ! [ -z ${2+x} ]; then local ldel=${2}; fi

    	write_tleap_merged_head "${pff}" "${lff}" "${wm}" "${inpfile}" "${lig1}" "${numnonstd}" "${lig2}" "${mdboxshape}" "${rbuf}" "${load}" "${ns}" "${nc}" "${boxbuild}" "${s}"
    	cat <<EOF >> tleap.in

savepdb x tmp.pdb
m = loadpdb tmp.pdb
setbox m centers
EOF

    	if [ ! -z "${fdel}" -a ! -z "${ldel}" ]; then
        	for res in $(seq ${ldel} -1 ${fdel}); do
            		echo "remove m m.${res}" >> tleap.in
        	done
    	fi
	unset fdel ; unset ldel

    	cat <<EOF >> tleap.in
saveamberparm m out.parm7 out.rst7
quit
EOF
    	tleap -s -f tleap.in >> log
}

function fix_solvent {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local lig1=$1; shift
        local numnonstd=$1; shift
        local lig2=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local nsod=$1; shift
        local ncl=$1; shift
        local nwat=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift



        for iter in $(seq 100); do
		RunTLEAP "${pff}" "${lff}" "${wm}" "${inpfile}" "${lig1}" "${numnonstd}" "${lig2}" "${mdboxshape}" "${rbuf}" "${load}" "${nsod}" "${ncl}" "${boxbuild}" "${s}"
                current_nwat=$(calcwaterinparm "out.parm7")
		current_totres=$(calctotalresinparm "out.parm7")
		current_sodiums=$(calcsodiuminparm "out.parm7"); #current_nions=$((current_nions*2))
		current_chlorides=$(calcchlorideinparm "out.parm7"); #current_nions=$((current_nions*2))
		current_nonwater=$(calcnonwaterinparm "out.parm7")

                excess_waters=$(bc -l <<< "${current_nwat} - ${nwat}")
                if [ "${excess_waters}" -ge 0 ]; then
                        break
                else
                        rbuf=$(bc -l <<< "${rbuf} + 1.")
                fi
        done
        if [ ${excess_waters} -gt 0 ]; then
                lastdel=$(calctotalresinparm "out.parm7")
                firstdel=$(( ${lastdel} - ${excess_waters} + 1 ))
		RunTLEAP "${pff}" "${lff}" "${wm}" "${inpfile}" "${lig1}" "${numnonstd}" "${lig2}" "${mdboxshape}" "${rbuf}" "${load}" "${nsod}" "${ncl}" "${boxbuild}" "${s}" "${firstdel}" "${lastdel}"
        fi

        current_nwat=$(calcwaterinparm "out.parm7")
        #echo "There are ${current_nwat} waters in out.parm7 built from ${inpfile}"
}


##################################################
## ASFE related
function write_tleap_asfe {
        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local addions=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift

        truncate -s0 tleap.in

        # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> tleap.in
        fi

        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> tleap.in
                printf "loadoff tip4pewbox.off\n" >> tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> tleap.in
                printf "loadamberparams frcmod.tip3p\n" >> tleap.in
                boxkey="TIP3PBOX"
        fi



        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
	printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in

        # load ligand parameter files
        printf "loadamberparams ${inpfile}_0.frcmod\n" >> tleap.in
        printf "loadoff ${inpfile}_0.lib\n" >> tleap.in


        # build box and neutralize with Na+ Cl-
        printf "${boxcmd} x ${boxkey} ${rbuf}\n" >> tleap.in
        printf "addions x Na+ 0\n" >> tleap.in
        printf "addions x Cl- 0\n" >> tleap.in

        # add additional Na+ Cl- if needed
        if [ "${addions}" != "0" ]; then
                printf "addions x Na+ ${addions}\n" >> tleap.in
                printf "addions x Cl- ${addions}\n" >> tleap.in
        fi

        # save parm and rst with -solv suffix
        printf "saveamberparm x merged.parm7 merged.rst7\n\n" >> tleap.in
        printf "quit\n" >> tleap.in
}




function write_tleap_head_asfe {
        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local nna=$1; shift
        local nch=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift

        truncate -s0 tleap.in

        # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> tleap.in
        fi


        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> tleap.in
                printf "loadoff tip4pewbox.off\n" >> tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> tleap.in
                boxkey="TIP3PBOX"
        fi


        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

	printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in

        # load ligand parameter files
        printf "loadamberparams ${inpfile}_0.frcmod\n" >> tleap.in
        printf "loadoff ${inpfile}_0.lib\n" >> tleap.in



        # build box and neutralize with Na+ Cl-
        printf "${boxcmd} x ${boxkey} ${rbuf}\n" >> tleap.in

        # add additional Na+ Cl- if needed
        if [ "${nna}" != "0" ]; then
                printf "addions x Na+ ${nna}\n" >> tleap.in
        fi
        if [ "${nch}" != "0" ]; then
                printf "addions x Cl- ${nch}\n" >> tleap.in
        fi

}

function RunTLEAP_asfe
{

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local ns=$1; shift
        local nc=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift

        if ! [ -z ${1+x} ]; then local fdel=${1}; fi
        if ! [ -z ${2+x} ]; then local ldel=${2}; fi

        write_tleap_head_asfe "${pff}" "${lff}" "${wm}" "${inpfile}" "${mdboxshape}" "${rbuf}" "${load}" "${ns}" "${nc}" "${boxbuild}" "${s}"
        cat <<EOF >> tleap.in

savepdb x tmp.pdb
m = loadpdb tmp.pdb
setbox m centers
EOF

        if [ ! -z "${fdel}" -a ! -z "${ldel}" ]; then
                for res in $(seq ${ldel} -1 ${fdel}); do
                        echo "remove m m.${res}" >> tleap.in
                done
        fi
        unset fdel ; unset ldel

        cat <<EOF >> tleap.in
saveamberparm m out.parm7 out.rst7
quit
EOF
        tleap -s -f tleap.in >> log
}

function fix_solvent_asfe {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local inpfile=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local load=$1; shift
        local nsod=$1; shift
        local ncl=$1; shift
        local nwat=$1; shift
        local boxbuild=$1; shift
        local s=$1; shift



        for iter in $(seq 100); do
                RunTLEAP_asfe "${pff}" "${lff}" "${wm}" "${inpfile}" "${mdboxshape}" "${rbuf}" "${load}" "${nsod}" "${ncl}" "${boxbuild}" "${s}"
                current_nwat=$(calcwaterinparm "out.parm7")
                current_totres=$(calctotalresinparm "out.parm7")
                current_sodiums=$(calcsodiuminparm "out.parm7"); #current_nions=$((current_nions*2))
                current_chlorides=$(calcchlorideinparm "out.parm7"); #current_nions=$((current_nions*2))
                current_nonwater=$(calcnonwaterinparm "out.parm7")

                excess_waters=$(bc -l <<< "${current_nwat} - ${nwat}")
                if [ "${excess_waters}" -ge 0 ]; then
                        break
                else
                        rbuf=$(bc -l <<< "${rbuf} + 1.")
                fi
        done
        if [ ${excess_waters} -gt 0 ]; then
                lastdel=$(calctotalresinparm "out.parm7")
                firstdel=$(( ${lastdel} - ${excess_waters} + 1 ))
                RunTLEAP_asfe "${pff}" "${lff}" "${wm}" "${inpfile}" "${mdboxshape}" "${rbuf}" "${load}" "${nsod}" "${ncl}" "${boxbuild}" "${s}" "${firstdel}" "${lastdel}"
        fi

        current_nwat=$(calcwaterinparm "out.parm7")
        #echo "There are ${current_nwat} waters in out.parm7 built from ${inpfile}"
}



##################################################

function calcMDions {

	local nwat=$1
	local ionconc=$2

        # calculate water and number of ions necessary to reach desired ion conc
        cat <<EOF > calcMDions.py
#!/usr/bin/env python

import numpy as np
import sys

density = 55.5

def get_waters():
        n_waters = sys.argv[1]
        try:
                n_waters = int(n_waters)
        except ValueError:
                print '### ERROR not an integer ###'
                n_waters = get_waters()
        return float(n_waters)

def get_ion_conc():
        ion_conc = sys.argv[2]
        try:
                ion_conc = float(ion_conc)
        except ValueError:
                print '### ERROR not a valid ion concentration ###'
                ion_conc = get_ion_conc()
        return ion_conc

n_waters = get_waters()
ion_conc = get_ion_conc()

ions_needed = round((n_waters/density)*ion_conc)

print ("{0:.0f}".format(ions_needed))
EOF
        chmod a+x calcMDions.py
	local nions=$(python2.7 ./calcMDions.py $nwat ${ionconc})
	echo "${nions}"
	rm -rf calcMDions.py
}

function get_index {
	local value=$1; shift
	local arr=("$@")
	for i in "${!arr[@]}"; do
		if [[ "${arr[$i]}" = "${value}" ]]; then echo "${i}"; fi
	done
}


function parse_pdb {
        local pdbfile=$1; shift
        local mapfile=$1; shift
        local ligmask=$1; shift
        local ligcopy=$1; shift
        local orgmol2=$1; shift
        local tarmol2=$1; shift
        local s=$1; shift
        local out=$1; shift
        local calc=$1; shift
        cat << EOF > parsePDB.py
#!/usr/bin/env python3


def CopyPDB( pdb ):
    import copy
    p = copy.copy( pdb )
    p.coordinates = copy.copy( pdb.coordinates )
    return p


def Strip( pdb, mask ):
    p = CopyPDB( pdb )
    p.strip( "%s"%(mask) )
    return p


def Extract( pdb, mask ):
    return Strip( pdb, "!(%s)"%(mask) )

def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]

def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


def Joincom( mol1, mol2, nonligand ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    elif nonligand.coordinates.shape[0] == 0:
        pq=nonligand
    else:
        pq = mol1 + mol2 + nonligand
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates, nonligand.coordinates) )
        return pq

def Joinaq( mol1, mol2 ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    else:
        pq = mol1 + mol2
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates) )
        return pq


if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    from sys import exit

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Takes as input - a PDB file, an atommap file, two ligand masks \"orglig\" and \"tarlig\", name of the two ligand mol2 files as \"orgmol2\" and \"tarmol2\", and output basename. Splits PDB into water.pdb, nonwater.pdb, mol1.pdb containing \"orglig\", and mol2.pdb containing a new residue \"tarlig\" which consists of a subset of \"orglig\" atoms that are mapped as per provided atommap file. Merges the different pdbs into a single pdb and generates the scmasks and timasks" )

    parser.add_argument("-p","--pdb",
                        help="PDB file",
                        type=str,
                        required=True)

    parser.add_argument("-m","--map",
                        help="Atom mapfile",
                        type=str,
                        required=True )

    parser.add_argument("-ol","--orglig",
                        help="name of original ligand in PDB file" ,
                        type=str,
                        required=False)

    parser.add_argument("-tl","--tarlig",
                        help="name of target ligand after mutation" ,
                        type=str,
                        required=False)

    parser.add_argument("-om","--orgmol2",
                        help="name of original ligand mol2 file" ,
                        type=str,
                        required=False)

    parser.add_argument("-tm","--tarmol2",
                        help="name of target ligand mol2 file" ,
                        type=str,
                        required=False)

    parser.add_argument("-s","--system",
                        help="sys either complex (com) or aqueous (aq)" ,
                        type=str,
                        required=False)

    parser.add_argument("-o","--output",
                        help="output basename of merged files, sc and ti masks" ,
                        type=str,
                        required=False)


    parser.add_argument("-cl","--calc",
                        help="type of calculation; rbfe or rsfe" ,
                        type=str,
                        required=False)


    args = parser.parse_args()


    if ".mol2" in args.pdb:
        q = parmed.load_file(args.pdb)
        q.save("{}.pdb".format(args.pdb.split('.mol2')[0]),overwrite=True)
        p = parmed.load_file("{}.pdb".format(args.pdb.split('.mol2')[0]))
    else:

        p = parmed.load_file(args.pdb)

    f1 = open ("{}".format(args.map),'r')
    lines = f1.readlines()

    if args.calc == "rbfe":
    	if args.orglig:
            p1 = CopyPDB (p)
            mol1 = Extract(p1, ":{}".format(args.orglig))
            mol2 = CopyPDB (mol1)

            p2 = CopyPDB (p)
            nl = Strip(p2, ":{}".format(args.orglig))
    else:
        mol1 = CopyPDB (p)
        mol2 = CopyPDB (p)


    mol2.residues[0].name = "{}".format(args.tarlig)

    mappedmol1atomnames=[]
    mappedmol2atomnames=[]
    for line in lines:
        atoms = re.split('=>|\n',line.replace(" ",""))
        mappedmol1atomnames.append(atoms[0])
        mappedmol2atomnames.append(atoms[1])

    mol1atomnames=[]
    for a in mol1.atoms:
        mol1atomnames.append(a.name)

    atomstodel=set (mol1atomnames) ^ set (mappedmol1atomnames)

    if len(atomstodel) > 0:
    	atoms=','.join([str(i) for i in atomstodel])
    	atoms="".join(('@',atoms))
    	mol2 = Strip(mol2,atoms)

    mol1tomol2map = []
    for j in mol2.atoms:
        count=0
        for i in mappedmol1atomnames:
            if ( j.name == i ):
                mol1tomol2map.append(count)
            count+=1

    count=0
    for j in mol2.atoms:
        j.name = mappedmol2atomnames[mol1tomol2map[count]]
        count+=1

    #parmed.tools.writeCoordinates(mol1,"mol1.pdb").execute()
    #parmed.tools.writeCoordinates(mol2,"mol2.pdb").execute()
    #parmed.tools.writeCoordinates(nl,"nonligand.pdb").execute()


    if args.system == "com":
        merge = Joincom(mol1,mol2,nl)
    elif args.system == "aq":
        merge = Joinaq(mol1,mol2)

    parmed.tools.writeCoordinates(merge,"merged_{}.pdb".format(args.output)).execute()

    fh = open("merged_%s.seq"%(args.output),"w")
    seq = GetResSeq( merge )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    fh.write(seqstr)
    fh.close()

    m1names = []
    m2names = []
    m1 = parmed.load_file(args.orgmol2)
    m2 = parmed.load_file(args.tarmol2)
    for j in m1.atoms:
        m1names.append(j.name)

    for j in m2.atoms:
        m2names.append(j.name)

    scatoms1 = set (m1names) ^ set (mappedmol1atomnames)
    scatoms2 = set (m2names) ^ set (mappedmol2atomnames)

    if len(scatoms1) > 0:
    	scmask1 = ','.join([str(i) for i in scatoms1])
    	scmask1 = "".join((":{}@".format(args.orglig),scmask1))
    else:
        scmask1 = '""'

    if len(scatoms2) > 0:
    	scmask2 = ','.join([str(i) for i in scatoms2])
    	scmask2 = "".join((":{}@".format(args.tarlig),scmask2))
    else:
        scmask2 = '""'

    timask1 = ":{}".format(args.orglig)
    timask2 = ":{}".format(args.tarlig)

    fh1 = open("{}.scmask1".format(args.output), "w")
    fh1.write(scmask1)
    fh1.close

    fh2 = open("{}.scmask2".format(args.output), "w")
    fh2.write(scmask2)
    fh2.close

    fh3 = open("{}.timask1".format(args.output), "w")
    fh3.write(timask1)
    fh3.close

    fh4 = open("{}.timask2".format(args.output), "w")
    fh4.write(timask2)
    fh4.close


EOF
        chmod a+x parsePDB.py
        ./parsePDB.py -p ${pdbfile} -m ${mapfile} -ol ${ligmask} -tl ${ligcopy} -om ${orgmol2} -tm ${tarmol2} -s ${s} -o ${out} -cl ${calc}
}

function parse_pdb_twostate {
        local pdbfile1=$1; shift
        local pdbfile2=$1; shift
        local mapfile=$1; shift
        local lig1mask=$1; shift
        local lig2mask=$1; shift
        local lig1mol2=$1; shift
        local lig2mol2=$1; shift
        local s=$1; shift
        local out=$1; shift
        cat << EOF > parsePDB.py
#!/usr/bin/env python3


def CopyPDB( pdb ):
    import copy
    p = copy.copy( pdb )
    p.coordinates = copy.copy( pdb.coordinates )
    return p


def Strip( pdb, mask ):
    p = CopyPDB( pdb )
    p.strip( "%s"%(mask) )
    return p


def Extract( pdb, mask ):
    return Strip( pdb, "!(%s)"%(mask) )

def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]

def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


def Joincom( mol1, mol2, nonligand ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    elif nonligand.coordinates.shape[0] == 0:
        pq=nonligand
    else:
        pq = mol1 + mol2 + nonligand
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates, nonligand.coordinates) )
        return pq

def Joinaq( mol1, mol2 ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    else:
        pq = mol1 + mol2
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates) )
        return pq


if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    from sys import exit

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Takes as input - a PDB file, an atommap file, two ligand masks \"orglig\" and \"tarlig\", name of the two ligand mol2 files as \"orgmol2\" and \"tarmol2\", and output basename. Splits PDB into water.pdb, nonwater.pdb, mol1.pdb containing \"orglig\", and mol2.pdb containing a new residue \"tarlig\" which consists of a subset of \"orglig\" atoms that are mapped as per provided atommap file. Merges the different pdbs into a single pdb and generates the scmasks and timasks" )

    parser.add_argument("-p1","--pdb1",
                        help="PDB file 1",
                        type=str,
                        required=True)

    parser.add_argument("-p2","--pdb2",
                        help="PDB file 2",
                        type=str,
                        required=True)

    parser.add_argument("-m","--map",
                        help="Atom mapfile",
                        type=str,
                        required=True )

    parser.add_argument("-l1","--lig1",
                        help="name of ligand 1 in PDB file 1" ,
                        type=str,
                        required=False)

    parser.add_argument("-l2","--lig2",
                        help="name of ligand 2 in PDB file 2" ,
                        type=str,
                        required=False)

    parser.add_argument("-m1","--lig1mol2",
                        help="name of mol2 file of ligand 1" ,
                        type=str,
                        required=False)

    parser.add_argument("-m2","--lig2mol2",
                        help="name of mol2 file of ligand 2" ,
                        type=str,
                        required=False)

    parser.add_argument("-s","--system",
                        help="sys either complex (com) or aqueous (aq)" ,
                        type=str,
                        required=False)

    parser.add_argument("-o","--output",
                        help="output basename of merged files, sc and ti masks" ,
                        type=str,
                        required=False)


    args = parser.parse_args()

    f1 = open ("{}".format(args.map),'r')
    lines = f1.readlines()

    mappedmol1atomnames=[]
    mappedmol2atomnames=[]
    for line in lines:
        atoms = re.split('=>|\n',line.replace(" ",""))
        mappedmol1atomnames.append(atoms[0])
        mappedmol2atomnames.append(atoms[1])

    if ".mol2" in args.pdb1:
        q = parmed.load_file(args.pdb1)
        q.save("{}.pdb".format(args.pdb1.split('.mol2')[0]),overwrite=True)
        pdb1 = parmed.load_file("{}.pdb".format(args.pdb1.split('.mol2')[0]))
    else:
        pdb1 = parmed.load_file(args.pdb1)

    if ".mol2" in args.pdb2:
        q = parmed.load_file(args.pdb2)
        q.save("{}.pdb".format(args.pdb2.split('.mol2')[0]),overwrite=True)
        pdb2 = parmed.load_file("{}.pdb".format(args.pdb2.split('.mol2')[0]))
    else:
        pdb2 = parmed.load_file(args.pdb2)


#### manipulate pdb1 
    if args.system == "com":
        pdb1copy1 = CopyPDB (pdb1)
        p1mol1 = Extract(pdb1copy1, ":{}".format(args.lig1))
        p1mol2 = CopyPDB (p1mol1)
        pdb1copy2 = CopyPDB (pdb1)
        p1nl = Strip(pdb1copy2, ":{}".format(args.lig1))
    else:
        p1mol1 = CopyPDB (pdb1)
        p1mol2 = CopyPDB (pdb1)

    p1mol2.residues[0].name = "{}".format(args.lig2)

    p1mol1atomnames=[]
    for a in p1mol1.atoms:
        p1mol1atomnames.append(a.name)

    atomstodel=set (p1mol1atomnames) ^ set (mappedmol1atomnames)

    if len(atomstodel) > 0:
    	atoms=','.join([str(i) for i in atomstodel])
    	atoms="".join(('@',atoms))
    	p1mol2 = Strip(p1mol2,atoms)

    p1mol1tomol2map = []
    for j in p1mol2.atoms:
        count=0
        for i in mappedmol1atomnames:
            if ( j.name == i ):
                p1mol1tomol2map.append(count)
            count+=1

    count=0
    for j in p1mol2.atoms:
        j.name = mappedmol2atomnames[p1mol1tomol2map[count]]
        count+=1
###################

#### manipulate pdb2
    if args.system == "com":
        pdb2copy1 = CopyPDB (pdb2)
        p2mol1 = Extract(pdb2copy1, ":{}".format(args.lig2))
        p2mol2 = CopyPDB (p2mol1)
        pdb2copy2 = CopyPDB (pdb2)
        p2nl = Strip(pdb2copy2, ":{}".format(args.lig2))
    else:
        p2mol1 = CopyPDB (pdb2)
        p2mol2 = CopyPDB (pdb2)

    p2mol2.residues[0].name = "{}".format(args.lig1)

    p2mol1atomnames=[]
    for a in p2mol1.atoms:
        p2mol1atomnames.append(a.name)

    atomstodel=set (p2mol1atomnames) ^ set (mappedmol2atomnames)

    if len(atomstodel) > 0:
    	atoms=','.join([str(i) for i in atomstodel])
    	atoms="".join(('@',atoms))
    	p2mol2 = Strip(p2mol2,atoms)

    p2mol1tomol2map = []
    for j in p2mol2.atoms:
        count=0
        for i in mappedmol2atomnames:
            if ( j.name == i ):
                p2mol1tomol2map.append(count)
            count+=1

    count=0
    for j in p2mol2.atoms:
        j.name = mappedmol1atomnames[p2mol1tomol2map[count]]
        count+=1
###################


    #parmed.tools.writeCoordinates(mol1,"mol1.pdb").execute()
    #parmed.tools.writeCoordinates(mol2,"mol2.pdb").execute()
    #parmed.tools.writeCoordinates(nl,"nonligand.pdb").execute()


    if args.system == "com":
        merge1 = Joincom(p1mol1,p1mol2,p1nl)
        merge2 = Joincom(p2mol2,p2mol1,p2nl)
    elif args.system == "aq":
        merge1 = Joinaq(p1mol1,p1mol2)
        merge2 = Joinaq(p2mol2,p2mol1)

    parmed.tools.writeCoordinates(merge1,"merged1_{}.pdb".format(args.output)).execute()
    parmed.tools.writeCoordinates(merge2,"merged2_{}.pdb".format(args.output)).execute()

    p1fh = open("merged1_%s.seq"%(args.output),"w")
    seq = GetResSeq( merge1 )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    p1fh.write(seqstr)
    p1fh.close()

    p2fh = open("merged2_%s.seq"%(args.output),"w")
    seq = GetResSeq( merge2 )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    p2fh.write(seqstr)
    p2fh.close()


    m1names = []
    m2names = []
    m1 = parmed.load_file(args.lig1mol2)
    m2 = parmed.load_file(args.lig2mol2)
    for j in m1.atoms:
        m1names.append(j.name)

    for j in m2.atoms:
        m2names.append(j.name)

    scatoms1 = set (m1names) ^ set (mappedmol1atomnames)
    scatoms2 = set (m2names) ^ set (mappedmol2atomnames)

    if len(scatoms1) > 0:
        scmask1 = ','.join([str(i) for i in scatoms1])
        scmask1 = "".join((":{}@".format(args.lig1),scmask1))
    else:
        scmask1 = '""'

    if len(scatoms2) > 0:
        scmask2 = ','.join([str(i) for i in scatoms2])
        scmask2 = "".join((":{}@".format(args.lig2),scmask2))
    else:
        scmask2 = '""'

    timask1 = ":{}".format(args.lig1)
    timask2 = ":{}".format(args.lig2)

    fh1 = open("{}.scmask1".format(args.output), "w")
    fh1.write(scmask1)
    fh1.close

    fh2 = open("{}.scmask2".format(args.output), "w")
    fh2.write(scmask2)
    fh2.close

    fh3 = open("{}.timask1".format(args.output), "w")
    fh3.write(timask1)
    fh3.close

    fh4 = open("{}.timask2".format(args.output), "w")
    fh4.write(timask2)
    fh4.close


EOF
        chmod a+x parsePDB.py
        ./parsePDB.py -p1 ${pdbfile1} -p2 ${pdbfile2} -m ${mapfile} -l1 ${lig1mask} -l2 ${lig2mask} -m1 ${lig1mol2} -m2 ${lig2mol2} -s ${s} -o ${out}
}


###################################################################
###################################################################
###################################################################
###################################################################
####################### MAIN FUNCTION #############################
###################################################################
###################################################################
###################################################################


function create_box_rbfe {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
	local boxbuild=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local ionconc=$1; shift
        local mapfile=$1; shift
        local s=$1; shift
	local ticalc=$1; shift
	local bidirection=$1; shift
	local translist=("$@")

        local mollist=(); local liglist=(); local numnonstdlist=()
        while read line; do
                IFS=' ' read -ra args <<< $line
                mollist+=(${args[0]}); liglist+=(${args[1]}); numnonstdlist+=(${args[2]})
        done < ${mapfile}

	# unnecessary now but may be needed later. 
        if [ "${s}" == "com" ]; then load=pdbseq; else load=pdbseq; fi

	###
	mergedpdbs=(); lig1s=(); lig2s=(); nnstds=(); mergedpdbsrev=(); nnstdsrev=()
	for i in "${!translist[@]}";do
	       stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
	       indA=$(get_index "${stA}" "${mollist[@]}"); indB=$(get_index "${stB}" "${mollist[@]}")
	       ligA=${liglist[${indA}]}; ligB=${liglist[${indB}]}
	       nnstdA=${numnonstdlist[${indA}]}; nnstdB=${numnonstdlist[${indB}]}
	       
	       parse_pdb "${stA}.pdb" "${stA}~${stB}.map.txt" "${ligA}" "${ligB}" "${stA}_0.mol2" "${stB}_0.mol2" "$s" "${stA}~${stB}" "${ticalc}"
	       mergedpdbs+=(merged_${stA}~${stB})
	       lig1s+=(${stA}); lig2s+=(${stB})
	       nnstds+=(${nnstdA})

	       if [ "${bidirection}" == "true" ]; then
		       parse_pdb "${stB}.pdb" "${stB}~${stA}.map.txt" "${ligB}" "${ligA}" "${stB}_0.mol2" "${stA}_0.mol2" "$s" "${stB}~${stA}" "${ticalc}"
		       mergedpdbsrev+=(merged_${stB}~${stA})
		       nnstdsrev+=(${nnstdB})
	       fi
	done
	###
        # write and run tleap to generate initial parm file
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[0]}" "${lig1s[0]}" "${nnstds[0]}" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "0" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log

        # calculate water and number of ions necessary to reach desired ion conc
        nwat=$(calcwaterinparm merged.parm7)
	nions=$(calcMDions "$nwat" "${ionconc}")
        #echo "No. of water and ions required in ${list[0]} : $nwat $nions"

        # generate parm file with calculated number of ions
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[0]}" "${lig1s[0]}" "${nnstds[0]}" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log
        mv merged.parm7 "${lig1s[0]}~${lig2s[0]}_${s}".parm7; mv merged.rst7 "${lig1s[0]}~${lig2s[0]}_${s}".rst7

        # calculate water and ions in final parm file for first system. Remaining systems will be built with these number of water and ions.
        nwat=$(calcwaterinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nions=$(calcionsinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nsodium=$(calcsodiuminparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nchloride=$(calcchlorideinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        #echo "final number of water and ions in ${lig1s[0]}_${s} : $nwat $nions $nsodium $nchloride"

	if [ "${#mergedpdbs[@]}" -gt 1 ]; then
        #	mergedpdbs=(${mergedpdbs[@]:1}); lig1s=(${lig1s[@]:1}); lig2s=(${lig2s[@]:1}); nnstds=(${nnstds[@]:1})

		if [ "${boxbuild}" != 2 ]; then

			# build all parm files with nions. each will have different number of water molecules.
			for m in "${!mergedpdbs[@]}";do
				if [ "${s}" == "com" ]; then load=pdbseq; else load=pdbseq; fi

        			write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        			tleap -s -f tleap.in > output
        			mv merged.parm7 "${lig1s[$m]}~${lig2s[$m]}_${s}".parm7; mv merged.rst7 "${lig1s[$m]}~${lig2s[$m]}_${s}".rst7

				if [ "${bidirection}" == "true" ]; then
					write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbsrev[$m]}" "${lig2s[$m]}" "${nnstdsrev[$m]}" "${lig1s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        				tleap -s -f tleap.in > output
        				mv merged.parm7 "${lig2s[$m]}~${lig1s[$m]}_${s}".parm7; mv merged.rst7 "${lig2s[$m]}~${lig1s[$m]}_${s}".rst7
				fi

			done

		else
			# build all parm files with nions and nwaters
			for m in "${!mergedpdbs[@]}";do
				fix_solvent "${pff}" "${lff}" "${wm}" "${mergedpdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
				mv out.parm7 ${lig1s[$m]}~${lig2s[$m]}_${s}.parm7; mv out.rst7 ${lig1s[$m]}~${lig2s[$m]}_${s}.rst7

				if [ "${bidirection}" == "true" ]; then
					fix_solvent "${pff}" "${lff}" "${wm}" "${mergedpdbsrev[$m]}" "${lig2s[$m]}" "${nnstdsrev[$m]}" "${lig1s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
					mv out.parm7 ${lig2s[$m]}~${lig1s[$m]}_${s}.parm7; mv out.rst7 ${lig2s[$m]}~${lig1s[$m]}_${s}.rst7
				fi
			done

		fi
	fi

	##########################################
        # setup H-mass repartitioning
        if [ "${hmr}" == "true" ]; then
		for m in "${!translist[@]}";do
			#mol="${translist[$m]}_${s}"
			mol="${lig1s[$m]}~${lig2s[$m]}_${s}"

                	if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                	cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                	parmed -i hmr.in -p ${mol}.parm7 -c ${mol}.rst7 >> output 2>&1
                	sleep 1
                	mv hmr.parm7 ${mol}.parm7
                	mv hmr.rst7  ${mol}.rst7
                	rm -rf hmr.in

			if [ "${bidirection}" == "true" ]; then
				#mol="${translistrev[$m]}_${s}"
				mol="${lig2s[$m]}~${lig1s[$m]}_${s}"
                        	cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                        	parmed -i hmr.in -p ${mol}.parm7 -c ${mol}.rst7 >> output 2>&1
                        	sleep 1
                        	mv hmr.parm7 ${mol}.parm7
                        	mv hmr.rst7  ${mol}.rst7
                        	rm -rf hmr.in
			fi
		done
    	fi
        ##########################################


	##########################################
        # double check prepared systems
        for m in "${!translist[@]}";do
                #mol="${translist[$m]}_${s}"
		mol="${lig1s[$m]}~${lig2s[$m]}_${s}"
                nwat=$(calcwaterinparm ${mol}.parm7); nsod=$(calcsodiuminparm ${mol}.parm7); ncl=$(calcchlorideinparm ${mol}.parm7)
                printf "\n${mol} has ${nwat} waters, ${nsod} Na+ ions, and ${ncl} Cl- ions\n"

		if [ "${bidirection}" == "true" ]; then
			#mol="${translistrev[$m]}_${s}"
			mol="${lig2s[$m]}~${lig1s[$m]}_${s}"
			nwat=$(calcwaterinparm ${mol}.parm7); nsod=$(calcsodiuminparm ${mol}.parm7); ncl=$(calcchlorideinparm ${mol}.parm7)
			printf "\n${mol} has ${nwat} waters, ${nsod} Na+ ions, and ${ncl} Cl- ions\n"

		fi
        done
        ##########################################



}


function create_box_twostate {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local boxbuild=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local ionconc=$1; shift
        local mapfile=$1; shift
        local s=$1; shift
        local translist=("$@")

	if [ "${s}" == "com" ]; then

        	local mollist=(); local liglist=(); local numnonstdlist=()
        	while read line; do
                	IFS=' ' read -ra args <<< $line
                	mollist+=(${args[0]}); liglist+=(${args[1]}); numnonstdlist+=(${args[2]})
        	done < ${mapfile}
		load=pdbseq

	else

                local mollist=(); local liglist=()
                while read line; do
                        IFS=' ' read -ra args <<< $line
                        mollist+=(${args[0]}); liglist+=(${args[1]})
                done < ${mapfile}
                load=pdbseq
	fi


	merged1pdbs=(); merged2pdbs=(); lig1s=(); lig2s=(); nnstds=()
	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		indA=$(get_index "${stA}" "${mollist[@]}"); indB=$(get_index "${stB}" "${mollist[@]}")
		ligA=${liglist[${indA}]}; ligB=${liglist[${indB}]}
		if [ "${s}" == "com" ]; then
			parse_pdb_twostate "${stA}.pdb" "${stB}.pdb" "${stA}~${stB}.map.txt" "${ligA}" "${ligB}" "${stA}_0.mol2" "${stB}_0.mol2" "$s" "${stA}~${stB}"
			nnstdA=${numnonstdlist[${indA}]}; nnstdB=${numnonstdlist[${indB}]}
			nnstds+=(${nnstdA})
		else
			parse_pdb_twostate "${stA}_0.mol2" "${stB}_0.mol2" "${stA}~${stB}.map.txt" "${ligA}" "${ligB}" "${stA}_0.mol2" "${stB}_0.mol2" "$s" "${stA}~${stB}"
			nnstds+=("1")
		fi

		merged1pdbs+=(merged1_${stA}~${stB}); merged2pdbs+=(merged2_${stA}~${stB})
		lig1s+=(${stA}); lig2s+=(${stB})
	done


        ###
        # write and run tleap to generate initial parm file
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${merged1pdbs[0]}" "${lig1s[0]}" "${nnstds[0]}" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "0" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log

        # calculate water and number of ions necessary to reach desired ion conc
        nwat=$(calcwaterinparm merged.parm7)
        nions=$(calcMDions "$nwat" "${ionconc}")
        #echo "No. of water and ions required in ${merged1pbds[0]} : $nwat $nions"

        # generate parm file with calculated number of ions
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${merged1pdbs[0]}" "${lig1s[0]}" "${nnstds[0]}" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log
        mv merged.parm7 "${lig1s[0]}~${lig2s[0]}-1_${s}".parm7; mv merged.rst7 "${lig1s[0]}~${lig2s[0]}-1_${s}".rst7

        # calculate water and ions in final parm file for first system. Remaining systems will be built with these number of water and ions.
        nwat=$(calcwaterinparm ${lig1s[0]}~${lig2s[0]}-1_${s}.parm7)
        nsodium=$(calcsodiuminparm ${lig1s[0]}~${lig2s[0]}-1_${s}.parm7)
        nchloride=$(calcchlorideinparm ${lig1s[0]}~${lig2s[0]}-1_${s}.parm7)
        #echo "final number of water, sodium, chloride ions in ${lig1s[0]}~${lig2s[0]}-1_${s} : $nwat $nsodium $nchloride"

	fix_solvent "${pff}" "${lff}" "${wm}" "${merged2pdbs[0]}" "${lig1s[0]}" "${nnstds[0]}" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
	mv out.parm7 ${lig1s[0]}~${lig2s[0]}-2_${s}.parm7; mv out.rst7 ${lig1s[0]}~${lig2s[0]}-2_${s}.rst7


        #nwat=$(calcwaterinparm ${lig1s[0]}~${lig2s[0]}-2_${s}.parm7)
        #nsodium=$(calcsodiuminparm ${lig1s[0]}~${lig2s[0]}-2_${s}.parm7)
        #nchloride=$(calcchlorideinparm ${lig1s[0]}~${lig2s[0]}-2_${s}.parm7)
        #echo "final number of water, sodium, chloride ions in ${lig1s[0]}~${lig2s[0]}-2_${s} : $nwat $nsodium $nchloride"

        if [ "${#merged1pdbs[@]}" -gt 1 ]; then
                merged1pdbs=(${merged1pdbs[@]:1}); merged2pdbs=(${merged2pdbs[@]:1}); lig1s=(${lig1s[@]:1}); lig2s=(${lig2s[@]:1}); nnstds=(${nnstds[@]:1})

                if [ "${boxbuild}" != 2 ]; then

                        # build all parm files with nions. each will have different number of water molecules.
                        for m in "${!merged1pdbs[@]}";do
                                if [ "${s}" == "com" ]; then load=pdbseq; else load=pdbseq; fi

                                write_tleap_merged "${pff}" "${lff}" "${wm}" "${merged1pdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"

                                tleap -s -f tleap.in > output
                                mv merged.parm7 "${lig1s[$m]}~${lig2s[$m]}-1_${s}".parm7; mv merged.rst7 "${lig1s[$m]}~${lig2s[$m]}-1_${s}".rst7
        			nwat=$(calcwaterinparm ${lig1s[$m]}~${lig2s[$m]}-1_${s}.parm7)
        			nsodium=$(calcsodiuminparm ${lig1s[$m]}~${lig2s[$m]}-1_${s}.parm7)
        			nchloride=$(calcchlorideinparm ${lig1s[$m]}~${lig2s[$m]}-1_${s}.parm7)
				fix_solvent "${pff}" "${lff}" "${wm}" "${merged2pdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
				mv out.parm7 ${lig1s[$m]}~${lig2s[$m]}-2_${s}.parm7; mv out.rst7 ${lig1s[$m]}~${lig2s[$m]}-2_${s}.rst7

                        done

                else

                        # build all parm files with nions and nwaters
                        for m in "${!merged1pdbs[@]}";do
                                fix_solvent "${pff}" "${lff}" "${wm}" "${merged1pdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
                                mv out.parm7 ${lig1s[$m]}~${lig2s[$m]}-1_${s}.parm7; mv out.rst7 ${lig1s[$m]}~${lig2s[$m]}-1_${s}.rst7
                                fix_solvent "${pff}" "${lff}" "${wm}" "${merged2pdbs[$m]}" "${lig1s[$m]}" "${nnstds[$m]}" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
                                mv out.parm7 ${lig1s[$m]}~${lig2s[$m]}-2_${s}.parm7; mv out.rst7 ${lig1s[$m]}~${lig2s[$m]}-2_${s}.rst7
                        done

                fi
        fi

        ##########################################
        # setup H-mass repartitioning
        for m in "${!translist[@]}";do
                mol1="${translist[$m]}-1_${s}"; mol2="${translist[$m]}-2_${s}"
                if [ "${hmr}" == "true" ]; then
                        if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                        cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                        parmed -i hmr.in -p ${mol1}.parm7 -c ${mol1}.rst7 >> output 2>&1
                        mv hmr.parm7 ${mol1}.parm7; mv hmr.rst7  ${mol1}.rst7
                        parmed -i hmr.in -p ${mol2}.parm7 -c ${mol2}.rst7 >> output 2>&1
                        mv hmr.parm7 ${mol2}.parm7; mv hmr.rst7  ${mol2}.rst7
                        rm -rf hmr.in
                fi
        done
        ##########################################


        ##########################################
        # double check prepared systems
        for m in "${!translist[@]}";do
                mol1="${translist[$m]}-1_${s}"; mol2="${translist[$m]}-2_${s}"
		nwat1=$(calcwaterinparm ${mol1}.parm7); nsod1=$(calcsodiuminparm ${mol1}.parm7); ncl1=$(calcchlorideinparm ${mol1}.parm7)
		nwat2=$(calcwaterinparm ${mol2}.parm7); nsod2=$(calcsodiuminparm ${mol2}.parm7); ncl2=$(calcchlorideinparm ${mol2}.parm7)
                printf "\n${mol1} has ${nwat1} waters, ${nsod1} Na+ ions, and ${ncl1} Cl- ions\n"
                printf "\n${mol2} has ${nwat2} waters, ${nsod2} Na+ ions, and ${ncl2} Cl- ions\n"
        done
        ##########################################

}


function create_box_rsfe {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local boxbuild=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local ionconc=$1; shift
        local mapfile=$1; shift
        local s=$1; shift
        local ticalc=$1; shift
	local bidirection=$1; shift
        local translist=("$@")

	local mollist=(); local liglist=()
        while read line; do
                IFS=' ' read -ra args <<< $line
                mollist+=(${args[0]}); liglist+=(${args[1]})
        done < ${mapfile}

        # unnecessary now but may be needed later.
        if [ "${s}" == "com" ]; then load=pdbseq; else load=pdbseq; fi

        ###
        mergedpdbs=(); lig1s=(); lig2s=(); mergedpdbsrev=()
        for i in "${!translist[@]}";do
               stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
               indA=$(get_index "${stA}" "${mollist[@]}"); indB=$(get_index "${stB}" "${mollist[@]}")
               ligA=${liglist[${indA}]}; ligB=${liglist[${indB}]}

               parse_pdb "${stA}_0.mol2" "${stA}~${stB}.map.txt" "${ligA}" "${ligB}" "${stA}_0.mol2" "${stB}_0.mol2" "$s" "${stA}~${stB}" "${ticalc}"
               mergedpdbs+=(merged_${stA}~${stB})
               lig1s+=(${stA}); lig2s+=(${stB})

	       if [ "${bidirection}" == "true" ]; then
		       parse_pdb "${stB}_0.mol2" "${stB}~${stA}.map.txt" "${ligB}" "${ligA}" "${stB}_0.mol2" "${stA}_0.mol2" "$s" "${stB}~${stA}" "${ticalc}"
		       mergedpdbsrev+=(merged_${stB}~${stA})
	       fi
        done
        ###
        # write and run tleap to generate initial parm file
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[0]}" "${lig1s[0]}" "1" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "0" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log

        # calculate water and number of ions necessary to reach desired ion conc
        nwat=$(calcwaterinparm merged.parm7)
        nions=$(calcMDions "$nwat" "${ionconc}")
        #echo "No. of water and ions required in ${list[0]} : $nwat $nions"

        # generate parm file with calculated number of ions
        write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[0]}" "${lig1s[0]}" "1" "${lig2s[0]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log
        mv merged.parm7 "${lig1s[0]}~${lig2s[0]}_${s}".parm7; mv merged.rst7 "${lig1s[0]}~${lig2s[0]}_${s}".rst7

        # calculate water and ions in final parm file for first system. Remaining systems will be built with these number of water and ions.
        nwat=$(calcwaterinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nions=$(calcionsinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nsodium=$(calcsodiuminparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)
        nchloride=$(calcchlorideinparm ${lig1s[0]}~${lig2s[0]}_${s}.parm7)

        #echo "final number of water and ions in ${translist[0]}_${s} : $nwat $nions"

        if [ "${#mergedpdbs[@]}" -gt 1 ]; then
                #mergedpdbs=(${mergedpdbs[@]:1}); lig1s=(${lig1s[@]:1}); lig2s=(${lig2s[@]:1}); nnstds=(${nnstds[@]:1})

                if [ "${boxbuild}" != 2 ]; then

                        # build all parm files with nions. each will have different number of water molecules.
                        for m in "${!mergedpdbs[@]}";do
                                write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbs[$m]}" "${lig1s[$m]}" 1 "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
                                tleap -s -f tleap.in > output
                                mv merged.parm7 "${lig1s[$m]}~${lig2s[$m]}_${s}".parm7; mv merged.rst7 "${lig1s[$m]}~${lig2s[$m]}_${s}".rst7

				
				if [ "${bidirection}" == "true" ]; then
					write_tleap_merged "${pff}" "${lff}" "${wm}" "${mergedpdbsrev[$m]}" "${lig2s[$m]}" 1 "${lig1s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
					tleap -s -f tleap.in > output
					mv merged.parm7 "${lig2s[$m]}~${lig1s[$m]}_${s}".parm7; mv merged.rst7 "${lig2s[$m]}~${lig1s[$m]}_${s}".rst7
				fi
                        done

                else

                        # build all parm files with nions and nwaters
                        for m in "${!mergedpdbs[@]}";do
                                fix_solvent "${pff}" "${lff}" "${wm}" "${mergedpdbs[$m]}" "${lig1s[$m]}" "1" "${lig2s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
                                mv out.parm7 ${lig1s[$m]}~${lig2s[$m]}_${s}.parm7; mv out.rst7 ${lig1s[$m]}~${lig2s[$m]}_${s}.rst7

				if [ "${bidirection}" == "true" ]; then
					fix_solvent "${pff}" "${lff}" "${wm}" "${mergedpdbsrev[$m]}" "${lig2s[$m]}" "1" "${lig1s[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
					mv out.parm7 ${lig2s[$m]}~${lig1s[$m]}_${s}.parm7; mv out.rst7 ${lig2s[$m]}~${lig1s[$m]}_${s}.rst7
				fi
                        done

                fi
        fi

        ##########################################
        # setup H-mass repartitioning
	if [ "${hmr}" == "true" ]; then
        	for m in "${!translist[@]}";do
                	#mol="${translist[$m]}_${s}"
			mol="${lig1s[$m]}~${lig2s[$m]}_${s}"
                        if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                        cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                        parmed -i hmr.in -p ${mol}.parm7 -c ${mol}.rst7 >> output 2>&1
                        sleep 1
                        mv hmr.parm7 ${mol}.parm7
                        mv hmr.rst7  ${mol}.rst7
                        rm -rf hmr.in

			if [ "${bidirection}" == "true" ]; then
				#mol="${translistrev[$m]}_${s}"
				mol="${lig2s[$m]}~${lig1s[$m]}_${s}"
				cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
				parmed -i hmr.in -p ${mol}.parm7 -c ${mol}.rst7 >> output 2>&1
				sleep 1
				mv hmr.parm7 ${mol}.parm7
				mv hmr.rst7  ${mol}.rst7
				rm -rf hmr.in
			fi
		done
	fi
        ##########################################


        ##########################################
        # double check prepared systems
        for m in "${!translist[@]}";do
                #mol="${translist[$m]}_${s}"
		mol="${lig1s[$m]}~${lig2s[$m]}_${s}"
                nwat=$(calcwaterinparm ${mol}.parm7)
                nsod=$(calcsodiuminparm ${mol}.parm7)
		ncl=$(calcchlorideinparm ${mol}.parm7)

                printf "\n${mol} has ${nwat} waters, ${nsod} Na+ ions, and ${ncl} Cl- ions"

		if [ "${bidirection}" == "true" ]; then
			#mol="${translistrev[$m]}_${s}"
			mol="${lig2s[$m]}~${lig1s[$m]}_${s}"
			nwat=$(calcwaterinparm ${mol}.parm7); nsod=$(calcsodiuminparm ${mol}.parm7); ncl=$(calcchlorideinparm ${mol}.parm7)
			printf "\n${mol} has ${nwat} waters, ${nsod} Na+ ions, and ${ncl} Cl- ions\n"
		fi
        done
        ##########################################

}


function create_box_asfe {

        local pff=$1; shift
        local lff=$1; shift
        local wm=$1; shift
        local boxbuild=$1; shift
        local mdboxshape=$1; shift
        local rbuf=$1; shift
        local ionconc=$1; shift
        local mapfile=$1; shift
        local s=$1; shift
        local ticalc=$1; shift
        local translist=("$@")

        local mollist=(); local liglist=()
        while read line; do
                IFS=' ' read -ra args <<< $line
                mollist+=(${args[0]}); liglist+=(${args[1]})
        done < ${mapfile}

        # unnecessary now but may be needed later.
        if [ "${s}" == "com" ]; then load=loadmol2; else load=loadmol2; fi

        ###
        ligs=()
        for i in "${!translist[@]}";do
		ind=$(get_index "${translist[$i]}" "${mollist[@]}")
		lig=${liglist[${ind}]}
		sed -i "s/LIG/${lig}/g" ${mollist[$i]}_0.mol2 ${mollist[$i]}_0.lib
               	ligs+=(${lig})
        done
        ###
        # write and run tleap to generate initial parm file
        write_tleap_asfe "${pff}" "${lff}" "${wm}" "${mollist[0]}" "${mdboxshape}" "${rbuf}" "${load}" "0" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log

        # calculate water and number of ions necessary to reach desired ion conc
        nwat=$(calcwaterinparm merged.parm7)
        nions=$(calcMDions "$nwat" "${ionconc}")
        #echo "No. of water and ions required in ${list[0]} : $nwat $nions"

        # generate parm file with calculated number of ions
        write_tleap_asfe "${pff}" "${lff}" "${wm}" "${mollist[0]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
        tleap -s -f tleap.in > output; rm -rf leap.log
        mv merged.parm7 "${mollist[0]}_${s}".parm7; mv merged.rst7 "${mollist[0]}_${s}".rst7

        # calculate water and ions in final parm file for first system. Remaining systems will be built with these number of water and ions.
        nwat=$(calcwaterinparm ${mollist[0]}_${s}.parm7)
        nions=$(calcionsinparm ${mollist[0]}_${s}.parm7)
	nsodium=$(calcsodiuminparm ${mollist[0]}_${s}.parm7)
        nchloride=$(calcchlorideinparm ${mollist[0]}_${s}.parm7)

        #echo "final number of water and ions in ${translist[0]}_${s} : $nwat $nions"

        if [ "${#mollist[@]}" -gt 1 ]; then
                mollist=(${mollist[@]:1}); ligs=(${ligs[@]:1})

                if [ "${boxbuild}" != 2 ]; then

                        # build all parm files with nions. each will have different number of water molecules.
                        for m in "${!mollist[@]}";do
                                write_tleap_asfe "${pff}" "${lff}" "${wm}" "${mollist[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nions}" "${boxbuild}" "${s}"
                                tleap -s -f tleap.in > output
                                mv merged.parm7 "${mollist[$m]}_${s}".parm7; mv merged.rst7 "${mollist[$m]}_${s}".rst7
                        done

                else

                        # build all parm files with nions and nwaters
                        for m in "${!mollist[@]}";do
                                fix_solvent_asfe "${pff}" "${lff}" "${wm}" "${mollist[$m]}" "${mdboxshape}" "${rbuf}" "${load}" "${nsodium}" "${nchloride}" "${nwat}" "${boxbuild}" "${s}"
                                mv out.parm7 ${mollist[$m]}_${s}.parm7; mv out.rst7 ${mollist[$m]}_${s}.rst7
                        done

                fi
        fi

        ##########################################
        # setup H-mass repartitioning
        for m in "${!translist[@]}";do
                mol="${translist[$m]}_${s}"
                if [ "${hmr}" == "true" ]; then
                        if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                        cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                        parmed -i hmr.in -p ${mol}.parm7 -c ${mol}.rst7 >> output 2>&1
                        sleep 1
                        mv hmr.parm7 ${mol}.parm7
                        mv hmr.rst7  ${mol}.rst7
                        rm -rf hmr.in
                fi
        done
        ##########################################


        ##########################################
        # double check prepared systems
        for m in "${!translist[@]}";do
                mol="${translist[$m]}_${s}"
                nwat=$(calcwaterinparm ${mol}.parm7)
                nions=$(calcionsinparm ${mol}.parm7)
                printf "${mol} has ${nwat} waters, ${nions} Na+ ions, and ${nions} Cl- ions \n"
        done
        ##########################################

}


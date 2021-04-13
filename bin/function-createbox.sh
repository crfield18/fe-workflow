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

function GetNumRes
{
    local parm="$1"
    local cols=($(grep -A 2 'SOLVENT_POINTERS' ${parm} | tail -n 1))
    echo ${cols[1]}
}

function RunTLEAP
{
    local PDB="$1"
    local SEQ="$2"
    local TLIGAND="$3"
    local NIONS="$4"
    local rbuf="$5"
    local firstdel="$6"
    local lastdel="$7"

    cat <<EOF > tleap.in

source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip4pew
loadamberparams frcmod.ff14SB
loadamberparams frcmod.tip4pew
loadAmberParams frcmod.ionsjc_tip4pew
loadoff tip4pewbox.off

loadamberparams ${TLIGAND}.frcmod
loadoff ${TLIGAND}.lib

m = loadPdbUsingSeq ${PDB} { $(cat ${SEQ}) }
solvateBox m TIP4PEWBOX ${rbuf}
addions m Na+ 0
addions m Cl- 0

addions m Na+ ${NIONS}
addions m Cl- ${NIONS}

savepdb m tmp.pdb
x = loadpdb tmp.pdb
setbox x centers
EOF

    if [ ! -z "${firstdel}" -a ! -z "${lastdel}" ]; then
        for res in $(seq ${lastdel} -1 ${firstdel}); do
            echo "remove x x.${res}" >> tleap.in
        done
    fi

    cat <<EOF >> tleap.in
saveamberparm x out.parm7 out.rst7
quit
EOF
    tleap -s -f tleap.in >> log
}

function fix_solvent {

        local args=$*; args=($args)
        local varlist=(nwat nions base rbuf)
        local i=0
        for a in "${varlist[@]}"; do
                declare -n arr="$a"
                arr=${args[$i]}
                i=$(($i+1))
        done

        for iter in $(seq 100); do
                RunTLEAP "${base}.pdb" "${base}.seq" "${base}" "${nions}" "${rbuf}"
                current_nwat=$(calcwaterinparm "out.parm7")
                excess_waters=$(bc -l <<< "${current_nwat} - ${nwat}")
                #echo "Current waters - ${current_nwat} Excess waters - ${excess_waters} rbuf - ${rbuf}"
                if [ "${excess_waters}" -ge 0 ]; then
                        break
                else
                        rbuf=$(bc -l <<< "${rbuf} + 1.")
                fi
        done
        if [ ${excess_waters} -gt 0 ]; then
                lastdel=$(GetNumRes "out.parm7")
                firstdel=$(( ${lastdel} - ${excess_waters} + 1 ))
                #echo "Total residues - ${lastdel} Excess waters - ${excess_waters} Delete residues from - ${firstdel}"
                RunTLEAP "${base}.pdb" "${base}.seq" "${base}" "${nions}" "${rbuf}" "${firstdel}" "${lastdel}"
        fi

        current_nwat=$(calcwaterinparm "out.parm7")
        #echo "There are ${current_nwat} waters in out.parm7"
}


function create_box_equal {
       	local s=$1
	shift
       	local rbuf=$1
	shift
       	local ionconc=$1
	shift
	local list=("$@")
	
	out=()
	for i in "${!list[@]}";do
                out+=($(echo ${list[$i]}|awk -F "_" '{print $1}'))
        done

	cat <<EOF > get_req_water_in_MDbox.py
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
	chmod a+x get_req_water_in_MDbox.py

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
complex = loadPdbUsingSeq ${list[0]}.pdb { $(cat ${list[0]}.seq) }

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

savepdb complex ${list[0]}-solv.pdb
saveamberparm complex ${list[0]}-solv.parm7 ${list[0]}-solv.rst7

quit
EOF
	tleap -s -f tleap.in > output

	nwat=$(calcwaterinparm ${list[0]}-solv.parm7)
	nions=$(python2.7 ./get_req_water_in_MDbox.py $nwat ${ionconc})

	#echo "No. of water and ions required in ${list[0]}_${s} : $nwat $nions"

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
complex = loadPdbUsingSeq ${list[0]}.pdb { $(cat ${list[0]}.seq) }

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

addions complex Na+ ${nions}
addions complex Cl- ${nions}

#savepdb complex ${list[0]}-solv.pdb
saveamberparm complex ${list[0]}-solv.parm7 ${list[0]}-solv.rst7

quit
EOF
	tleap -s -f tleap.in > output


	
	mv ${list[0]}-solv.parm7 ${out[0]}_${s}.parm7; mv ${list[0]}-solv.rst7  ${out[0]}_${s}.rst7

	nwat=$(calcwaterinparm ${out[0]}_${s}.parm7)
	nions=$(calcionsinparm ${out[0]}_${s}.parm7)
	#echo "final number of water and ions in ${list[0]}_${s} : $nwat $nions"

	listnew=(${list[@]:1}); outnew=(${out[@]:1})

	for i in "${!listnew[@]}";do
		mol=${listnew[$i]}; outfile=${outnew[$i]}
		fix_solvent $nwat $nions ${mol} ${rbuf}
		mv out.parm7 ${outfile}_${s}.parm7; mv out.rst7 ${outfile}_${s}.rst7
	done

	for i in "${!list[@]}";do
		mol="${out[$i]}_${s}"
        	nwat=$(calcwaterinparm ${mol}.parm7)
        	nions=$(calcionsinparm ${mol}.parm7)
        	echo "${mol} has ${nwat} waters and ${nions} Na+ ions"
	done
	rm -rf get_req_water_in_MDbox.py

}

function create_box {
        local s=$1
        shift
        local rbuf=$1
        shift
        local ionconc=$1
        shift
        local list=("$@")

	out=()
	for i in "${!list[@]}";do
		out+=($(echo ${list[$i]}|awk -F "_" '{print $1}'))
	done

        cat <<EOF > get_req_water_in_MDbox.py
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
        chmod a+x get_req_water_in_MDbox.py

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
complex = loadPdbUsingSeq ${list[0]}.pdb { $(cat ${list[0]}.seq) }

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

#savepdb complex ${list[0]}-solv.pdb
saveamberparm complex ${list[0]}-solv.parm7 ${list[0]}-solv.rst7

quit
EOF
        tleap -s -f tleap.in > output

        nwat=$(calcwaterinparm ${list[0]}-solv.parm7)
        nions=$(python2.7 ./get_req_water_in_MDbox.py $nwat ${ionconc})

        #echo "No. of water and ions required in ${list[0]} : $nwat $nions"

	for i in "${!list[@]}";do
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
loadamberparams ${list[$i]}.frcmod
loadoff ${list[$i]}.lib

# load the coordinates and create the complex
complex = loadPdbUsingSeq ${list[$i]}.pdb { $(cat ${list[$i]}.seq) }

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

addions complex Na+ ${nions}
addions complex Cl- ${nions}

#savepdb complex ${list[$i]}-solv.pdb
saveamberparm complex ${list[$i]}-solv.parm7 ${list[$i]}-solv.rst7

quit
EOF
 		tleap -s -f tleap.in > output
		

		if [ ! -e "${list[$i]}-solv.parm7" ]; then
		    echo "tleap failed to write ${list[$i]}-solv.parm7"
		    exit 1
		fi
		if [ ! -e "${list[$i]}-solv.rst7" ]; then
		    echo "tleap failed to write ${list[$i]}-solv.rst7"
		    exit 1
		fi
		
        	mv ${list[$i]}-solv.parm7 ${out[$i]}_${s}.parm7; mv ${list[$i]}-solv.rst7  ${out[$i]}_${s}.rst7
	done
	rm -rf get_req_water_in_MDbox.py

}

function create_box_noseq {
        local s=$1
        shift
        local rbuf=$1
        shift
        local ionconc=$1
        shift
        local list=("$@")

        out=()
        for i in "${!list[@]}";do
                out+=($(echo ${list[$i]}|awk -F "_" '{print $1}'))
        done

        cat <<EOF > get_req_water_in_MDbox.py
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
        chmod a+x get_req_water_in_MDbox.py

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

        nwat=$(calcwaterinparm ${list[0]}-solv.parm7)
        nions=$(python2.7 ./get_req_water_in_MDbox.py $nwat ${ionconc})

        #echo "No. of water and ions required in ${list[0]} : $nwat $nions"

        for i in "${!list[@]}";do
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
loadamberparams ${list[$i]}.frcmod
loadoff ${list[$i]}.lib

# load the coordinates and create the complex
complex = loadPdb ${list[$i]}.pdb 

# create complex in solution for vdw+bonded transformation
solvateBox complex TIP4PEWBOX ${rbuf}
addions complex Na+ 0
addions complex Cl- 0

addions complex Na+ ${nions}
addions complex Cl- ${nions}

#savepdb complex ${list[$i]}-solv.pdb
saveamberparm complex ${list[$i]}-solv.parm7 ${list[$i]}-solv.rst7

quit
EOF
               tleap -s -f tleap.in > output

                mv ${list[$i]}-solv.parm7 ${out[$i]}_${s}.parm7; mv ${list[$i]}-solv.rst7  ${out[$i]}_${s}.rst7
        done
        rm -rf get_req_water_in_MDbox.py

}


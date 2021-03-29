function rename_in_pdb {
	local molname=$1
	local path=$2
	local ligname=$3
	local ligid=$4

	cat << EOF > renameinpdb.py
#!/usr/bin/env python3

if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Rename residue name in PDB file" )

    parser.add_argument("-p","--pdb",
                        help="PDB file",
                        type=str,
                        required=True)

    parser.add_argument("-r","--resid",
                        help="RESID in PDB to be renamed",
                        type=int,
                        required=True )

    parser.add_argument("-n","--resname",
                        help="Target residue name" ,
                        type=str,
                        default="XX",
                        required=False)

    parser.add_argument("-o","--output",
                        help="Output file basename" ,
                        type=str,
                        default="Renamed",
                        required=False)


    args = parser.parse_args()

    p = parmed.load_file(args.pdb)
    p.residues[args.resid].name = args.resname
    parmed.tools.writeCoordinates(p,"{}.pdb".format(args.output)).execute()
EOF

	chmod a+x renameinpdb.py
	./renameinpdb.py -p ${path}/${molname}.pdb -r ${ligid} -n ${ligname} -o ${molname}
	rm -rf renameinpdb.py
}

function mol2_atomname_fix {
	# in some mol2 files, antechamber changes names of atoms even with -an set to no. 
	# renaming mol2 atom names to pdb atom names for consistency

	local molname=$1

	cat << EOF > mol2atomnamefix.py
#!/usr/bin/env python3

if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Rename residue name in PDB file" )

    parser.add_argument("-p","--pdb",
                        help="PDB file",
                        type=str,
                        required=True)

    parser.add_argument("-m","--mol2",
                        help="RESID in PDB to be renamed",
                        type=str,
                        required=True )

    parser.add_argument("-o","--output",
                        help="Output file basename" ,
                        type=str,
                        default="Renamed",
                        required=False)


    args = parser.parse_args()

    p = parmed.load_file(args.pdb)
    m = parmed.load_file(args.mol2, structure=False)


    c=0
    for i in m.atoms:
        i.name=p.atoms[c].name
        c=c+1
    parmed.formats.Mol2File.write(m, "{}.mol2".format(args.output))
EOF
	chmod a+x mol2atomnamefix.py
	./mol2atomnamefix.py -p ${molname}.pdb -m ${molname}.mol2 -o ${molname} 

	rm -rf mol2atomnamefix.py

}


function preparePDBs {
	local larr=$1
	shift
	local path=$1
	shift
	local arr=("$@")

	local uniqueligs=($(echo "${arr[@]:0:${larr}}"))
	local uniquechgs=($(echo "${arr[@]:$(($larr)):$larr}"))

	truncate -s0 molname-ligname.mapping
	for i in "${!uniqueligs[@]}";do
		molname=${uniqueligs[$i]}
		ligname=$(printf "L%02d\n" $i)
		ligcharge=${uniquechgs[$i]}
		ligid=$(grep "${molname}\|L1" ${path}/${molname}.pdb  | awk -F " " '{print $5}'|sort --unique); ligid=$((ligid-1))
		rename_in_pdb "${molname}" "${path}" "${ligname}" "${ligid}"
		if ! command -v pdb4amber &> /dev/null; then printf "\n\n pdb4amber not present in path. pdb4amber is required for setup using inputformat=pdb" && exit 0; fi
		pdb4amber -i ${molname}.pdb -o ${molname}_com_dry.pdb -d >> output 2>&1
		mv ${molname}_com_dry_nonprot.pdb ${molname}_lig_dry.pdb

		if ! command -v antechamber &> /dev/null; then printf "\n\n antechamber not present in path. antechamber is required for setup using inputformat=pdb" && exit 0; fi
		antechamber -i ${molname}_lig_dry.pdb -fi pdb -o ${molname}_lig_dry.mol2 -fo mol2 -at "gaff" -pf y -an n -c bcc -nc ${ligcharge} -rn ${ligname} >> output 2>&1; sleep 1

		if grep -Fq 'Please check the total charge (-nc flag) and spin multiplicity (-m flag)' output; then
			printf "\n\n Check the charge provided in chargelist corresponding to ligand in pdb file ${path_to_input}/${system}/${molname}.pdb\n\n" && exit 0
                fi
		mol2_atomname_fix "${molname}_lig_dry"

		if ! command -v parmchk2 &> /dev/null; then printf "\n\n parmchk2 not present in path. parmchk2 is required for setup using inputformat=pdb" && exit 0; fi
                parmchk2 -i ${molname}_lig_dry.mol2 -o ${molname}_lig_dry.frcmod -f mol2 -s 1

		cat <<EOF > tleap.in

source leaprc.constph
set default PBradii mbondi3

source leaprc.protein.ff14SB
source leaprc.water.tip4pew
source leaprc.gaff2

loadamberparams ${molname}_lig_dry.frcmod
${ligname} = loadmol2 ${molname}_lig_dry.mol2
saveoff ${ligname} ${molname}_lig_dry.lib
quit

EOF
		if ! command -v tleap &> /dev/null; then printf "\n\n tleap not present in path. tleap is required for setup using inputformat=pdb" && exit 0; fi
		tleap -s -f tleap.in >> output 2>&1; sleep 1

		cp ${molname}_lig_dry.frcmod ${molname}_com_dry.frcmod
		cp ${molname}_lig_dry.lib    ${molname}_com_dry.lib
		echo "${molname} ${ligname}" >> molname-ligname.mapping
	done

}


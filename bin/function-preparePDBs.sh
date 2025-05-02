function rename_in_pdb {
	local molname=$1
	local path=$2
	local tligname=$3
	local iligname=$4

	cat << EOF > renameinpdb.py
#!/usr/bin/env python3

def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]

def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser     ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Rename residue name in PDB file" )

    parser.add_argument("-p","--pdb",
                        help="PDB file",
                        type=str,
                        required=True)

    parser.add_argument("-i","--iname",
                        help="RESNAME in PDB to be renamed",
                        type=str,
                        required=True )

    parser.add_argument("-t","--tname",
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
    print(f"renameinpdb.py is working on pdb {args.pdb}")
    for res in p.residues:
        if res.name == args.iname:
            res.name = args.tname
    #p.residues[args.iname].name = args.tname
    fh = open("%s.seq"%(args.output),"w") 
    seq = GetResSeq( p )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    fh.write(seqstr)
    fh.close()
    parmed.tools.writeCoordinates(p,"{}.pdb".format(args.output)).execute()

EOF

	chmod a+x renameinpdb.py
	./renameinpdb.py -p ${path}/${molname}.pdb -i ${iligname} -t ${tligname} -o ${molname}
	rm -rf renameinpdb.py
}


function preparePDBs {
	local larr=$1
	shift
	local path=$1
	shift
	local boxbuild=$1
	shift
	local uniqueligs=("$@")

	#local uniqueligs=($(echo "${arr[@]:0:${larr}}"))
	#local uniquechgs=($(echo "${arr[@]:$(($larr)):$larr}"))


	truncate -s0 molname-ligname.mapping
	for i in "${!uniqueligs[@]}";do

		molname=${uniqueligs[$i]}
		ligname=$(printf "L%02d\n" $i)

		rename_in_pdb "${molname}" "${path}" "${ligname}" "LIG"

		#### use pdb4amber to get s-s links and renumbered residues
		if [ "${boxbuild}" == 0 ]; then
			pdb4amber -i ${molname}.pdb -o tmp.pdb >> output 2>&1
		else
			pdb4amber -i ${molname}.pdb -o tmp.pdb -d >> output 2>&1
		fi

		mv tmp_sslink ${molname}_sslinks
		mv tmp.pdb ${molname}.pdb
		mv tmp_renum.txt ${molname}_res-renum
		rm tmp_nonprot.pdb
		####
		cp ${path}/${molname}_?.mol2 .
		cp ${path}/${molname}_?.frcmod .
		cp ${path}/${molname}_?.lib .
		sed -i "s/LIG/${ligname}/g" ${molname}_0.mol2 ${molname}_0.frcmod ${molname}_0.lib
		# number of non-standard residues estimated from the number of *_?.frcmod/*_?.lib files provided
		numnonstd=$(ls -l ${molname}_?.frcmod | wc -l)

		printf "${molname} ${ligname} ${numnonstd}\n" >> molname-ligname.mapping
	done

}


function getambermask {
	parm=$1; mask=$2
	cat << EOF > getmask.py
#!/usr/bin/python
import parmed
import argparse

def OpenParm( fname, xyz=None ):
    import parmed
    from parmed.constants import IFBOX
    if ".mol2" in fname:
        param = parmed.load_file( fname, structure=True )
        #help(param)
    else:
        param = parmed.load_file( fname,xyz=xyz )
        if xyz is not None:
            if ".rst7" in xyz:
                param.load_rst7(xyz)
    if param.box is not None:
        if abs(param.box[3]-109.471219)<1.e-4 and \\
           abs(param.box[4]-109.471219)<1.e-4 and \\
           abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param


def GetSelectedAtomIndices(param,maskstr):
    sele = []
    if len(maskstr) > 0:
        newmaskstr = maskstr.replace("@0","!@*")
        if len(newmaskstr) > 0:
            sele = [ param.atoms[i].idx for i in parmed.amber.mask.AmberMask( param, newmaskstr ).Selected() ]
    return sele


def ListToSelection(atomlist):
    alist = list(sorted(set(atomlist)))
    rs=[]
    if len(alist) > 0:
        rs = [ (alist[0],alist[0]) ]
        for a in alist[1:]:
            if a == rs[-1][1]+1:
                rs[-1] = ( rs[-1][0], a )
            else:
                rs.append( (a,a) )
    sarr = []
    for r in rs:
        if r[0] != r[1]:
            sarr.append( "%i-%i"%(r[0]+1,r[1]+1) )
        else:
            sarr.append( "%i"%(r[0]+1) )
    sele = "@0"
    if len(sarr) > 0:
        sele = "@" + ",".join(sarr)
    return sele



if __name__ == "__main__":
    parser = argparse.ArgumentParser \\
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""Converts an arbitrary amber mask to an atom selection mask
""" )
    
    parser.add_argument \\
        ("-p","--parm",
         help="parm7 file",
         type=str,
         required=True )


    parser.add_argument \\
        ("-w","--whole-molecule",
         help="if present, then selecting any portion of a molecule selects the whole molecule",
         action='store_true' )

    
    parser.add_argument \\
        ("-r","--byresidue",
         help="if present, then return a list of matching residues",
         action='store_true' )


    parser.add_argument \\
        ("-n","--norange",
         help="if present, then list each match individually -- as opposed to using dashes to represent ranges",
         action='store_true' )

    
    parser.add_argument \\
        ("-c","--crd",
         help="input restart file",
         type=str,
         default=None,
         required=False )
    
    parser.add_argument \\
        ("-m","--mask",
         help="amber mask selection of atoms to polarize",
         type=str,
         required=True )
    
    args = parser.parse_args()
    

    parmfile = args.parm
    crd7file = args.crd
    maskstr  = args.mask
    maskstr = maskstr.replace('\\\',"").replace('\\\',"").replace('\\\',"").replace('\\\',"").replace('\\\',"")

    if crd7file is not None:
        try:
            parm = OpenParm( parmfile, xyz=crd7file )
        except:
            parm = OpenParm( parmfile )
    else:
        parm = OpenParm( parmfile )


        
    atoms = GetSelectedAtomIndices( parm, maskstr )

    if args.whole_molecule:
       smols=[]
       iat = 0
       for imol,mnat in enumerate(parm.parm_data["ATOMS_PER_MOLECULE"]):
           for i in range(mnat):
               if iat in atoms:
                  smols.append( imol )
               iat += 1
       smols=list(set(smols))
       smols.sort()
       atoms = []
       iat=0
       for imol,mnat in enumerate(parm.parm_data["ATOMS_PER_MOLECULE"]):
           for i in range(mnat):
               if imol in smols:
                  atoms.append(iat)
               iat += 1 

    res = []
    for a in atoms:
        res.append( parm.atoms[a].residue.idx )
    res = list(set(res))
    res.sort()

    if args.byresidue:
        if args.norange:
            s = ":%s"%( ",".join( ["%i"%(x+1) for x in res ] ) )
        else:
            s = ListToSelection( res )
            s = s.replace("@",":")
    else:
        if args.norange:
            s = "@%s"%( ",".join( ["%i"%(x+1) for x in atoms ] ) )
        else:
            s = ListToSelection( atoms )
            
    print(s)

EOF
	chmod a+x getmask.py
	echo $(./getmask.py -p ${parm} -m ${mask})

}



if [ "${setupmode}" == 0 ]; then

        if [ "${protocol}" == "unified" ]; then
                cd $path/$system/setup
			if [ "${ticalc}" != "asfe" ]; then
                        	for i in "${!translist[@]}";do
                                	stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                	for s in ${slist[@]}; do
                                        	mkdir -p ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
						if [ "${twostate}" == "false" ]; then
							cp ${stA}~${stB}_${s}.parm7  ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7
							cp ${stA}~${stB}_${s}.rst7   ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/stateA.rst7
						else
                                        		cp ${stA}~${stB}-1_${s}.parm7  ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7
                                        		cp ${stA}~${stB}-1_${s}.rst7   ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/stateA.rst7
                                        		cp ${stA}~${stB}-2_${s}.rst7   ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/stateB.rst7
						fi

						scmask1=$(cat "${stA}~${stB}".scmask1)
						scmask2=$(cat "${stA}~${stB}".scmask2)
						timask1=$(cat "${stA}~${stB}".timask1)
						timask2=$(cat "${stA}~${stB}".timask2)
						noshakemask="${timask1},${timask2:1:${#timask2}}"
						#printf "\n ${timask1} ${timask2} ${scmask1} ${scmask2} ${noshakemask} \n"
						# get succint ambermasks
						if [ "${scmask1}" != '""' ]; then
							scmask1=$(getambermask "${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7" "${scmask1}")
							scmask1="'${scmask1}'"
						fi
						if [ "${scmask2}" != '""' ]; then
							scmask2=$(getambermask "${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7" "${scmask2}")
							scmask2="'${scmask2}'"
						fi
						timask1=$(getambermask "${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7" "${timask1}"); timask1="'${timask1}'"
						timask2=$(getambermask "${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7" "${timask2}"); timask2="'${timask2}'"
						noshakemask=$(getambermask "${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/unisc.parm7" "${noshakemask}"); noshakemask="'${noshakemask}'"

						#printf "\n ${timask1} ${timask2} ${scmask1} ${scmask2} ${noshakemask} \n"

                                        	cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
                                                	if [ "${ticalc}" == "rbfe" ]; then
                                                        	writetemplate_rbfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate}
                                                	else
                                                        	writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${boxbufaq}
                                                	fi

							sh TEMPLATE.sh; sleep 1
						cd $path/${system}/setup



						if [ "${s}" == "aq" ]; then bidirection=${bidirection_aq}; else bidirection=${bidirection_com}; fi
						if [ "${bidirection}" == "true" ] && [ "${twostate}" == "false" ]; then
							mkdir -p ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}
							if [ "${s}" == "aq" ] || [ "${twostate}" == "false" ]; then
								cp ${stB}~${stA}_${s}.parm7  ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7
								cp ${stB}~${stA}_${s}.rst7   ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/stateA.rst7
							fi
							scmask1=$(cat "${stB}~${stA}".scmask1)
							scmask2=$(cat "${stB}~${stA}".scmask2)
							timask1=$(cat "${stB}~${stA}".timask1)
							timask2=$(cat "${stB}~${stA}".timask2)
							noshakemask="${timask1},${timask2:1:${#timask2}}"
							# get succint ambermasks
							if [ "${scmask1}" != '""' ]; then
								scmask1=$(getambermask "${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7" "${scmask1}")
								scmask1="'${scmask1}'"
							fi
							if [ "${scmask2}" != '""' ]; then
								scmask2=$(getambermask "${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7" "${scmask2}")
								scmask2="'${scmask2}'"
							fi
							timask1=$(getambermask "${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7" "${timask1}"); timask1="'${timask1}'"
							timask2=$(getambermask "${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7" "${timask2}"); timask2="'${timask2}'"
							noshakemask=$(getambermask "${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/unisc.parm7" "${noshakemask}"); noshakemask="'${noshakemask}'"

							cd ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}
                                                        	if [ "${ticalc}" == "rbfe" ]; then
                                                                	writetemplate_rbfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate}
                                                        	else
                                                                	writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${boxbufaq}
                                                        	fi

                                                        	sh TEMPLATE.sh; sleep 1
                                                	cd $path/${system}/setup


						fi
					done
				done
			else
				for i in "${!translist[@]}";do
					for s in ${slist[@]}; do
						mkdir -p ${path}/${system}/${protocol}/run/${translist[$i]}/${s}
						cp ${translist[$i]}_${s}.parm7 ${path}/${system}/${protocol}/run/${translist[$i]}/${s}/unisc.parm7
						cp ${translist[$i]}_${s}.rst7  ${path}/${system}/${protocol}/run/${translist[$i]}/${s}/stateA.rst7

						scmask1="':1'"
						scmask2="''"
						timask1="':1'"
						timask2="''"
						noshakemask="':1'"

						cd ${path}/${system}/${protocol}/run/${translist[$i]}/${s}
							 writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${boxbufaq} 
							 sh TEMPLATE.sh; sleep 1
						cd $path/${system}/setup
					done
				done
			fi



			for i in "${!translist[@]}";do
                                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
				for s in ${slist[@]}; do
					if [ "${ticalc}" != "asfe" ]; then
						folders=("${stA}~${stB}" "${stB}~${stA}")
					else
						folders=("${translist[$i]}")
					fi
					for dir in "${folders[@]}"; do
						if [ -d "${path}/${system}/${protocol}/run/${dir}/${s}" ]; then
							cd ${path}/${system}/${protocol}/run/${dir}/${s}

                                				# if hmr=true set timestep to 4fs
                                				if [ "${hmr}" == "true" ]; then
                                        				sed -i '/dt.*.=.*.*/c\dt              = 0.004' inputs/*_ti.mdin
                                				fi

								# if notrajecory=true, set ntwx=0 in input files
                                				if [ "${notrajectory}" == "true" ]; then
									sed -i 's/ntwx.*/ntwx            = 0/' inputs/*.mdin
                                				fi

                                				# if repex=false, alter input files and slurm files
                                				if [ "${repex}" == "false" ]; then
                                        				sed -i 	-e '/numexchg/d' \
										-e '/gremd_acyc/d' \
										-e 's/ntwx .*/ntwx            = 10000/' \
										-e 's/ntwr .*/ntwx            = 5000/' \
										-e 's/ntpr .*/ntpr            = 1000/' \
										-e 's/bar_intervall .*/bar_intervall   = 1000/' \
										inputs/*_ti.mdin
                                        				sed -i 	-e 's/ -rem 3 -remlog remt${trial}.log//g' \
										-e 's/running replica ti/running regular ti/g' \
										run_alltrials.slurm
                                				fi

								for(( t=1;t<=${ntrials};t++));do
									mkdir -p t${t}
									cp current/*_init.rst7 t${t}/
									if [ "${twostate}" == "true" ]; then
										sed "s/current/t${t}/g" inputs/eqpre1P0.groupfile       > inputs/t${t}_eqpre1P0.groupfile
										sed "s/current/t${t}/g" inputs/eqpre2P0.groupfile       > inputs/t${t}_eqpre2P0.groupfile
										sed "s/current/t${t}/g" inputs/eqP0.groupfile           > inputs/t${t}_eqP0.groupfile
										sed "s/current/t${t}/g" inputs/eqNTP4.groupfile         > inputs/t${t}_eqNTP4.groupfile
										sed "s/current/t${t}/g" inputs/eqV.groupfile            > inputs/t${t}_eqV.groupfile
										sed "s/current/t${t}/g" inputs/eqP.groupfile            > inputs/t${t}_eqP.groupfile
										sed "s/current/t${t}/g" inputs/eqA.groupfile            > inputs/t${t}_eqA.groupfile
										if [ "${s}" == "com" ]; then
										        sed "s/current/t${t}/g" inputs/eqProt2.groupfile        > inputs/t${t}_eqProt2.groupfile
										        sed "s/current/t${t}/g" inputs/eqProt1.groupfile        > inputs/t${t}_eqProt1.groupfile
										        sed "s/current/t${t}/g" inputs/eqProt05.groupfile       > inputs/t${t}_eqProt05.groupfile
										        sed "s/current/t${t}/g" inputs/eqProt025.groupfile      > inputs/t${t}_eqProt025.groupfile
										        sed "s/current/t${t}/g" inputs/eqProt01.groupfile       > inputs/t${t}_eqProt01.groupfile
										        sed "s/current/t${t}/g" inputs/eqProt0.groupfile        > inputs/t${t}_eqProt0.groupfile
										fi
									fi
									sed "s/current/t${t}/g" inputs/eqATI.groupfile          > inputs/t${t}_eqATI.groupfile
									sed "s/current/t${t}/g" inputs/preTI.groupfile          > inputs/t${t}_preTI.groupfile
									sed "s/current/t${t}/g" inputs/ti.groupfile             > inputs/t${t}_ti.groupfile
								
								done


                                        		cd $path/${system}/setup
						fi
                                	done
					if [ -d ${path}/${system}/${protocol}/run/${stA}~${stB}/${s} ] && [ -d ${path}/${system}/${protocol}/run/${stB}~${stA}/${s} ]; then
						mkdir -p ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/forward
						mv ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/* ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/forward/ 2>/dev/null
						mkdir -p ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/reverse
						mv ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}/* ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/reverse/ 2>/dev/null
						rm -rf ${path}/${system}/${protocol}/run/${stB}~${stA}/${s}
					fi
				done
				printf "Done with ${translist[$i]}...\n"
				if [ -d "${path}/${system}/${protocol}/run/${stB}~${stA}" ] && [ -z "$(ls -A ${path}/${system}/${protocol}/run/${stB}~${stA})" ]; then rm -rf ${path}/${system}/${protocol}/run/${stB}~${stA}; fi
                        done
                cd ${path}

        fi

fi
# END of setupmode=0


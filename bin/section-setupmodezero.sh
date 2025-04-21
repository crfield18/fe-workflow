function getambermask {
	parm=$1; mask=$2
	cat << EOF > getmask.py
#!/usr/bin/env python3
import parmed
import argparse

def OpenParm( fname, xyz=None ):
    import parmed
    try:
        from parmed.constants import PrmtopPointers
        IFBOX = PrmtopPointers.IFBOX
    except: 
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
	echo $(python3 getmask.py -p ${parm} -m ${mask})

}

printf "\n ***************************************\n"
printf " *                                     *\n"
printf " *  Entering SetupZero                 *\n"
printf " *                                     *\n"
printf " ***************************************\n\n"

if [ "${setupmode}" == 0 ]; then

        if [ "${protocol}" == "unified" ]; then
                if [ ${equil_type} == 1 ]; then

                cat << EOF_runalltrials > run_all_equil.slurm
translist=(${translist[@]})
logfile=\${PWD}/${system}/submitted_equil_jobs.log

touch \${logfile}

for i in "\${!translist[@]}";do
        for dir in aq com; do
                cd \${PWD}/${system}/unified/run/\${translist[\$i]}/\${dir}
                job_id="\${translist[\$i]}_\${dir}_run_equil.slurm"
                if ! grep -q "\${job_id}" \${logfile}; then
                        output=\$(sbatch run_equilibration.slurm)
                        if [[ \${output} == *"Submitted batch job"* ]]; then
                                echo "---->Submitted job: \${job_id}"
                                echo "\${job_id}" >> \${logfile}
                        else
                                echo "Failed to submit job \${job_id}:"
                                echo "**********REASON************"
                                echo "\${output}"
                                echo "******************************"
                                exit 1
                        fi
                else
                        echo "Job \${job_id} has already been submitted."
                fi
                cd -
        done 
done
EOF_runalltrials


                cat << EOF_runalltrials > run_all_prod.slurm
translist=(${translist[@]})
logfile=\${PWD}/${system}/submitted_production_jobs.log

touch \${logfile}

for i in "\${!translist[@]}";do
        for dir in aq com; do
                cd \${PWD}/${system}/unified/run/\${translist[\$i]}/\${dir}
                for j in \$(seq 1 1 ${ntrials}); do
                        job_id="\${translist[\$i]}_\${dir}_run_trial_\${j}.slurm"
                        if ! grep -q "\${job_id}" \${logfile}; then
                                output=\$(sbatch run_production_trial_\${j}.slurm)
                                if [[ \${output} == *"Submitted batch job"* ]]; then
                                        echo "---->Submitted job: \${job_id}"
                                        echo "\${job_id}" >> \${logfile}
                                else
                                        echo "Failed to submit job \${job_id}:"
                                        echo "**********REASON************"
                                        echo "\${output}"
                                        echo "******************************"
                                        exit 1
                                fi
                        else
                                echo "Job \${job_id} has already been submitted."
                        fi

                done
                cd -
        done 
done
EOF_runalltrials

                elif [ ${equil_type} == 2 ]; then
                cat << EOF_runalltrials > run_all_trials.slurm
translist=(${translist[@]})
logfile=\${PWD}/${system}/submitted_jobs.log

touch \${logfile}

for i in "\${!translist[@]}";do
EOF_runalltrials
if [ ${combine_aq} == "false" ]; then
        printf "for dir in aq com; do\n" >> run_all_trials.slurm
else
        printf "for dir in com; do\n" >> run_all_trials.slurm
fi
cat << EOF_runalltrials >> run_all_trials.slurm
                cd \${PWD}/${system}/unified/run/\${translist[\$i]}/\${dir}
                for j in \$(seq 1 1 ${ntrials}); do
                        job_id="\${translist[\$i]}_\${dir}_run_trial_\${j}.slurm"
                        if ! grep -q "\${job_id}" \${logfile}; then
                                output=\$(sbatch run_trial_\${j}.slurm)
                                if [[ \${output} == *"Submitted batch job"* ]]; then
                                        echo "---->Submitted job: \${job_id}"
                                        echo "\${job_id}" >> \${logfile}
                                else
                                        echo "Failed to submit job \${job_id}:"
                                        echo "**********REASON************"
                                        echo "\${output}"
                                        echo "******************************"
                                        exit 1
                                fi
                        else
                                echo "Job \${job_id} has already been submitted."
                        fi

                done
                cd -
        done 
done
EOF_runalltrials
                fi
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
                                                                if [[ -v override_lambda ]]; then 
                                                                        lams=()
                                                                        read_lambda_schedule "${path}/${override_lambda}/${stA}~${stB}_${s}_ar_${nlambda}.txt" lams
                                                                fi
                                                        	writetemplate_rbfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${ntwx_equil} ${ntwx} ${ntwr} ${ntpr} ${equil_type} ${source_header} ${max_dt} ${nnodes} ${ntwx_ep} ${combine_aq} ${autoimage}
                                                	else
                                                        	writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${nnodes} ${ntwx_equil} ${ntwx} ${ntwr} ${ntpr} ${source_header} ${max_dt}
                                                	fi

							bash TEMPLATE.sh; sleep 1
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
                                                                        if [[ -v override_lambda ]]; then 
                                                                                lams=()
                                                                                read_lambda_schedule "${path}/${override_lambda}/${stB}~${stA}_${s}_ar_${nlambda}.txt" lams
                                                                        fi
                                                                	writetemplate_rbfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${ntwx_equil} ${ntwx} ${ntwr} ${ntpr} ${equil_type} ${source_header} ${max_dt} ${nnodes} ${ntwx_ep} ${combine_aq} ${autoimage}
                                                        	else
                                                                	writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${ntwx_equil} ${ntwx} ${ntwr} ${ntpr} ${source_header} ${max_dt}
                                                        	fi

                                                        	bash TEMPLATE.sh; sleep 1
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
							 writetemplate_rsfe $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp ${translist[$i]} $s ${twostate} ${ntwx_equil} ${ntwx} ${ntwr} ${ntpr} ${source_header} ${max_dt}
							 bash TEMPLATE.sh; sleep 1
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
                                                                if [[ "${full_sc}" == "true" ]]; then 
                                                                    # We assume that all mdinputs have the same timasks and scmasks (as they should)
                                                                    timask1=`awk -F\' '/timask1\s*=/{print $2}' inputs/${lams[0]}_ti.mdin`
                                                                    timask2=`awk -F\' '/timask2\s*=/{print $2}' inputs/${lams[0]}_ti.mdin`
								    for f in inputs/*mdin; do
									if [ -e "${f}" ]; then
									    sed -i -e "/^\s*scmask1\s*=/ s/'[^']*'/'$timask1'/" \
										-e "/^\s*scmask2\s*=/ s/'[^']*'/'$timask2'/" "${f}"
									fi
								    done
								    for f in inputs/*mdin.template; do
									if [ -e "${f}" ]; then
									    sed -i -e "/^\s*scmask1\s*=/ s/'[^']*'/'$timask1'/" \
										-e "/^\s*scmask2\s*=/ s/'[^']*'/'$timask2'/" "${f}"
									fi
								    done
                                                                fi
                                                                if [ "${equil_type}" == "1" ] ; then 
                                                                    mkdir -p equil
								    for f in inputs/eqp*.groupfile \
										 inputs/eqP*.groupfile \
										 inputs/eqN*.groupfile \
										 inputs/eqV*.groupfile \
										 inputs/eqB*.groupfile \
										 inputs/eqA*.groupfile; do
									if [ -e "${f}" ]; then
									    dname=$(dirname $f)
									    bname=$(basename $f)
									    sed "s/current/equil/g" ${f} > ${dname}/equil_${bname}
									fi
								    done
                                                                        # sed "s/current/equil/g" inputs/eqpre1P0.groupfile       > inputs/equil_eqpre1P0.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqpre2P0.groupfile       > inputs/equil_eqpre2P0.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqP0.groupfile           > inputs/equil_eqP0.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqP1.groupfile           > inputs/equil_eqP1.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqP2.groupfile           > inputs/equil_eqP2.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqNTP4.groupfile         > inputs/equil_eqNTP4.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqV.groupfile            > inputs/equil_eqV.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqP.groupfile            > inputs/equil_eqP.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqA.groupfile            > inputs/equil_eqA.groupfile
                                                                        # if [ "${s}" == "com" ]; then
                                                                        #         sed "s/current/equil/g" inputs/eqProt2.groupfile        > inputs/equil_eqProt2.groupfile
                                                                        #         sed "s/current/equil/g" inputs/eqProt1.groupfile        > inputs/equil_eqProt1.groupfile
                                                                        #         sed "s/current/equil/g" inputs/eqProt05.groupfile       > inputs/equil_eqProt05.groupfile
                                                                        #         sed "s/current/equil/g" inputs/eqProt025.groupfile      > inputs/equil_eqProt025.groupfile
                                                                        #         sed "s/current/equil/g" inputs/eqProt01.groupfile       > inputs/equil_eqProt01.groupfile
                                                                        #         sed "s/current/equil/g" inputs/eqProt0.groupfile        > inputs/equil_eqProt0.groupfile
                                                                        # fi
                                                                        # sed "s/current/equil/g" inputs/eqATI.groupfile          > inputs/equil_eqATI.groupfile
                                                                        # sed "s/current/equil/g" inputs/eqBTI.groupfile          > inputs/equil_eqBTI.groupfile
                                                                elif [ "${equil_type}" == "2" ] ; then 
                                                                        for (( t=1;t<=${ntrials};t++));do
                                                                            mkdir -p equil${t}
									    for f in inputs/eqp*.groupfile \
											 inputs/eqP*.groupfile \
											 inputs/eqN*.groupfile \
											 inputs/eqV*.groupfile \
											 inputs/eqB*.groupfile \
											 inputs/eqA*.groupfile; do
										if [ -e "${f}" ]; then
										    dname=$(dirname $f)
										    bname=$(basename $f)
										    sed "s/current/equil${t}/g" ${f} > ${dname}/equil${t}_${bname}
										fi
									    done
                                                                            # sed "s/current/equil${t}/g" inputs/eqpre1P0.groupfile       > inputs/equil${t}_eqpre1P0.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqpre2P0.groupfile       > inputs/equil${t}_eqpre2P0.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqP0.groupfile           > inputs/equil${t}_eqP0.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqP1.groupfile           > inputs/equil${t}_eqP1.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqP2.groupfile           > inputs/equil${t}_eqP2.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqNTP4.groupfile         > inputs/equil${t}_eqNTP4.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqV.groupfile            > inputs/equil${t}_eqV.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqP.groupfile            > inputs/equil${t}_eqP.groupfile
                                                                            # sed "s/current/equil${t}/g" inputs/eqA.groupfile            > inputs/equil${t}_eqA.groupfile
                                                                            # if [ "${s}" == "com" ]; then
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt2.groupfile        > inputs/equil${t}_eqProt2.groupfile
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt1.groupfile        > inputs/equil${t}_eqProt1.groupfile
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt05.groupfile       > inputs/equil${t}_eqProt05.groupfile
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt025.groupfile      > inputs/equil${t}_eqProt025.groupfile
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt01.groupfile       > inputs/equil${t}_eqProt01.groupfile
                                                                            #     sed "s/current/equil${t}/g" inputs/eqProt0.groupfile        > inputs/equil${t}_eqProt0.groupfile
                                                                            # fi
                                                                            # sed "s/current/equil${t}/g" inputs/eqATI.groupfile          > inputs/equil${t}_eqATI.groupfile
                                                                        done
                                                                fi
                                                                

								for(( t=1;t<=${ntrials};t++));do
								    mkdir -p t${t}
								    if [ -e inputs/preTI.groupfile ]; then
									sed "s/current/t${t}/g" inputs/preTI.groupfile          > inputs/t${t}_preTI.groupfile
								    fi
                                                                    if [ "${equil_type}" == "2" ]; then
									if [ -e inputs/t${t}_preTI.groupfile ]; then
                                                                            sed "s/equil/equil${t}/g"  inputs/t${t}_preTI.groupfile         > inputs/t${t}_preTI.groupfiletmp
                                                                            cp inputs/t${t}_preTI.groupfiletmp inputs/t${t}_preTI.groupfile
									fi
                                                                    fi
								    if [ -e inputs/ti.groupfile ]; then
									sed "s/current/t${t}/g" inputs/ti.groupfile             > inputs/t${t}_ti.groupfile
								    fi
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
################################################
#starting to fix the box size to be the same as
#a post-process after generating all the box
#after the workflow. The step is to find the biggest box
#dimension in each environment (biggest dimension, not the biggest box)
#and replace all the box with the biggest dimension
#It's better to build all the edges (including self-ACES to 
#create the reservoirs) at once, to avoid manual change
#to the eqNTP4 stage. Recommendation: build all the
#self-ACES transformations along with the edges
#you want, which will definitely have exact same
#box for future use.
################################################

      echo "fix box"
      cd $path/$system/setup
     if [ "${ticalc}" == "asfe" ]; then
     if [ -f "fix_box_aq_size.txt" ]
     then
             read -r max_a < fix_box_aq_size.txt
             read -r max_b < <(sed -n 2p fix_box_aq_size.txt)
             read -r max_c < <(sed -n 3p fix_box_aq_size.txt)
     else     
                        max_a=$(awk '{print $1}' *_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
                        max_b=$(awk '{print $2}' *_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
                        max_c=$(awk '{print $3}' *_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
                        echo $max_a >> fix_box_aq_size.txt
                        echo $max_b >> fix_box_aq_size.txt
                        echo $max_c >> fix_box_aq_size.txt
     fi
  for aq_file in *_aq*rst7;do

     FILENAME_WITHOUT_EXT=$(echo "$aq_file" | sed 's/_aq.rst7$//')
     
     for num in 1;do
     cat <<EOF > fix_box_aq_cpptraj.in
     parm ${FILENAME_WITHOUT_EXT}_aq.parm7
     trajin ${FILENAME_WITHOUT_EXT}_aq.rst7
     trajout ${FILENAME_WITHOUT_EXT}_aq.pdb pdb include_ep
     go
     quit
EOF

     cpptraj -i fix_box_aq_cpptraj.in

           # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "source leaprc.phosaa14SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        elif [ "${pff}" == "ff19SB" ]; then
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "source leaprc.phosaa19SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadamberparams frcmod.ff19SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        elif [ "${pff}" == "nucleic" ]; then   
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "source leaprc.RNA.OL3\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "source leaprc.DNA.OL21\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        fi


        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadoff tip4pewbox.off\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                boxkey="TIP3PBOX"
        elif [ "${wm}" == "opc" ]; then
                printf "source leaprc.water.opc\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadamberparams frcmod.opc\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadamberparams frcmod.ionslm_hfe_opc\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                printf "loadoff opcbox.off\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
                boxkey="OPC3BOX"
        fi


        printf "loadamberparams ${FILENAME_WITHOUT_EXT}_0.frcmod\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        printf "loadoff ${FILENAME_WITHOUT_EXT}_0.lib\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in

        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
#        if [ "${load}" == "pdb" ]; then
                printf "x = loadPdb ${FILENAME_WITHOUT_EXT}_aq.pdb\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
#        elif [ "${load}" == "pdbseq" ]; then
#                printf "x = loadPdbUsingSeq ${inpfile}.pdb { $(cat ${inpfile}.seq) }\n" >> tleap.in
#        else
#                printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in
#        fi

        # add S-S cysteine linkkages if present
#        if [ -f ${inpfile}_sslinks ] && [ "$(cat ${inpfile}_sslinks | wc -l)" -gt 0 ]; then
#                while read line; do
#                        IFS=' ' read -ra args <<< $line
#                        printf "bond x.${args[0]}.SG x.${args[1]}.SG\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
#                done < ${inpfile}_sslinks
#        fi

        # build box and neutralize with Na+ Cl-
#        if [ "${boxbuild}" == 0 ] && [ "${s}" == "com" ]; then
#                printf "setbox x vdw \n" >> tleap.in
#        else
        printf "set x box {$max_a $max_b $max_c}\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        printf "saveamberparm x fix_box_${FILENAME_WITHOUT_EXT}_aq.parm7 fix_box_${FILENAME_WITHOUT_EXT}_aq.rst7\n\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in
        printf "quit\n" >> fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in

        tleap -s -f fix_box_${FILENAME_WITHOUT_EXT}_aq_tleap.in >> fix_box_aq_log

	# Reperform HMR if needed
	if [ "${hmr}" == "true" ]; then
		if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
			cat <<EOFM > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOFM
	        parmed -i hmr.in -p fix_box_${FILENAME_WITHOUT_EXT}_aq.parm7 -c fix_box_${FILENAME_WITHOUT_EXT}_aq.rst7 >> output 2>&1
		mv hmr.parm7 fix_box_${FILENAME_WITHOUT_EXT}_aq.parm7; mv hmr.rst7  fix_box_${FILENAME_WITHOUT_EXT}_aq.rst7
	fi

        # Redo fix box if mdboxshape is oct 
        if [ "${mdboxshape}" == "oct" ]; then
                mol="fix_box_${FILENAME_WITHOUT_EXT}_aq"
                boxinfo=$(tail -n 1 ${mol}.rst7)
                xval=$(echo "$boxinfo" | awk '{print $1}')
                mv ${mol}.parm7 temp_${mol}.parm7
                mv ${mol}.rst7 temp_${mol}.rst7
                cat > fix_parm7_oct.in <<EOF 
parm temp_${mol}.parm7
parmbox truncoct x ${xval} 
parmwrite out ${mol}.parm7
run 
EOF

                cat > fix_rst7_oct.in <<EOF 
parm temp_${mol}.parm7
trajin temp_${mol}.rst7
box x ${xval} y ${xval} z ${xval} truncoct
trajout ${mol}.rst7
run 
EOF

                cpptraj -i fix_parm7_oct.in
                cpptraj -i fix_rst7_oct.in
                rm temp_${mol}.parm7
                rm temp_${mol}.rst7

        fi

  done

              cp fix_box_${FILENAME_WITHOUT_EXT}_aq.parm7  ${path}/${system}/${protocol}/run/${FILENAME_WITHOUT_EXT}/aq/unisc.parm7
              cp fix_box_${FILENAME_WITHOUT_EXT}_aq.rst7   ${path}/${system}/${protocol}/run/${FILENAME_WITHOUT_EXT}/aq/stateA.rst7
        for(( t=1;t<=${ntrials};t++));do
        cp ${path}/${system}/${protocol}/run/${FILENAME_WITHOUT_EXT}/aq/stateA.rst7 ${path}/${system}/${protocol}/run/${FILENAME_WITHOUT_EXT}/aq/t${t}/0.00000000_init.rst7
	done
done

else #######


      #FIXME, read in the last line of *rst7 file and find the maximum box
     if [ -f "fix_box_aq_size.txt" ]
     then
	     read -r max_a < fix_box_aq_size.txt
	     read -r max_b < <(sed -n 2p fix_box_aq_size.txt)
	     read -r max_c < <(sed -n 3p fix_box_aq_size.txt)
     else
			max_a=$(awk '{print $1}' *-1_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
			max_b=$(awk '{print $2}' *-1_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
			max_c=$(awk '{print $3}' *-1_aq*rst7 | tail -n 1 | sort -n | tail -n 1)
			echo $max_a >> fix_box_aq_size.txt
			echo $max_b >> fix_box_aq_size.txt
			echo $max_c >> fix_box_aq_size.txt
     fi
   for aq_file in *-1_aq*rst7;do
    if [[ "${aq_file}" == "fix_box"* ]]; then
        continue
    fi
     FILENAME_WITHOUT_EXT=$(echo "$aq_file" | sed 's/-[0-9]\+_aq.rst7$//')
     BEFORE_TILDE=$(echo "$FILENAME_WITHOUT_EXT" | cut -d'~' -f1)
     AFTER_TILDE=$(echo "$FILENAME_WITHOUT_EXT" | cut -d'~' -f2)
     
     for num in 1 2;do
     cat <<EOF > fix_box_aq_cpptraj.in
     parm ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.parm7
     trajin ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.rst7
     trajout ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.pdb pdb include_ep
     go
     quit
EOF
     
     cpptraj -i fix_box_aq_cpptraj.in
     
        # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "source leaprc.phosaa14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        elif [ "${pff}" == "ff19SB" ]; then
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "source leaprc.phosaa19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadamberparams frcmod.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        elif [ "${pff}" == "nucleic" ]; then   
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "source leaprc.RNA.OL3\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "source leaprc.DNA.OL21\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        fi
	

        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadoff tip4pewbox.off\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                boxkey="TIP3PBOX"
        elif [ "${wm}" == "opc" ]; then
                printf "source leaprc.water.opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadamberparams frcmod.opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadamberparams frcmod.ionslm_hfe_opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                printf "loadoff opcbox.off\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
                boxkey="OPCBOX"
        fi

        # check and load non-standard residue parameter files
        #i=0
#	numnonstd=$(($numnonstd+0))
#        while [ "$i" -lt "${numnonstd}" ]; do
#                printf "loadamberparams ${lig1}_${i}.frcmod\n" >> tleap.in
#                printf "loadoff ${lig1}_${i}.lib\n" >> tleap.in
#                i=$(($i+1))
#        done

        # load ligand 2 parameter files
#        printf "loadamberparams ${lig2}_0.frcmod\n" >> tleap.in
#        printf "loadoff ${lig2}_0.lib\n" >> tleap.in

        printf "loadamberparams ${BEFORE_TILDE}_0.frcmod\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        printf "loadoff ${BEFORE_TILDE}_0.lib\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        printf "loadamberparams ${AFTER_TILDE}_0.frcmod\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        printf "loadoff ${AFTER_TILDE}_0.lib\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in

        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
#        if [ "${load}" == "pdb" ]; then
                printf "x = loadPdb ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.pdb\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
#        elif [ "${load}" == "pdbseq" ]; then
#                printf "x = loadPdbUsingSeq ${inpfile}.pdb { $(cat ${inpfile}.seq) }\n" >> tleap.in
#        else
#                printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in
#        fi

        # add S-S cysteine linkkages if present
#        if [ -f ${inpfile}_sslinks ] && [ "$(cat ${inpfile}_sslinks | wc -l)" -gt 0 ]; then
#                while read line; do
#                        IFS=' ' read -ra args <<< $line
#                        printf "bond x.${args[0]}.SG x.${args[1]}.SG\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
#                done < ${inpfile}_sslinks
#        fi

        # build box and neutralize with Na+ Cl-
#        if [ "${boxbuild}" == 0 ] && [ "${s}" == "com" ]; then
#                printf "setbox x vdw \n" >> tleap.in
#        else
        printf "set x box {$max_a $max_b $max_c}\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
	printf "saveamberparm x fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.parm7 fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.rst7\n\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        printf "quit\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in
        
        tleap -s -f fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq_tleap.in >> fix_box_aq_log

	# Reperform HMR if needed
	if [ "${hmr}" == "true" ]; then
		if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
			cat <<EOFM > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOFM
	        parmed -i hmr.in -p fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.parm7 -c fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.rst7 >> output 2>&1
		mv hmr.parm7 fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.parm7; mv hmr.rst7  fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq.rst7
                #parmed -i hmr.in -p fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_aq.parm7 -c fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_aq.rst7 >> output 2>&1
                #mv hmr.parm7 fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_aq.parm7; mv hmr.rst7  fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_aq.rst7
	fi

        # Redo fix box if mdboxshape is oct 
        if [ "${mdboxshape}" == "oct" ]; then
		mol="fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_aq"
                boxinfo=$(tail -n 1 ${mol}.rst7)
                xval=$(echo "$boxinfo" | awk '{print $1}')
                mv ${mol}.parm7 temp_${mol}.parm7
                mv ${mol}.rst7 temp_${mol}.rst7
		cat > fix_parm7_oct.in <<EOF 
parm temp_${mol}.parm7
parmbox truncoct x ${xval} 
parmwrite out ${mol}.parm7
run 
EOF

		cat > fix_rst7_oct.in <<EOF 
parm temp_${mol}.parm7
trajin temp_${mol}.rst7
box x ${xval} y ${xval} z ${xval} truncoct
trajout ${mol}.rst7
run 
EOF

		cpptraj -i fix_parm7_oct.in 
		cpptraj -i fix_rst7_oct.in 
		rm temp_${mol}.parm7
                rm temp_${mol}.rst7

        fi


    done


              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-1_aq.parm7  ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/unisc.parm7
              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-1_aq.rst7   ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateA.rst7
              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_aq.rst7   ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateB.rst7
        if [ ${equil_type} == 1 ]; then
                cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateA.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/equil/0.00000000_init.rst7
                cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateB.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/equil/1.00000000_init.rst7
        elif [ ${equil_type} == 2 ]; then
                for(( t=1;t<=${ntrials};t++));do
                        cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateA.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/equil${t}/0.00000000_init.rst7
                        cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/stateB.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/aq/equil${t}/1.00000000_init.rst7
                done
        fi

done
fi

if [ "${ticalc}" == "rbfe" ]; then                         
  if [ -f "fix_box_com_size.txt" ]
     then
             read -r max_a < fix_box_com_size.txt
	     read -r max_b < <(sed -n 2p fix_box_com_size.txt)
	     read -r max_c < <(sed -n 3p fix_box_com_size.txt)
     else

                        max_a=$(awk '{print $1}' *-1_com*rst7 | tail -n 1 | sort -n | tail -n 1)
                        max_b=$(awk '{print $2}' *-1_com*rst7 | tail -n 1 | sort -n | tail -n 1)
                        max_c=$(awk '{print $3}' *-1_com*rst7 | tail -n 1 | sort -n | tail -n 1)
   			echo $max_a >> fix_box_com_size.txt
                        echo $max_b >> fix_box_com_size.txt
                        echo $max_c >> fix_box_com_size.txt
  fi

for com_file in *-1_com*rst7;do
    # Skip files that are from previous generations.
    if [[ "${com_file}" == "fix_box"* ]]; then
        continue
    fi

     
     FILENAME_WITHOUT_EXT=$(echo "$com_file" | sed 's/-[0-9]\+_com.rst7$//')
     BEFORE_TILDE=$(echo "$FILENAME_WITHOUT_EXT" | cut -d'~' -f1)
     AFTER_TILDE=$(echo "$FILENAME_WITHOUT_EXT" | cut -d'~' -f2)
     

     printf "**********************************************\n"
     printf "com file: $com_file\n"
     printf "Before tilde: $BEFORE_TILDE\n"
     printf "After tilde: $AFTER_TILDE\n"
     
     for num in 1 2;do
     cat <<EOF > fix_box_com_cpptraj.in
     parm ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.parm7
     trajin ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.rst7
     trajout ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.pdb pdb include_ep
     go
     quit
EOF
     
     cpptraj -i fix_box_com_cpptraj.in
     
           # assign protein forcefield
        if [ "${pff}" == "ff14SB" ]; then
                printf "source leaprc.protein.ff14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "source leaprc.phosaa14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadamberparams frcmod.ff14SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        elif [ "${pff}" == "ff19SB" ]; then
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "source leaprc.phosaa19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadamberparams frcmod.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        elif [ "${pff}" == "nucleic" ]; then   
                printf "source leaprc.protein.ff19SB\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "source leaprc.RNA.OL3\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "source leaprc.DNA.OL21\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        fi
	

        # assign ligand forcefield
        if [ "${lff}" == "gaff2" ]; then
                printf "source leaprc.gaff2\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        elif [ "${lff}" == "gaff" ]; then
                printf "source leaprc.gaff\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        fi

        # assign water model
        if [ "${wm}" == "tip4pew" ]; then
                printf "source leaprc.water.tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadamberparams frcmod.tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadAmberParams frcmod.ionsjc_tip4pew\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadoff tip4pewbox.off\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                boxkey="TIP4PEWBOX"
        elif [ "${wm}" == "tip3p" ]; then
                printf "source leaprc.water.tip3p\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                boxkey="TIP3PBOX"
        elif [ "${wm}" == "opc" ]; then
                printf "source leaprc.water.opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadamberparams frcmod.opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadamberparams frcmod.ionslm_hfe_opc\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                printf "loadoff opcbox.off\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                boxkey="OPCBOX"
        fi

        # check and load non-standard residue parameter files
        #i=0
#	numnonstd=$(($numnonstd+0))
#        while [ "$i" -lt "${numnonstd}" ]; do
#                printf "loadamberparams ${lig1}_${i}.frcmod\n" >> tleap.in
#                printf "loadoff ${lig1}_${i}.lib\n" >> tleap.in
#                i=$(($i+1))
#        done

        # load ligand 2 parameter files
#        printf "loadamberparams ${lig2}_0.frcmod\n" >> tleap.in
#        printf "loadoff ${lig2}_0.lib\n" >> tleap.in

        printf "loadamberparams ${BEFORE_TILDE}_0.frcmod\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        printf "loadoff ${BEFORE_TILDE}_0.lib\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        printf "loadamberparams ${AFTER_TILDE}_0.frcmod\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        printf "loadoff ${AFTER_TILDE}_0.lib\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in

        # assign MD box
        if [ "${mdboxshape}" == "cubic" ]; then
                boxcmd="solvateBox"
        elif [ "${mdboxshape}" == "oct" ]; then
                boxcmd="solvateOct"
        fi

        # load pdb, pdb with sequence, or mol2
#        if [ "${load}" == "pdb" ]; then
                printf "x = loadPdb ${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.pdb\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
#        elif [ "${load}" == "pdbseq" ]; then
#                printf "x = loadPdbUsingSeq ${inpfile}.pdb { $(cat ${inpfile}.seq) }\n" >> tleap.in
#        else
#                printf "x = loadmol2  ${inpfile}_0.mol2\n" >> tleap.in
#        fi

        # add S-S cysteine linkkages if present
        if [ -f ${inpfile}_sslinks ] && [ "$(cat ${inpfile}_sslinks | wc -l)" -gt 0 ]; then
                while read line; do
                        IFS=' ' read -ra args <<< $line
                        printf "bond x.${args[0]}.SG x.${args[1]}.SG\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
                done < ${inpfile}_sslinks
        fi

        # build box and neutralize with Na+ Cl-
#        if [ "${boxbuild}" == 0 ] && [ "${s}" == "com" ]; then
                printf "setbox x vdw \n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
#        else
        printf "set x box {$max_a $max_b $max_c}\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        printf "saveamberparm x fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.parm7 fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.rst7\n\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        printf "quit\n" >> fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in
        
        tleap -s -f fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com_tleap.in >> fix_box_com_log

        # Reperform HMR if needed
        if [ "${hmr}" == "true" ]; then
                printf "Applying HMR on ${BEFORE_TILDE}~${AFTER_TILDE}\n"
                if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                        cat <<EOFM > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOFM
                parmed -i hmr.in -p fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.parm7 -c fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.rst7 >> output 2>&1
                mv hmr.parm7 fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.parm7; mv hmr.rst7  fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com.rst7
        fi

        # Redo fix box if mdboxshape is oct 
        if [ "${mdboxshape}" == "oct" ]; then
                mol="fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-${num}_com"
                boxinfo=$(tail -n 1 ${mol}.rst7)
                xval=$(echo "$boxinfo" | awk '{print $1}')
                mv ${mol}.parm7 temp_${mol}.parm7
                mv ${mol}.rst7 temp_${mol}.rst7
                cat > fix_parm7_oct.in <<EOF 
parm temp_${mol}.parm7
parmbox truncoct x ${xval} 
parmwrite out ${mol}.parm7
run 
EOF

                cat > fix_rst7_oct.in <<EOF 
parm temp_${mol}.parm7
trajin temp_${mol}.rst7
box x ${xval} y ${xval} z ${xval} truncoct
trajout ${mol}.rst7
run 
EOF

                cpptraj -i fix_parm7_oct.in
                cpptraj -i fix_rst7_oct.in
                rm temp_${mol}.parm7
                rm temp_${mol}.rst7

        fi

  done
              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-1_com.parm7  ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/unisc.parm7
              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-1_com.rst7   ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateA.rst7
              cp fix_box_${BEFORE_TILDE}~${AFTER_TILDE}-2_com.rst7   ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateB.rst7
        if [ ${equil_type} == 1 ]; then
                cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateA.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/equil/0.00000000_init.rst7
                cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateB.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/equil/1.00000000_init.rst7
        elif [ ${equil_type} == 2 ]; then
                for(( t=1;t<=${ntrials};t++));do
                        cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateA.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/equil${t}/0.00000000_init.rst7
                        cp ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/stateB.rst7 ${path}/${system}/${protocol}/run/${BEFORE_TILDE}~${AFTER_TILDE}/com/equil${t}/1.00000000_init.rst7
                done
        fi
done
fi

cd ${path}
echo "Running: Image_Writer on ${system}"
python3 ${pathTOWFToolKit}/bin/fewf-image_writer.py --sys ${system} --image_dir results_${system}/imgdir --sub_dir aq --showidxs --showlabels
echo "Finished Image_Writer"

echo "Copying initial files to results directory."
mkdir -p results_${system}/inputs
cp -r ${path_to_input} results_${system}/inputs/
cp input results_${system}/inputs/
if [[ -v override_lambda ]]; then 
        cp -r ${override_lambda} results_${system}/inputs/
fi
if [[ -v read_map ]]; then 
        cp -r ${read_map} results_${system}/inputs/
fi
echo "Done copying initial files."

# END of setupmode=0


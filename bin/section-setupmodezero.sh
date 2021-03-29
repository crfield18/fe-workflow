        if [ "${setupmode}" == 0 ]; then
                # setup FE simulations using 1-state model
                ##########################
                ##########################
                ##########################

                #UNIFIED protocol
                if [ "${protocol}" == "unified" ]; then
                        if [ "${twostate}" == "false" ]; then
                                cd $path/$system/setup
                                        for i in "${!translist[@]}";do
                                                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                                lig=$(grep -w ${stA} molname-ligname.mapping |awk -F " " '{print $2}')
                                                for s in ${slist[@]}; do
                                                        parmutils-timutate.py -p ${stA}_${s}.parm7 -c ${stA}_${s}.rst7 --target ":${lig}" --mol2 ${stB}_lig_dry.mol2 --lib ${stB}_lig_dry.lib --frcmod ${stB}_lig_dry.frcmod --uniti --nlambda ${nlambda} --map ${translist[$i]}.map.txt >> output
                                                        mkdir -p ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}
                                                        for ticopyfile in ticopy.*; do
                                                                if [ -e "${ticopyfile}" ]; then mv ${ticopyfile} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s} ;fi
                                                        done
                                                        for fffiles in ${stA}_lig_dry.mol2 ${stA}_lig_dry.lib ${stA}_lig_dry.frcmod ${stB}_lig_dry.mol2 ${stB}_lig_dry.lib ${stB}_lig_dry.frcmod; do
                                                                if [ -e "${fffiles}" ]; then cp ${fffiles} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/; fi
                                                        done

                                                        cd ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}
                                                                sh ticopy.sh >> output 2>&1
                                                                sleep 1
                                                                timask1=$(grep molmask2 ticopy.sh|awk -F " " '{print $2}'|awk -F "=" '{print $2}')
                                                                scmask1=$(grep molmask2 ticopy.sh|awk -F " " '{print $3}'|awk -F "=" '{print $2}')
                                                                timask2=$(grep molmask2 ticopy.sh|awk -F " " '{print $5}'|awk -F "=" '{print $2}')
                                                                scmask2=$(grep molmask2 ticopy.sh|awk -F " " '{print $6}'|awk -F "=" '{print $2}')
                                                                noshakemask="${timask1:0:-1},${timask2:2:${#timask2}}"
                                                                if [ "${ticalc}" == "rbfe" ]; then
                                                                        writetemplate_rbfe $nlambda $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB $s ${twostate}
                                                                else
                                                                        writetemplate_rsfe $nlambda $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB ${twostate}
                                                                fi
                                                                ##########################################
                                                                # setup H-mass repartitioning
                                                                if [ "${hmr}" == "true" ]; then
                                                                        if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                                                                        cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                                                                        parmed -i hmr.in -p unisc/template/unisc.parm7 -c unisc/template/unisc.rst7 >> output 2>&1
                                                                        sleep 1
                                                                        mv hmr.parm7 unisc/template/unisc.parm7
                                                                        mv hmr.rst7  unisc/template/unisc.rst7
                                                                        rm -rf hmr.in
                                                                fi
                                                                ##########################################
                                                        cd $path/${system}/setup
                                                done
                                        done
                                cd ${path}
                        fi

                        # setup FE simulations using 2-state model
                        if [ "${twostate}" == "true" ]; then
                                cd $path/$system/setup
                                        for i in "${!translist[@]}";do
                                                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                                ligA=$(grep -w ${stA} molname-ligname.mapping |awk -F " " '{print $2}'); ligB=$(grep -w ${stB} molname-ligname.mapping |awk -F " " '{print $2}')
                                                awk ' { t = $1; $1 = $3; $3 = t; print; } ' ${stA}~${stB}.map.txt | column -t > ${stB}~${stA}.map.txt

                                                for s in ${slist[@]}; do
                                                        parmutils-timutate.py -p ${stA}_${s}.parm7 -c ${stA}_${s}.rst7 --target ":${ligA}" --mol2 ${stB}_lig_dry.mol2 --lib ${stB}_lig_dry.lib --frcmod ${stB}_lig_dry.frcmod --uniti --nlambda ${nlambda} --map ${stA}~${stB}.map.txt >> output 2>&1
                                                        mkdir -p        ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stA}-${stB}
                                                        for ticopyfile in ticopy.*; do
                                                                if [ -e "${ticopyfile}" ]; then mv ${ticopyfile} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stA}-${stB}/${ticopyfile}; fi
                                                        done

                                                        parmutils-timutate.py -p ${stB}_${s}.parm7 -c ${stB}_${s}.rst7 --target ":${ligB}" --mol2 ${stA}_lig_dry.mol2 --lib ${stA}_lig_dry.lib --frcmod ${stA}_lig_dry.frcmod --uniti --nlambda ${nlambda} --map ${stB}~${stA}.map.txt >> output 2>&1
                                                        mkdir -p        ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}
                                                        for ticopyfile in ticopy.*; do
                                                                if [ -e "${ticopyfile}" ]; then mv ${ticopyfile} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/${ticopyfile}; fi
                                                        done


                                                        for fffiles in ${stA}_lig_dry.mol2 ${stA}_lig_dry.lib ${stA}_lig_dry.frcmod ${stB}_lig_dry.mol2 ${stB}_lig_dry.lib ${stB}_lig_dry.frcmod; do
                                                                if [ -e "${fffiles}" ]; then
                                                                        cp ${fffiles} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stA}-${stB}/
                                                                        cp ${fffiles} ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/
                                                                fi
                                                        done

                                                        sed -i  -e "s/mol0 = loadPdbUsingSeq ticopy.mol0.pdb { ${ligB} }/mol0 = loadPdbUsingSeq ticopy.mol0.pdb { ${ligA} }/g" \
                                                                -e "s/mol1 = loadPdbUsingSeq ticopy.mol1.pdb { ${ligA} }/mol1 = loadPdbUsingSeq ticopy.mol1.pdb { ${ligB} }/g" \
                                                                ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.sh

                                                        masks=$(grep molmask ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stA}-${stB}/ticopy.sh)
                                                        sed -i "/molmask/c\\${masks}\\" ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.sh

                                                        mv  ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.mol0.pdb tmp.pdb
                                                        mv  ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.mol1.pdb ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.mol0.pdb
                                                        mv  tmp.pdb ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}/ticopy.mol1.pdb

                                                        cd ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stA}-${stB}
                                                                sh ticopy.sh >> output 2>&1; sleep 1
                                                                cp unisc/template/unisc.parm7 ../merged-lam0_${s}.parm7
                                                                cp unisc/template/unisc.rst7  ../merged-lam0_${s}.rst7
                                                        cd ${path}
                                                        cd ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/${stB}-${stA}
                                                                sh ticopy.sh >> output 2>&1; sleep 1
                                                                cp unisc/template/unisc.parm7 ../merged-lam1_${s}.parm7
                                                                cp unisc/template/unisc.rst7  ../merged-lam1_${s}.rst7
                                                        cd ${path}

                                                        if [ "$s" == "com" ]; then rbuf=${boxbufcom}; else rbuf=${boxbufaq}; fi
                                                        cd ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/
                                                                mergeligs=(merged-lam0_${s} merged-lam1_${s})

                                                                parmutils-Organize.py -p ${mergeligs[0]}.parm7 -c ${mergeligs[0]}.rst7 -e1 '!:WAT,K+,Na+,Cl-' -o1 "${mergeligs[0]}"
                                                                parmutils-Organize.py -p ${mergeligs[1]}.parm7 -c ${mergeligs[1]}.rst7 -e1 '!:WAT,K+,Na+,Cl-' -o1 "${mergeligs[1]}"
                                                                create_box_equal "${s}" "${rbuf}" "${ionconc}" "${mergeligs[@]}"
                                                                timask1=$(grep molmask2 ${stA}-${stB}/ticopy.sh|awk -F " " '{print $2}'|awk -F "=" '{print $2}')
                                                                scmask1=$(grep molmask2 ${stA}-${stB}/ticopy.sh|awk -F " " '{print $3}'|awk -F "=" '{print $2}')
                                                                timask2=$(grep molmask2 ${stA}-${stB}/ticopy.sh|awk -F " " '{print $5}'|awk -F "=" '{print $2}')
                                                                scmask2=$(grep molmask2 ${stA}-${stB}/ticopy.sh|awk -F " " '{print $6}'|awk -F "=" '{print $2}')
                                                                noshakemask="${timask1:0:-1},${timask2:2:${#timask2}}"
                                                                if [ "${ticalc}" == "rbfe" ]; then
                                                                        writetemplate_rbfe $nlambda $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB $s ${twostate}
                                                                else
                                                                        writetemplate_rsfe $nlambda $cutoff $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $noshakemask $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB ${twostate}
                                                                fi

                                                                ##########################################
                                                                for lm in lam0 lam1; do
                                                                        if [ "${hmr}" == "true" ]; then
                                                                                # setup H-mass repartitioning
                                                                                if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                                                                                cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                                                                                parmed -i hmr.in -p merged-${lm}_${s}.parm7 -c merged-${lm}_${s}.rst7 >> output 2>&1
                                                                                sleep 1
                                                                                mv hmr.parm7 merged-${lm}_${s}.parm7
                                                                                mv hmr.rst7  merged-${lm}_${s}.rst7
                                                                                rm -rf hmr.in
                                                                        fi
                                                                        mkdir -p unisc/template/
                                                                        cp ${stA}-${stB}/unisc/template/* unisc/template/
                                                                        cp merged-${lm}_${s}.parm7 unisc/template/merged-${lm}.parm7
                                                                        cp merged-${lm}_${s}.rst7  unisc/template/merged-${lm}.rst7
                                                                done
                                                                ##########################################
                                                        cd ${path}/$system/setup
                                                done
                                        done
                                cd ${path}
                        fi
                fi
                ##########################
                ##########################
                ##########################

                # Setup input file infrastructure
                ##########################
                ##########################
                ##########################
                #UNIFIED protocol
                if [ "${protocol}" == "unified" ]; then
                        for i in "${!translist[@]}";do
                                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                for s in ${slist[@]}; do

                                        cd ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}
                                                sh TEMPLATE.sh; sleep 1
                                        cd ${path}

                                        # if hmr=true set timestep to 4fs
                                        if [ "${hmr}" == "true" ]; then
                                                sed -i '/dt.*.=.*.*/c\dt              = 0.004' ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/inputs/*_ti.mdin
                                        fi

                                        # if repex=false, alter input files and slurm files
                                        if [ "${repex}" == "false" ]; then
                                                sed -i -e '/numexchg/d' -e '/gremd_acyc/d' ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/inputs/*_ti.mdin
                                                sed -i -e 's/ -rem 3//g' -e 's/running replica ti/running regular ti/g' ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/prod.slurm
                                        fi

                                        rm -rf ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/current ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/inputs

                                        mkdir -p ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
                                        cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/
                                                mv -f ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/unisc.parm7        .
                                                mv -f ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/current            .
                                                mv -f ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/inputs             .
                                                mv -f ${path}/${system}/${protocol}/build/${stA}~${stB}/${s}/*slurm             .

                                                truncate -s0 sub-eq.sh sub-ti.sh
                                                for(( t=1;t<=${ntrials};t++));do
                                                        mkdir -p t${t}
                                                        cp current/*_init.rst7 t${t}/
                                                        if [ "${twostate}" == true ]; then
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

                                                        sed     -e "s/current/t${t}/g" \
                                                                -e "s/EOF2/EOF/g" equilgroup.slurm > t${t}_eq.slurm
                                                        sed     -e "s/current/t${t}/g" \
                                                                -e "s/ti.groupfile/t${t}_ti.groupfile/g" prod.slurm > t${t}_ti.slurm

                                                        echo "sbatch t${t}_eq.slurm" >> sub-eq.sh
                                                        echo "sbatch t${t}_ti.slurm" >> sub-ti.sh
                                                done
                                        cd ${path}
                                done
                                echo "Done with ${stA}~${stB}..."
                        done
                fi

                #rm -rf ${path}/*mol2 ${path}/genmol2.in ${path}/*lib  ${path}/*frcmod ${path}/*parm7 ${path}/*rst7 ${path}/*pdb ${path}/map.txt
        fi
        # END of setupmode=0


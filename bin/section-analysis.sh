
# analysis section 
printf "\n\n Running in analysis mode... \n\n" 
mkdir -p results/data
write_analyze
for i in "${translist[@]}";do
	for s in ${slist[@]};do
		if [ "${s}" == "com" ]; then bidirection=${bidirection_com}; else bidirection=${bidirection_aq}; fi

		if [ -d "results/data/${i}/${s}" ]; then
			printf "\n Folder results/data/${i}/${s} already present. Skipping initial extraction.\n"
			continue
		else
			if [ "${bidirection}" == "false" ]; then
				mkdir -p results/data/${i}/${s}
				for trial in $(seq 1 1 ${ntrials});do
					path_out=results/data/${i}/${s}/${trial}
					path_data=${path_to_data}/${i}/${s}/t${trial}
					if [ -d "${path_data}" ]; then 
						./mdouts2dats.py --odir=${path_out} $(ls ${path_data}/*ti.mdout)
					fi
				done
			else
				mkdir -p results/data/${i}/${s} results/data/${i}/${s}tmp
				for trial in $(seq 1 1 ${ntrials});do
					path_out=results/data/${i}/${s}/${trial}
					path_data=${path_to_data}/${i}/${s}/forward/t${trial}
					if [ -d "${path_data}" ]; then 
						./mdouts2dats.py --odir=${path_out} $(ls ${path_data}/*ti.mdout)
					fi

					path_out=results/data/${i}/${s}tmp/${trial}
					path_data=${path_to_data}/${i}/${s}/reverse/t${trial}
					if [ -d "${path_data}" ]; then
						./mdouts2dats.py --odir=${path_out} $(ls ${path_data}/*ti.mdout)
					fi
				done
				for trial in $(seq 1 1 ${ntrials});do
					fktrial=$((${trial} + ${ntrials}))
					path_tar=results/data/${i}/${s}/${fktrial}
					path_src=results/data/${i}/${s}tmp/${trial}
					mkdir -p ${path_tar}
					lamsreversed=($(for lr in ${lams[@]}; do echo $lr; done | sort -r))

					for li in ${!lams[@]}; do
						lamri=${lams[${li}]}; lamfi=${lamsreversed[${li}]}
						mv ${path_src}/dvdl_${lamri}.dat ${path_tar}/dvdl_${lamfi}.dat
						for lj in ${!lams[@]}; do
							lamrj=${lams[${lj}]}; lamfj=${lamsreversed[${lj}]}
							mv ${path_src}/efep_${lamri}_${lamrj}.dat ${path_tar}/efep_${lamfi}_${lamfj}.dat
						done
					done
				done
				rm -rf results/data/${i}/${s}tmp
			fi
		fi
	done
done
rm -rf mdouts2dats.py


if [ "${exptdatafile}" == "skip" ]; then
        truncate -s0 Expt.dat
        for i in ${uniqueligs[@]};do
                printf "${i}    0 \n" >> Expt.dat
        done
        exptdatafile="Expt.dat"
fi
cp ${exptdatafile} results/


cd results
	write_gmbar ${ticalc}

	#######################################
	export PATH="$PATH:${pathTOFEToolKit}/local/bin"
	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${pathTOFEToolKit}/local/lib"
	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${pathTOFEToolKit}/local/lib64"
	export PYTHONPATH="${pathTOFEToolKit}/local/lib/python3.8/site-packages"
	echo "${PYTHONPATH}"
	source ${pathTOFEToolKit}/local/bin/activate
	#######################################

	if [ "${bidirection}" == "true" ]; then ntrials=$(($ntrials*2)); fi

	gmbarargs=("-f ${exptdatafile}")
	graphmbarargs=(--nboot=1 --fc=10 --guess=2)
	#graphmbarargs=(--nboot=200 --fc=10 --guess=2)

	if [ "${ccc}" == "true" ] && [ "${ccc_ddG}" == "true" ]; then gmbarargs=("--ddG" "${gmbarargs[@]}"); fi
	if [ "${bar}" == "true" ]; then gmbarargs=("--bar" "${gmbarargs[@]}"); fi
	if [ "${showallcycles}" == "true" ]; then gmbarargs=("--super" "${gmbarargs[@]}"); fi
	if [ "${ccc}" == "false" ]; then gmbarinpargs=("--nocyc" "${gmbarargs[@]}"); else gmbarinpargs=("${gmbarargs[@]}"); fi

	#echo "python3 ./gmbar.py $(echo ${gmbarargs[*]}) -t $(seq 1 ${ntrials} | xargs -n ${ntrials} echo) > graphmbar.inp" 
	python3 ./gmbar.py $(echo ${gmbarinpargs[*]}) -t $(seq 1 ${ntrials} | xargs -n ${ntrials} echo) > graphmbar.inp


	if [ "${check_convergence}" == "false" ]; then 
		graphmbar $(echo ${graphmbarargs[*]}) --start ${start} --stop ${stop} graphmbar.inp > graphmbar.out
		python3 ./gmbar.py $(echo ${gmbarargs[*]}) -t $(seq 1 ${ntrials} | xargs -n ${ntrials} echo) -a graphmbar.out > out
	else
		for i in $(seq ${start} 20 ${stop}); do
        		if [ "${start}" != ${i} ]; then
				graphmbar $(echo ${graphmbarargs[*]}) --start ${start} --stop ${i} graphmbar.inp > graphmbar_${start}-${i}.out
				python3 ./gmbar.py $(echo ${gmbarargs[*]}) -t $(seq 1 ${ntrials} | xargs -n ${ntrials} echo) -a graphmbar_${start}-${i}.out > ${start}-${i}.out
        		fi

        		if [ "${i}" != ${stop} ]; then
				graphmbar $(echo ${graphmbarargs[*]}) --start ${i} --stop ${stop} graphmbar.inp > graphmbar_${i}-${stop}.out
				python3 ./gmbar.py $(echo ${gmbarargs[*]}) -t $(seq 1 ${ntrials} | xargs -n ${ntrials} echo) -a graphmbar_${i}-${stop}.out > ${i}-${stop}.out
        		fi

		done


	fi


	conda deactivate

cd $path



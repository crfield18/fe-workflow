
# analysis section 
printf "\n\n Running in analysis mode... \n\n" 
mkdir -p results/data

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
						${pathTOFEToolKit}/edgembar/src/python/bin/edgembar-amber2dats.py -o ${path_out} $(ls ${path_data}/*ti.mdout)
					fi
				done
			else
				mkdir -p results/data/${i}/${s} results/data/${i}/${s}tmp
				for trial in $(seq 1 1 ${ntrials});do
					path_out=results/data/${i}/${s}/${trial}
					path_data=${path_to_data}/${i}/${s}/forward/t${trial}
					if [ -d "${path_data}" ]; then 
						edgembar-amber2dats.py --odir=${path_out} $(ls ${path_data}/*ti.mdout)
					fi

					path_out=results/data/${i}/${s}tmp/${trial}
					path_data=${path_to_data}/${i}/${s}/reverse/t${trial}
					if [ -d "${path_data}" ]; then
						${pathTOFEToolKit}/edgembar/src/python/bin/edgembar-amber2dats.py -o ${path_out} $(ls ${path_data}/*ti.mdout)
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


if [ "${exptdatafile}" == "skip" ]; then
        truncate -s0 Expt.dat
        for i in ${uniqueligs[@]};do
                printf "${i}    0 \n" >> Expt.dat
        done
        exptdatafile="Expt.dat"
fi
cp ${exptdatafile} results/


cd results
        write_edgembar ${ticalc}

if [ "${bidirection}" == "true" ]; then ntrials=$(($ntrials*2)); fi

echo "Running: python3 DiscoverEdges.py"
python3 ./DiscoverEdges.py
echo "Finished DiscoverEdges.py"

for xml in analysis/*.xml; do
    if [ -e "${xml}" ]; then
        echo ""
        echo "Running: time OMP_NUM_THREADS=4 edgembar.OMP --halves --fwdrev ${xml}"
        time OMP_NUM_THREADS=4 edgembar.OMP --halves --fwdrev ${xml}
        echo "Finished creating ${xml%.xml}.py"
    fi
done

for xml in analysis/*.xml; do
    if [ "${ticalc}" == "asfe" ]; then
        cp ${xml%.xml}.py ${xml%.xml}~.py
        py=${xml%.xml}~.py
    else
        py=${xml%.xml}.py
    fi
    if [ -e "${py}" ]; then
        echo ""
        echo "Running: python3 ${py}"
        python3 ${py}
        echo "Finished creating ${xml%.xml}.html"
    fi
done


echo ""
echo "Running: edgembar-WriteGraphHtml.py -o analysis/Graph.html \$(ls analysis/*~*.py)"
edgembar-WriteGraphHtml.py -o analysis/Graph.html $(ls analysis/*~*.py)
echo "Finished creating analysis/Graph.html"

if  [ "${ticalc}" != "asfe" ]; then
 
echo ""
echo "Running: edgembar-WriteGraphHtml.py -o analysis/GraphWithExpt.html -x ExptVals.txt \$(ls analysis/*~*.py)"
edgembar-WriteGraphHtml.py -o analysis/GraphWithExpt.html -x ${exptdatafile} $(ls analysis/*~*.py)
echo "Finished creating analysis/GraphWithExpt.html"
fi

cd $path




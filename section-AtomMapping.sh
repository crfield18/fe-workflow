# generation of atom-mapping between ligands
##########################
##########################
##########################
if [ "${ticalc}" != "asfe" ]; then

	if [ "${mapinspect}" -ne 2 ]; then
        	cd ${system}/setup
                	rm -rf *.map.txt
                	if [ "${mapnetwork}" == true ]; then
                        	for i in "${!translist[@]}";do
                                	l1=$(basename ${translist[$i]}); l2="${l1##*~}"; l1="${l1%~*}"
                                	listl1+=("${l1}"); listl2+=("${l2}")
                        	done
                        	l1uniques=($(echo "${listl1[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
                        	truncate -s0 map-network
                        	for l1 in "${!l1uniques[@]}";do
                                	echo -n "${l1uniques[$l1]}_0.mol2" >> map-network
                                	for j in "${!translist[@]}";do
                                        	l11=$(basename ${translist[$j]}); l22="${l11##*~}"; l11="${l11%~*}"
                                        	if [ "${l1uniques[$l1]}" == "${l11}" ]; then
                                                	echo -n " ${l22}_0.mol2 " >> map-network
                                        	fi
                                	done
                                	echo "" >> map-network
                        	done
                        	cat map-network |column -t > tmp && mv tmp map-network
                        	parmutils-scmapper.py --graph map-network -t ${mapmethod} >> output 2>&1

                        	for map in *map.txt; do
                                	mv -f ${map} $(echo ${map}|awk -F "_0"  '{print $1$2$3}')
                        	done

                	else
                        	for i in ${!listA[@]};do
                                	l1=${listA[$i]}; l2=${listB[$i]}
                               		parmutils-scmapper.py -a ${l1}_0.mol2 -b ${l2}_0.mol2 -o ${l1}~${l2}.map.txt -t ${mapmethod} >> output 2>&1
                        	done

                	fi
        	cd ${path}
	fi

	if [ "${mapinspect}" -eq 1 ]; then
        	printf "\n\ncheck *.map.txt files in ${system}/setup. Once done, proceed by re-running script with mapinspect = 2\n\n" && exit 0
	fi

	if [ "${mapinspect}" -eq 2 ];then
        	cd ${system}/setup
                	mapmissing=0
                	for i in "${!translist[@]}";do
                        	if [ ! -f ${translist[$i]}.map.txt ]; then printf "\n\n${translist[$i]}.map.txt missing in working directory.\n\n" && mapmissing=1; fi
                	done
                	if [ "${mapmissing}" -eq 1 ]; then exit 0; fi
        	cd ${path}
	fi
else

	printf "\n\n\"ticalc\" is set to \"asfe\". Skipping atom mapping.\n\n"

fi
##########################
##########################
##########################


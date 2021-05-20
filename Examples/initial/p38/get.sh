#!/bin/bash

srcpath="/run/media/abir/hdd1/work/afe/stx/auto/work-in-progress/initial/p38/com"
path=`pwd`

for i in ${srcpath}/*parm7; do
	f=$(basename $i); f="${f/_com.parm7/}"
	python ~/bin/parmutils/parmutils-Organize.py -p ${srcpath}/${f}_com.parm7 -c ${srcpath}/${f}_com.rst7 -m :${f} -n LIG -e1 :LIG -o1 ${f}_0 -e2 '!:WAT,Na+,Cl-' -o2 ${f}
	rm None.* ${f}.mol2 ${f}.lib ${f}.seq ${f}.frcmod
	echo "Done with ${f}..."
done
#./parmutils-Organize.py -p 2AA_com.parm7 -c 2AA_com.rst7 -m :2AA -n LIG -e1 :LIG -f 2AA_0 -l 2AA_0 -o1 2AA_0 -e2 '!:WAT,Na+,Cl-' -o2 2AA

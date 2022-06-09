#!/bin/bash
#SBATCH --job-name="eq_methanol.slurm"
#SBATCH --output="eq_methanol.slurm.slurmout"
#SBATCH --error="eq_methanol.slurm.slurmerr"
#SBATCH --partition=general-long-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --gres=gpu:8
#SBATCH --time=1-00:00:00

top=${PWD}
endstates=(0.00000000)
lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)
eqstage=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA minTI eqpre1P0TI eqpre2P0TI eqP0TI eqATI preTI)


# check if AMBERHOME is set
if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi
# check if cpptraj is present
if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi

EXE=${AMBERHOME}/bin/pmemd.cuda

for trial in 4; do
	
	if [ ! -d t${trial} ];then mkdir t${trial}; fi

	count=-1; alllams=0
	for stage in ${eqstage[@]}; do
        	count=$(($count+1))
        	lastcount=$(($count-1))
		if [ "${stage}" == "init" ] || [ "${stage}" == "eqpre1P0TI" ] || [ "${stage}" == "eqpre2P0TI" ] || [ "${stage}" == "eqP0TI" ]; then continue; fi
        	laststage=${eqstage[${lastcount}]}

		if [ "${stage}" == "minTI" ];then alllams=1; fi

		
        	if [ ${alllams} -eq 0 ];then


			# check if pmemd.cuda is present
                        if ! command -v ${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
			export LAUNCH="srun"
			export EXE=${AMBERHOME}/bin/pmemd.cuda

                	lam=0.00000000
                	echo "Running $stage for lambda ${lam}..."
                	${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
                	cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                	# check if cpptraj is present
                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                	cpptraj < center.in
                	sleep 1
                	mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7

		elif [ ${alllams} -eq 1 ] && [ "${stage}" == "minTI" ];then
			# check if pmemd.cuda is present
			if ! command -v ${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
			export LAUNCH="srun"
			export EXE=${AMBERHOME}/bin/pmemd.cuda
			for i in ${!lams[@]}; do
				lam=${lams[$i]}
				if [ "${i}" -eq 0 ]; then
					init=${endstates[0]}_eqA.rst7
 				else
					init=${lams[$(($i-1))]}_eqP0TI.rst7
				fi			

				echo "Running $stage for lambda ${lam}..."

				stage=minTI
                                ${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${init} -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${init}
                                sleep 1

                                laststage=minTI; stage=eqpre1P0TI
                                ${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
                                sleep 1

                                laststage=eqpre1P0TI; stage=eqpre2P0TI
                                ${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
                                sleep 1

                                laststage=eqpre2P0TI; stage=eqP0TI
                                ${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
                                sleep 1

                                cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                cpptraj < center.in
                                sleep 1
                                mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7
                        done
                        laststage=eqP0TI
		else
                        # check if pmemd.cuda.MPI is present
                        if ! command -v ${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                        export LAUNCH="mpirun -np ${#lams[@]}"
                        export EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
                        export MV2_ENABLE_AFFINITY=0
                        ${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/t${trial}_${stage}.groupfile

			for lam in ${lams[@]};do
				cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                cpptraj < center.in
                                sleep 1
                                mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7
                        done
		fi
	done

        # run production
        EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
        echo "running replica ti"
        mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/t${trial}_ti.groupfile


	cat << EOFP > extract.py
#!/usr/bin/env python3

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
        if abs(param.box[3]-109.471219)<1.e-4 and            abs(param.box[4]-109.471219)<1.e-4 and            abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param

def CopyParm( parm ):
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p

def Strip( parm, mask ):
    p = CopyParm( parm )
    p.strip( "%s"%(mask) )
    return p

def Extract( parm, mask ):
    return Strip( parm, "!(%s)"%(mask) )

def SaveParmRst( param, fname ):
    from parmed.constants import IFBOX
    for a in param.atoms:
        param.parm_data["CHARGE"][ a.idx ] = a.charge
    if param.box is not None:
       if abs(param.box[3]-109.471219)<1.e-4 and           abs(param.box[4]-109.471219)<1.e-4 and           abs(param.box[5]-109.471219)<1.e-4:
           param.parm_data["POINTERS"][IFBOX]=2
           param.pointers["IFBOX"]=2
    try:
        param.save( "{}.parm7".format(fname), overwrite=True )
        #param.save( fname, overwrite=True )
    except:
        param.save( "{}.parm7".format(fname) )
        #param.save( fname )
    rst = parmed.amber.Rst7(natom=len(param.atoms),title="BLAH")
    rst.coordinates = param.coordinates
    rst.box = [param.box[0], param.box[1], param.box[2], param.box[3], param.box[4], param.box[5]]
    rst.write( "{}.rst7".format(fname) )

if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser     ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Extracts Amber parameters from a parm7 file" )

    parser.add_argument("-p","--parm",
                        help="Amber parm7 file",
                        type=str,
                        required=True)

    parser.add_argument("-c","--crd",
                        help="Amber rst7 file",
                        type=str,
                        required=True)

    parser.add_argument("-m","--mask",
                        help="Amber selection mask. Only residues within this mask will be extracted",
                        type=str,
                        default="@*",
                        required=False )

    parser.add_argument("-o","--output",
                        help="Output file basename if -n/--name is specified",
                        type=str,
                        required=False)

    args = parser.parse_args()

    p    = OpenParm(args.parm,xyz=args.crd)

    if args.mask is not None:
        if len(args.mask) > 0:
            q    = Extract(p,args.mask)
            SaveParmRst(q, args.output)

EOFP
	chmod a+x extract.py

	mkdir -p ../vac/t${trial} ../vac/inputs
	for lam in ${lams[@]};do
        	./extract.py -p unisc.parm7 -c t${trial}/${lam}_preTI.rst7 -m '!:WAT,Na+,K+,Cl-' -o ../vac/t${trial}/${lam}_init
        	sed -e 's/nstlim.*/nstlim          = 500000/g' -e 's/restraint_wt.*/restraint_wt    = 0/g' -e 's/nmropt.*/nmropt          = 0/g' -e '59,65d' -e "s/clambda.*/clambda         = /g" -e 's/irest.*/irest           = 0/g' -e 's/ntx.*/ntx             = 1/g' inputs/0.00000000_eqA.mdin > ../vac/inputs/${lam}_preTI.mdin
        	sed -e 's/ntb.*/ntb             = 1/g' -e '/barostat.*/d' -e '/ntp.*/d' -e '/pres0.*/d' -e '/taup.*/d' inputs/${lam}_ti.mdin > ../vac/inputs/${lam}_ti.mdin
                sed -i -e "s/ntwr.*.=*.*/ntwr            = 5000\nntpr            = 1000/" ../vac/inputs/${lam}_ti.mdin
	done

	cd ../vac
        	cd t${trial}
                	mv 0.00000000_init.parm7 ../unisc.parm7
                        for lam in ${lams[@]}; do
                            cp 0.00000000_init.rst7 ${lam}_init.rst7
                        done
                	rm *.parm7
        	cd ../

        	stage=preTI; laststage=init
        	EXE=${AMBERHOME}/bin/pmemd.cuda
        	LAUNCH="srun --exclusive -N 1 -n 1 -c 1 --gres=gpu:1"
        	for lam in ${lams[@]};do
                	${LAUNCH} ${EXE} -O -p unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7 &
        	done
        	wait
	cd ../

	for dir in vac; do
		cd ${dir}
			stage=ti; laststage=preTI
                        EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
                        echo "running replica ti"
                        mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/t${trial}_ti.groupfile
###			for lam in ${lams[@]};do
###				${LAUNCH} ${EXE} -O -p unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7 &
###			done
		cd ../
	done
	wait

done


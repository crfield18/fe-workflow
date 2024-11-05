##########################################
# read input file
function read_input {

while read line; do
        varlist=(system translist path_to_input nlambda lamschedule lams protocol mapmethod mapinspect mapnetwork boxbuild boxbufcom boxbufaq ionconc pff lff wm mdboxshape ntrials cutoff repex nstlimti numexchgti hmr notrajectory scalpha scbeta gti_add_sc gti_scale_beta gti_cut gti_cut_sc_on gti_cut_sc_off gti_lam_sch gti_ele_sc gti_vdw_sc gti_cut_sc gti_ele_exp gti_vdw_exp stage setupmode twostate bidirection_aq bidirection_com ticalc partition nnodes ngpus wallclock path_to_data exptdatafile bar ccc start stop check_convergence lead_ligand) 
        IFS=$'\t| |=' read -ra args <<< $line
        if [[ "${args[0]}" =~ ^#.* ]]; then continue; fi
        keyword=${args[0]}; value=${args[1]}
        if [ "${args[0]}" == "translist" ] || [ "${args[0]}" == "lams" ]; then
                IFS='[=()]' read -ra args <<< $line
                keyword=${args[0]}; value=${args[2]}
        fi
        for var in ${varlist[@]}; do
                if [ "$var" == "$keyword" ]; then
                        if [ "$var" == "translist" ] || [ "$var" == "lams" ]; then
                                eval "$var=($value)"
                        else
                                eval "$var=$value"
                        fi
                fi
        done
done < $1
}
##########################################


#!/bin/bash

function write_analyze {

	cat << EOF > mdouts2dats.py
#!/usr/bin/env python3
import sys,os
import re

import argparse




def extract_traditional_ti( fname, write=False, odir="" ):
    import os
    from collections import defaultdict as ddict

    fh = open(fname,"r")
    if not fh:
        raise Exception("Could not open %s\n"%(fname))

    numexchg=0
    nstlim=None
    ntpr=None
    dt=None
    irest=0
    mbar_states=0
    mbar_lambda=[]
    lam_values=[]
    for line in fh:
        cmdstr,sepstr,comstr = line.partition("!")
        cmdstr = cmdstr.strip()
        if "mbar_lambda" in cmdstr:
            lams=[]
            cols = cmdstr.replace("="," ").replace(","," ").strip().split()
            for icol in range(len(cols)-1):
                if "mbar_lambda" in cols[icol]:
                    fcol=icol+1
                    break
            cs = [ float(x) for x in cols[fcol:] ]
            m = re.search(r"mbar_lambda\( *([0-9]+) *: *([0-9]+) *\).*",cmdstr)
            if m:
                i0 = int( m.group(1) )-1
                i1 = int( m.group(2) )-1
                while len(mbar_states) < i1+1:
                    mbar_lambda.append(-1.)
                for ii,i in enumerate(range(i0,i1+1)):
                    mbar_lambda[i] = lams[ii]
            else:
                #m = re.search(r"mbar_lambda\( *(%i) *\).*",cmdstr)
                m = re.search(r"mbar_lambda\( *([0-9]+) *\).*",cmdstr)
                if m:
                    i0 = int( m.group(1) )-1
                    while len(mbar_lambda) < i0+1:
                        mbar_lambda.append(-1.)
                    mbar_lambda[i0] = cs[0]
                else:
                    mbar_lambda=cs
        if "mbar_states" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "mbar_states":
                    mbar_states = int( cols[icol+1] )
                    break
        if "lambda values considered:" in line:
            while True:
                line = next(fh)
                if "Extra" in line:
                    break
                cmdstr,sepstr,comstr = line.partition("!")
                if "total:" in cmdstr:
                    cs = cmdstr.split()[2:]
                else:
                    cs = cmdstr.split()
                lam_values.extend( [float(x) for x in cs] )

        if "ntpr" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            #print(cols)
            for icol in range(len(cols)-1):
                if cols[icol] == "ntpr":
                    ntpr = int( cols[icol+1] )
                    break
        if "dt" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "dt":
                    dt = float( cols[icol+1] )
                    break
        if "numexchg" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "numexchg":
                    numexchg = int( cols[icol+1] )
                    break
        if "nstlim" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "nstlim":
                    nstlim = int( cols[icol+1] )
                    break
        if "irest" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "irest":
                    irest = int( cols[icol+1] )
                    break

    if ntpr is None:
        raise Exception("Could not determine ntpr from %s"%(fname))

    if dt is None:
        raise Exception("Could not determine dt from %s"%(fname))

    if nstlim is None:
        raise Exception("Could not determine nstlim from %s"%(fname))

    if numexchg < 1:
        numexchg = 1

    #print(len(mbar_lambda),mbar_states,len(lam_values))

    if len(mbar_lambda) > 0 or mbar_states > 0:
        if len(mbar_lambda) != mbar_states:
            if len(lam_values) == mbar_states:
                mbar_lambda = lam_values
            else:
                print("len(mbar_lambda) != mbar_states: %i vs %i"%(len(mbar_lambda),mbar_states))


    #print(mbar_lambda)


    dt = dt
    nstep_per_sim = nstlim * numexchg
    nframe_per_sim = nstep_per_sim / ntpr

    if nstep_per_sim % ntpr != 0:
        print("num md steps per simulation is not a multiple of ntpr. Unclear how the simulation time works")

    t_per_frame = dt * ntpr
    t_per_sim = t_per_frame * nframe_per_sim


    fh = open(fname,"r")


    efeps = []
    dvdls = []
    efep = []
    reading_region_1 = False

    lam = None
    nlam = 0

    for line in fh:
        if "A V E R A G E S" in line:
            break
        if "MBAR Energy analysis:" in line:
            efep = []
        if "clambda" in line:
            if lam is None:
                cols = line.replace("="," ").replace(","," ").split()
                for i in range(len(cols)):
                    if cols[i] == "clambda":
                        lam = float(cols[i+1])
                        break
        elif "Energy at " in line:
            #print line
            val = line.strip().split()[-1]
            if "****" in val:
                val = 10000.00
                #if len(efep) > 0:
                #   if efep[-1] < 0:
                #       val = -val
            else:
                val = float(val)
            efep.append( val )
        elif "TI region  1" in line:
            reading_region_1 = True
            dvdl = 0
        elif "| TI region  2" in line:
            #print line
            reading_region_1 = False
            #print dvdl
            dvdls.append( dvdl )
            if len( efep ) > 0:
                efeps.append( efep )
                nlam = len(efep)
        elif "TI region " in line:
            reading_region_1 = False

        if "DV/DL  =" in line and reading_region_1:
            #print line
            cols = line.strip().split()
            #print cols
            dvdl = float( cols[-1] )
            #dvdls.append( float( cols[-1] ) )
            #if len( efep ) > 0:
            #    efeps.append( efep )
            #    nlam = len(efep)
    if write:
        #lams = [ float(i) / ( nlam-1. ) for i in range(nlam) ]
        lams = mbar_lambda
        for l in lams:
            if abs(l-lam) < 0.001:
                lam = l
                break
        head, tail = os.path.split(fname)
        if len(odir) > 0:
            head = odir
            if not os.path.isdir(head):
                os.makedirs(head)
        dvdl_fname = os.path.join( head, "dvdl_%.8f.dat"%( lam ) )

        if irest == 0:
           dvdls=dvdls[1:]

        fh = open(dvdl_fname,"w")
        for i in range(len(dvdls)):
            fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,dvdls[i]))
        fh.close()
        for ilam,plam in enumerate(lams):
            efep_fname = os.path.join( head, "efep_%.8f_%.8f.dat"%( lam, plam ) )
            fh = open(efep_fname,"w")
            for i in range(len(efeps)):
                fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,efeps[i][ilam]))
            fh.close()

    return dvdls,efeps



parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""
Extracts DVDL and MBAR data from 1-or-more mdout files and writes
the data into timeseries files
""" )

parser.add_argument \
    ("-o","--odir",
     help="Output directory. If empty, then same dir as mdout",
     type=str,
     required=False,
     default="" )

parser.add_argument \
    ('mdout',
     metavar='mdout',
     type=str,
     nargs='+',
     help='Amber mdout file')

args = parser.parse_args()


odir = None
if len(args.odir) > 0:
    odir = args.odir
for arg in args.mdout:
    if os.path.isfile( arg ):
        if ".mdout" in arg or ".out" in arg:
            extract_traditional_ti( arg, write=True, odir=odir )
        else:
            print("File does not end in .mdout nor .out: %s"%(arg))
    else:
        print("File not found: %s"%(arg))

EOF
	chmod a+x mdouts2dats.py
	
}

function gen_lambdas {
        cat <<EOF > gen_lambda.py
#!/usr/bin/env python3

import numpy as np
import sys

nlam = sys.argv[1]; nlam = float (nlam)
a = np.linspace(0, 1, num = nlam, dtype = float)
for x in range(len(a)):
    print("{:.8f}".format(a[x])),

EOF

local nlambda=$1
local lams=$(python2.7 gen_lambda.py ${nlambda})
echo $lams
rm -rf gen_lambda.py
}


function write_gmbar {

	#nlambda=$1; shift
	ticalc=$1; shift


        cat << EOF > gmbar.py
#!/usr/bin/env python3

import graphmbar as gmbar


class MySetup(gmbar.Analyzer):

    def __init__(self,args):
        super().__init__(args.expt,args.archive,args.bar,args.trials)
        self.SetRefMol(args.lead)
        self.SetExptEdgeConstraints(args.exptcon)


    # ######################################################################
    #
    # These routines are specific to the directory structure of this project
    #
    # ######################################################################


    def getdir( self, transform, env, stage, trial ):
        """
        transform -- string, name of the A->B transformation
        env -- string, the environment. This is either: "solvated" or "complex"
        stage -- string, the transformation step; e.g., "recharge", 
                 "decharge", "vdw"
        trial -- string, the independent trial
        
        returns a directory name where the efep_*_*.dat files are located
        """
EOF
	if [ "${ticalc}" == "rbfe" ]; then
		cat << EOF >> gmbar.py
        if env == "complex":
            myenv = "com"
        elif env == "solvated":
            myenv = "aq"
        return "data/%s/%s/%i/"%(transform,myenv,int(trial))
EOF
	elif [ "${ticalc}" == "rsfe" ]; then
		cat << EOF >> gmbar.py
        if env == "complex":
            myenv = "aq"
        elif env == "solvated":
            myenv = "vac"
        return "data/%s/%s/%i/"%(transform,myenv,int(trial))
EOF

        fi


	cat << EOF >> gmbar.py

    def getedges( self, arcfiles ):
        """
        arcfiles -- dict, keys are filenames and values are file sizes
                If an archive was not read, then this dict is empty
                and the files should be examined on the OS filesystem

        returns a list (strings) of transformation names using the format: A~B
        for the transformation A->B
        """

        import glob
        import os.path
        import re

        edges=[]

        if len([fname for fname in arcfiles]) > 0:
            #
            # Look at files in archive
            #
            p = re.compile(r"dats/1/free_energy/(.*)_ambest/.*")
            for fname in arcfiles:
                m = p.match(fname)
                if m:
                    if m.group(1) not in edges:
                        edges.append(m.group(1))
        else:
            #
            # Look at files on OS filesystem
            #
            fs = glob.glob( "data/*~*" )
            for f in fs:
                cs = f.split("/")
                c = cs[-1]
                c = c.replace("-","~")
                edges.append(c)
        edges.sort()
        return edges



    def getstagedlambdas(self):
         """
         returns a dict. The dict value is a list. The list elements are
         strings. The strings are values/labels for the lambda values used
         for the transformation. The keys you choose here are the "stage"
         values in getdir(), above
         """
         lams = {}
EOF


#lams=($(gen_lambdas ${nlambda}))

printf '         lams["unified"] = [ ' >> gmbar.py

c=0

for lam in ${lams[@]};do
        if [ "${lam}" == "${lams[-1]}" ]; then
                printf "\"$lam\" ]\n" >> gmbar.py
        else
                printf "\"$lam\"," >> gmbar.py
        fi
        if [ "$c" -eq 2 ]; then
                printf "\n" >> gmbar.py
                printf '                        ' >> gmbar.py
                c=-1
        fi
        c=$(($c+1))
done
printf "\n" >> gmbar.py

cat << EOF >> gmbar.py
         return lams



# ######################################################################
#
# The rest of the file should work for any directory structure
#
# ######################################################################



if __name__ == "__main__":

    import os
    import sys
    import argparse
    import pickle
    from collections import OrderedDict as odict
    from collections import defaultdict as ddict

    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""
graphmbar input file generator and output file analyzer
""" )

    parser.add_argument \
        ("-f","--expt",
         help="Filename containing experimental dG values",
         type=str,
         required=True )

    parser.add_argument \
        ("-x","--exptcon",
         help="If present, then apply experimental ddG constraint between the two specified molecule names. Each instance of --exptcon takes two arguments.",
         nargs=2,
         type=str,
         action='append' )

    parser.add_argument \
        ("-l","--lead",
         help="string, name of lead ligand. If empty "", then the molecule with the most negative experimental free energy is chosen as the lead ligand.",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("-a","--analyze",
         help="If present, then analyze the output file, otherwise, write an input file",
         type=str,
         default="",
         required=False)

    parser.add_argument \
        ("-t","--trials",
         help="a list of trials (strings) to analyze. Default: [\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\"]",
         type=str,
         nargs='+',
         default=["1","2","3","4","5","6","7","8","9","10"],
         required=False)

    parser.add_argument \
        ("-b","--bar",
         help="If present, then only check for files necessary for BAR calculations, otherwise make sure all files necessary for MBAR are available",
         default=False,
         action='store_true',
         required=False)

    parser.add_argument \
        ("-z","--archive",
         help="If present, then read files from the compressed archive",
         type=str,
         default="",
         required=False)

    #################################

    parser.add_argument \
        ("-r","--rescyc",
         help="If present, then explicitly restraint minimum-length cycle constraints",
         action='store_true',
         required=False )

    parser.add_argument \
        ("-c","--concyc",
         help="If present, then implicitly constrain all cycles using slack variables",
         action='store_true',
         required=False )


    # parser.add_argument \
    #     ("-p","--pickle",
    #      help="If present, then read (if file exists) or write (if file does not exist) the information necessary to avoid examining the OS filesystem nor archive. Only useful if you plan on running the script on the same data many times",
    #      type=str,
    #      default="",
    #      required=False)


    # parser.add_argument \
    #     ("-W","--wang",
    #      help="If present, then output analysis should also include constraint enforcement based on the maximum likelihood optimization described in Wang, L.; Deng, Y.; Knight, J. L.; Wu, Y.; Kim, B.; Sherman, W.; Shelly, J. C.; Lin, T.; Abel, R.; J. Chem. Theory Comput. 2013, 9, 1282-1293.",
    #      default=False,
    #      action='store_true',
    #      required=False)


    args = parser.parse_args()

    ana = MySetup(args)

    if args.analyze:
        ana.AnalyzeOutput( args.analyze )
    else:
        ana.WriteInput(sys.stdout,args.rescyc,args.concyc)


#print(getstagedlambdas())
#print(getedges({}))
#print(getdir( "T", "e", "s", "t" ))

EOF
chmod a+x gmbar.py

}

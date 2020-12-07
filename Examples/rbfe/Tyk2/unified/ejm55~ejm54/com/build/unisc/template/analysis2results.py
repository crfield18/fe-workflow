#!/usr/bin/env python2.7
import os
from collections import defaultdict as ddict
prefix="unisc"
lams = [ "0.00000000","0.09090909","0.18181818","0.27272727","0.36363636","0.45454545","0.54545455","0.63636364","0.72727273","0.81818182","0.90909091","1.00000000" ]
merge_gaps=False

mdin = "TEMPLATE.mdin"
fh = file(mdin,"r")
numexchg=0
nstlim=None
ntwx=None
dt=None
for line in fh:
    cmdstr,sepstr,comstr = line.partition("!")
    if "ntpr" in cmdstr:
        cols = cmdstr.replace("="," ").replace(",","").strip().split()
        for icol in range(len(cols)-1):
            if cols[icol] == "ntpr":
                ntwx = int( cols[icol+1] )
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
        cols = cmdstr.replace("=","").replace(",","").strip().split()
        for icol in range(len(cols)-1):
            if cols[icol] == "nstlim":
                nstlim = int( cols[icol+1] )
                break

if ntwx is None:
    raise Exception("Could not determine ntwx from %s"%(mdin))

if dt is None:
    raise Exception("Could not determine dt from %s"%(mdin))

if nstlim is None:
    raise Exception("Could not determine nstlim from %s"%(mdin))

if numexchg < 1:
    numexchg = 1

dt = dt
nstep_per_sim = nstlim * numexchg
nframe_per_sim = nstep_per_sim / ntwx

if nstep_per_sim % ntwx != 0:
    print "num md steps per simulation is not a multiple of ntwx. Unclear how the simulation time works"

t_per_frame = dt * ntwx
t_per_sim = t_per_frame * nframe_per_sim

dvdl_data = ddict( lambda: ddict( float ) )
efep_data = ddict( lambda: ddict( lambda: ddict( float ) ) )

last_read_sim=0
missing_dirs=[]
for isim in range(1,100001):
    
    dirstr = "production/%06i"%(isim)
    if not os.path.isdir(dirstr):
        missing_dirs.append(dirstr)
        continue
    
    if last_read_sim != (isim-1):
        for d in missing_dirs:
            print "TIME GAP! Missing directory: %s"%(d)
        missing_dirs=[]

    t0 = (isim-1) * t_per_sim + t_per_frame

    missing_files=False
    error_msgs=[]
    error=False
    data = ddict(list)
    for lam in lams:
        nframe = 0
        dat = "%s/dvdl_%s.dat"%(dirstr,lam)
        if not os.path.isfile(dat):
            error=True
            missing_files=True
        else:
            fh = file(dat,"r")
            for line in fh:
                cols = line.strip().split()
                if len(cols) == 2:
                    nframe += 1
                    data[lam].append( cols[-1] )
        if nframe != nframe_per_sim and nframe-1 != nframe_per_sim:
            msg="%s expected %i frames, but found %i"%(dat,nframe_per_sim,nframe)
            error_msgs.append(msg)
            error=True
    if not error:
        for iframe in range(nframe_per_sim):
            t = t0 + iframe * t_per_frame
            for lam in lams:
                dvdl_data[t][lam] = data[lam][iframe]

    if missing_files and len(error_msgs) == len(lams):
        print "%s doesn't appear to have been analyzed yet"%(dirstr)
    else:
        for msg in error_msgs:
            print msg
    if len(error_msgs) > 0:
        missing_dirs.append(dirstr)
        continue
            

    
    missing_files=False
    error_msgs=[]
    error=False
    data = ddict(lambda: ddict(list))
    for tlam in lams:
        for plam in lams:
            nframe = 0
            dat = "%s/efep_%s_%s.dat"%(dirstr,tlam,plam)
            if not os.path.isfile(dat):
                error=True
                missing_files=True
            else:
                fh = file(dat,"r")
                for line in fh:
                    cols = line.strip().split()
                    if len(cols) == 2:
                        nframe += 1
                        data[tlam][plam].append( cols[-1] )
            if nframe != nframe_per_sim and nframe-1 != nframe_per_sim:
                msg="%s expected %i frames, but found %i"%(dat,nframe_per_sim,nframe)
                error_msgs.append(msg)
                error=True
    if not error:
        for iframe in range(nframe_per_sim):
            t = t0 + iframe * t_per_frame
            for tlam in lams:
                for plam in lams:
                    efep_data[t][tlam][plam] = data[tlam][plam][iframe]
        

    for msg in error_msgs:
        print msg
    if len(error_msgs) > 0:
        missing_dirs.append(dirstr)
        continue

    
if not os.path.exists("results/data"):
    os.makedirs("results/data")
    
ts=[ t for t in dvdl_data ]
if len( ts ) > 0:
    for lam in lams:
        dat = "results/data/dvdl_%s.dat"%(lam)
        fh = file(dat,"w")
        for i,t in enumerate(sorted(dvdl_data)):
            time=t
            if merge_gaps:
                time = (i+1)*t_per_frame
            fh.write("%12.1f %s\n"%(t,dvdl_data[t][lam]))
        fh.close()
            
ts=[ t for t in dvdl_data ]
if len( ts ) > 0: 
    for tlam in lams:
        for plam in lams:
            dat = "results/data/efep_%s_%s.dat"%(tlam,plam)
            fh = file(dat,"w")
            for i,t in enumerate(sorted(efep_data)):
                time=t
                if merge_gaps:
                    time = (i+1)*t_per_frame
                fh.write("%12.1f %s\n"%(t,efep_data[t][tlam][plam]))
            fh.close()


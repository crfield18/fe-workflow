#!/usr/bin/env python2.7

import parmed
#import parmutils

##############################################################
##############################################################
##############################################################
##############################################################
##############################################################



def OpenParm( fname, xyz=None ):
    import parmed
    from parmed.constants import IFBOX
    param = parmed.load_file( fname,xyz=xyz )
    if param.box is not None:
        if abs(param.box[3]-109.471219)<1.e-4 and \
           abs(param.box[4]-109.471219)<1.e-4 and \
           abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param


def SaveParm( param, fname ):
    from parmed.constants import IFBOX
    for a in param.atoms:
        param.parm_data["CHARGE"][ a.idx ] = a.charge
    if param.box is not None:
       if abs(param.box[3]-109.471219)<1.e-4 and \
          abs(param.box[4]-109.471219)<1.e-4 and \
          abs(param.box[5]-109.471219)<1.e-4:
           param.parm_data["POINTERS"][IFBOX]=2
           param.pointers["IFBOX"]=2
    try:
        param.save( fname, overwrite=True )
    except:
        param.save( fname )


def CopyParm( parm ):
    import copy
    parm.remake_parm()
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    p.hasbox = copy.copy( parm.hasbox )
    return p


def Strip( parm, mask ):
    p = CopyParm( parm )
    p.strip( "%s"%(mask) )
    return p


def Extract( parm, mask ):
    return Strip( parm, "!(%s)"%(mask) )


def Join( p, q ):
    import numpy
    import copy
    if q.coordinates.shape[0] == 0:
        pq=p
    elif p.coordinates.shape[0] == 0:
        pq=q
    else:
        pq = p + q    
        pq.coordinates = numpy.concatenate( (p.coordinates, q.coordinates) )
        pq.box = copy.copy( p.box )
        pq.hasbox = copy.copy( p.hasbox )
    return pq


def ReindexFromDeletions( idxs, deletions ):
    newidxs = []
    for idx in idxs:
        if idx not in deletions:
            newidx = idx
            for d in deletions:
                if d < idx:
                    newidx -= 1
            newidxs.append(newidx)
    return newidxs

def mark_path(atom,dist,path,scsel):
    #print("marking path ",atom.idx,path)
    mark_partners=True
    if atom.idx in path:
        if path[atom.idx] <= dist:
            mark_partners=False
        else:
            path[atom.idx] = dist
    else:
        path[atom.idx] = dist
        
    if mark_partners:
        for b in atom.bond_partners:
            if b.idx in scsel:
                if b.atomic_number != 1:
                    mark_path(b,dist+1,path,scsel)
                
def find_atoms_to_delete( p, scsel ):
    from collections import defaultdict as ddict
    import copy
    ncpath = ddict( int )
    hatoms = ddict( list )
    for idx in scsel:
        a = p.atoms[idx]
        if a.atomic_number != 1:
            hatoms[idx] = [idx]
            for b in a.bond_partners:
                if b.idx in scsel:
                    if b.atomic_number == 1:
                        hatoms[idx].append( b.idx )
                else:
                    mark_path(a,0,ncpath,scsel)

    maxval=0
    for i in ncpath:
        maxval = max(ncpath[i],maxval)
    for i in hatoms:
        if i not in ncpath:
            mark_path(p.atoms[i],maxval,ncpath,scsel)

    maxval=0
    for i in ncpath:
        maxval = max(ncpath[i],maxval)
    servals = []
    seenvals = []
    for target in range(maxval,-1,-1):
        highvals = []
        for i in hatoms:
            if ncpath[i] == target:
                #print("hatoms",hatoms[i])
                highvals.extend( hatoms[i] )
        if len(highvals) > 0:
            servals.append(highvals)
            seenvals.extend(highvals)
    highvals=[]
    for i in scsel:
        if i not in seenvals:
            highvals.append(i)
    if len(highvals) > 0:
        servals.append(highvals)
    return servals

                    
def GetSelectedAtomIndices(param,maskstr):
    #param = parmed.load_file(parmfile)
    #mask = parmed.amber.mask.AmberMask( param, maskstr )
    #aidxs = mask.Selected()
    #for aidx in aidxs:
    #    atom = param.atoms[aidx]
    #    res  = atom.residue
    sele = []
    if len(maskstr) > 0:
        newmaskstr = maskstr.replace("@0","!@*")
        # newmaskstr = maskstr.replace("@0","").replace("()","")
        # for c in ["|","&"]:
        #     if len(newmaskstr) > 0:
        #         if newmaskstr[0] == c:
        #             newmaskstr = newmaskstr[1:]
                    
        # for c in ["|","&"]:
        #     if len(newmaskstr) > 0:
        #         if newmaskstr[-1] == c:
        #             newmaskstr = newmaskstr[:-1]
                
        #print "newmaskstr=",newmaskstr
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


def AddNewBondType(parm,sel1,sel2):
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from parmed.topologyobjects import Bond,Angle,Dihedral
    import copy
    
#    sel1 = parmed.amber.mask.AmberMask(parm,mask1).Selection()
#    sel2 = parmed.amber.mask.AmberMask(parm,mask2).Selection()
#    foo = GetSelectedAtomIndices(parm,mask1)
    if len(sel1) != len(sel2):
        raise Exception('parmutils.AddNewBondType: Each mask must select the same '
                        'number of atoms!')

    # If no atoms, nothing to do
    if len(sel1) == 0: return

    new_bnd_typ = None
    #atnum1, atnum2 = -1, -1
    # Loop through all of the selected atoms
    for atnum1,atnum2 in zip(sel1,sel2):
        # Collect the atoms involved
        #atnum1 = sel1.index(1, atnum1+1)
        #atnum2 = sel2.index(1, atnum2+1)
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        # See if the bond exists in the first place, and if so, replace its
        # bond type with our new bond type (new_bnd)
        if atm2 in atm1.bond_partners and atm1 in atm2.bond_partners:
            for bond in atm1.bonds:
                if atm2 in bond:
                    # print "found %4s@%-4s - %4s@%-4s : %7.4f %7.4f"%\
                    #     ( atm1.residue.name,atm1.name,
                    #       atm2.residue.name,atm2.name,
                    #       bond.type.req,bond.type.k)
                    new_bnd_typ = copy.copy(bond.type)
                    break

    
    exists = False
    # If the bond is new, add it to the type list
    if not exists:
        parm.bond_types.append(new_bnd_typ)
        new_bnd_typ.list = parm.bond_types

    #atnum1, atnum2 = -1, -1
    # Loop through all of the selected atoms
    for atnum1,atnum2 in zip(sel1,sel2):
        # Collect the atoms involved
        #atnum1 = sel1.index(1, atnum1+1)
        #atnum2 = sel2.index(1, atnum2+1)
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]

        # See if the bond exists in the first place, and if so, replace its
        # bond type with our new bond type (new_bnd)
        if atm2 in atm1.bond_partners and atm1 in atm2.bond_partners:
            for bond in atm1.bonds:
                if atm2 in bond:
                    bond.type = new_bnd_typ
                    parm.bonds.changed = True
                    break
            # Otherwise, it doesn't exist, so we just create a new one
        else:
            parm.bonds.append(Bond(atm1, atm2, new_bnd_typ))
    # Make sure we update 1-4 exception handling if we created any rings
    parm.update_dihedral_exclusions()

    


def AddNewAngleType(parm,sel1,sel2,sel3):
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from parmed.topologyobjects import Bond,Angle,Dihedral
    import copy
    
    #sel1 = parmed.amber.mask.AmberMask(parm,mask1).Selection()
    #sel2 = parmed.amber.mask.AmberMask(parm,mask2).Selection()
    #sel3 = parmed.amber.mask.AmberMask(parm,mask3).Selection()

    if len(sel1) != len(sel2):
        raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
                        'number of atoms!')

    if len(sel1) != len(sel3):
        raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
                        'number of atoms!')

    # If no atoms, nothing to do
    if len(sel1) == 0: return


    new_ang_typ = None
    #atnum1, atnum2, atnum3 = -1, -1, -1
    # Loop through all of the selections
    for ii in range(len(sel1)):
        # Collect the atoms involved
        atnum1 = sel1[ii]
        atnum2 = sel2[ii]
        atnum3 = sel3[ii]
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        atm3 = parm.atoms[atnum3]

        # See if the angle exists in the first place, and if so, replace its
        # angle type with our new angle type (new_ang)
        found = False
        if atm1 in atm3.angle_partners:
            for ang in atm1.angles:
                if atm2 in ang and atm3 in ang:
                    new_ang_typ = copy.copy(ang.type)
                    break
    
    exists = False
    # If the angle is new, add it to the type list
    if not exists:
        parm.angle_types.append(new_ang_typ)
        new_ang_typ.list = parm.angle_types

    # Loop through all of the selections
    for ii in range(len(sel1)):
        # Collect the atoms involved
        atnum1 = sel1[ii]
        atnum2 = sel2[ii]
        atnum3 = sel3[ii]
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        atm3 = parm.atoms[atnum3]

        # See if the angle exists in the first place, and if so, replace its
        # angle type with our new angle type (new_ang)
        found = False
        if atm1 in atm3.angle_partners:
            for ang in atm1.angles:
                if atm2 in ang and atm3 in ang:
                    ang.type = new_ang_typ
                    parm.angles.changed = True
                    found = True
                    break
        # If not found, create a new angle
        if not found:
            parm.angles.append(Angle(atm1, atm2, atm3, new_ang_typ))
    # Make sure we update 1-4 exception handling if we created any rings
    parm.update_dihedral_exclusions()
    



##############################################################
##############################################################
##############################################################
##############################################################
##############################################################



class AmberCmd(object):
    #
    # these are default values
    #
    exe    = "${AMBERHOME}/bin/sander.MPI"
    parm7  = "prmtop"
    mdin   = "mdin"
    mdout  = "mdout"
    crd7   = "inpcrd"
    rst7   = "restrt"
    refc   = "refc"
    nc     = "mdcrd"
    mdinfo = "mdinfo"

    def __init__(self,base = None, inpcrd = None):
        self.SetDefaults()
        if base != None: self.SetBaseName( base )
        if inpcrd != None: self.CRD7 = inpcrd

    def SetDefaults(self):
        self.reference = False
        self.EXE    = AmberCmd.exe
        self.PARM7  = AmberCmd.parm7
        self.MDIN   = AmberCmd.mdin
        self.MDOUT  = AmberCmd.mdout
        self.CRD7   = AmberCmd.crd7
        self.RST7   = AmberCmd.rst7
        self.REFC   = AmberCmd.refc
        self.NC     = AmberCmd.nc
        self.MDINFO = AmberCmd.mdinfo
        self.INPNC  = None
        self.gsize  = None
        self.gmmsize = None

    def SetBaseName(self,base):
        self.BASE   = base
        self.MDIN   = base + ".mdin"
        self.MDOUT  = base + ".mdout"
        self.RST7   = base + ".rst7"
        self.NC     = base + ".nc"
        self.MDINFO = base + ".mdinfo"
        self.DISANG = base + ".disang"
        self.DUMPAVE= base + ".dumpave"

    def SetGroupSize(self,gsize,gmmsize=None):
        self.gsize = gsize
        if gmmsize is not None:
            self.gmmsize = gmmsize
            
    def CmdString(self):
        crd = self.CRD7
        if crd == self.RST7:
            crd = self.RST7.replace(".rst7",".crd7")
        if self.gsize is not None:
            cmd = "-gsize %4i"%(self.gsize)
            if self.gmmsize is not None:
                cmd += " -gmmsize %4i"%(self.gmmsize)
            else:
                cmd += " -gmmsize %4i"%(self.gsize)
        else:
            cmd = ""
        cmd += " -O -p " + self.PARM7 \
               + " -i " + self.MDIN + " -o " + self.MDOUT \
               + " -c " + crd + " -r " + self.RST7 \
               + " -x " + self.NC \
               + " -inf " + self.MDINFO
        if self.REFC != "refc":
            cmd += " -ref " + self.REFC
        if self.INPNC is not None:
            cmd += " -y " + self.INPNC
            
        return cmd

class Disang1D(object):
    def __init__(self,template,r1,k1):
        self.TEMPLATE = template
        self.R1 = r1
        self.K1 = k1

    def WriteDisang(self,fname):
        cout = open(fname,"w")
        cin = file(self.TEMPLATE,"r")
        for line in cin:
            line = re.sub(r'R1',"%.2f"%(self.R1),line)
            line = re.sub(r'K1',"%.2f"%(self.K1),line)
            cout.write(line)
        cin.close()
        cout.close()
    
class Disang2D(object):
    def __init__(self,template,r1,k1,r2,k2):
        self.TEMPLATE = template
        self.R1 = r1
        self.K1 = k1
        self.R2 = r2
        self.K2 = k2

    def WriteDisang(self,fname):
        cout = open(fname,"w")
        cin = open(self.TEMPLATE,"r")
        for line in cin:
            line = re.sub(r'R1',"%.2f"%(self.R1),line)
            line = re.sub(r'K1',"%.2f"%(self.K1),line)
            line = re.sub(r'R2',"%.2f"%(self.R2),line)
            line = re.sub(r'K2',"%.2f"%(self.K2),line)
            cout.write(line)
        cin.close()
        cout.close()
    

class Mdin(AmberCmd):
    def __init__(self,base = None, inpcrd = None):
        from collections import defaultdict as ddict
        AmberCmd.__init__(self,base,inpcrd)
        self.title="title"
        self.cntrl = ddict( str )
        self.qmmm = ddict( str )
        self.ewald = ddict( str )
        self.dlfind = ddict( str )
        self.shake = False # used for dlfind
        self.restraints = None
        self.DUMPFREQ = 25
        self.cntrl["ntxo"]=1
        self.cntrl["cut"]=12.0
        self.cntrl["nstlim"]=1000000
        self.cntrl["dt"]=0.001
        self.cntrl["iwrap"]=1
        self.cntrl["ntf"]=2
        self.cntrl["ntc"]=2
        self.cntrl["ioutfm"]=1
        self.cntrl["ig"]=-1
        self.Set_Restart(True)
        self.Set_PrintFreq(1000)


    def Set_DumpFreq(self,dumpfreq):
        self.DUMPFREQ = dumpfreq
        
    def Set_PrintFreq(self,freq):
        self.cntrl["ntpr"]=freq
        self.cntrl["ntwr"]=freq
        self.cntrl["ntwx"]=freq

    def Set_Restart(self,restart=True):
        if not restart:
            self.cntrl["imin"]=0
            self.cntrl["ntx"]=1
            self.cntrl["irest"]=0
        else:
            self.cntrl["imin"]=0
            self.cntrl["ntx"]=5
            self.cntrl["irest"]=1
            
    def Set_NVT(self,temp0=298.0):
        self.cntrl["ntt"]=3
        self.cntrl["ntb"]=1
        self.cntrl["ntp"]=0
        self.cntrl["temp0"]=temp0
        self.cntrl["gamma_ln"]=5.0
        self.cntrl.pop("pres0",None)
        self.cntrl.pop("taup",None)
        self.cntrl.pop("nscm",None)

    def Set_NPT(self,temp0=298.0):
        self.cntrl["ntt"] = 3
        self.cntrl["ntb"] = 2
        self.cntrl["ntp"] = 1
        self.cntrl["temp0"] = temp0
        self.cntrl["gamma_ln"] = 5.0
        self.cntrl["barostat"] = 2
        self.cntrl["pres0"] = 1.013
        self.cntrl["taup"] = 2.0
        self.cntrl.pop("nscm",None)

    def Set_NVE(self):
        self.cntrl["ntt"]=0
        self.cntrl["ntb"]=1
        self.cntrl["ntp"]=0
        self.cntrl.pop("temp0",None)
        self.cntrl.pop("gamma_ln",None)
        self.cntrl.pop("pres0",None)
        self.cntrl.pop("taup",None)
        self.cntrl["nscm"]=0

        
    def Set_QMMM_AM1(self,qmmask=None,qmcharge=None,tight=False):
        self.qmmm["qm_theory"] = "'AM1D'"
        self.qmmm["qmmm_switch"]=1
        self.qmmm["qm_ewald"]=1
        self.qmmm["qmshake"]=0
        if tight:
            self.qmmm["diag_routine"] = 6
            self.qmmm["tight_p_conv"] = 1
            self.qmmm["scfconv"]=1.e-11

        if qmmask is not None:
            self.qmmm["qmmask"] = qmmask
        if qmcharge is not None:
            self.qmmm["qmcharge"] = qmcharge
            
        self.qmmm.pop("hfdf_theory",None)
        self.qmmm.pop("hfdf_basis",None)
        self.qmmm.pop("hfdf_mempercore",None)
        self.qmmm.pop("hfdf_ewald",None)
        self.qmmm.pop("hfdf_mm_percent",None)
        self.qmmm.pop("hfdf_qmmm_wswitch",None)


    def Set_QMMM_PBE0(self,qmmask=None,qmcharge=None,basis=None,tight=False):
        if qmmask is not None:
            self.qmmm["qmmask"] = qmmask
        if qmcharge is not None:
            self.qmmm["qmcharge"] = qmcharge
            
        self.qmmm["qmshake"]=0
        self.qmmm["qm_theory"] = "'HFDF'"
        self.qmmm["hfdf_theory"] = "'PBE0'"
        if basis is not None:
            self.qmmm["hfdf_basis"] = basis
        else:
            self.qmmm["hfdf_basis"] = "'6-31G*'"
        self.qmmm["scfconv"] = 1.e-7
        if tight:
            self.qmmm["scfconv"] = 1.e-8
            
        self.qmmm["hfdf_mempercore"] = 2000
        self.qmmm["hfdf_ewald"] = "T"
        self.qmmm["qm_ewald"] = 1
        self.qmmm["hfdf_mm_percent"] = 0.0
        self.qmmm["hfdf_qmmm_wswitch"] = 0.0
        self.qmmm["qmmm_switch"] = 0

        self.qmmm.pop("tight_p_conv",None)
        self.qmmm.pop("diag_routine",None)

    def Set_DLFIND_Minimize(self):
        self.dlfind["crdrep"] = "'HDLC'"
        self.dlfind["optalg"] = "'LBFGS'"
        self.dlfind["trustrad"] = "'GRAD'"
        self.dlfind["hessini"] = "'IDENTITY'"
        self.dlfind["hessupd"] = "'BFGS'"
        self.dlfind["neb"] = "''"
        self.dlfind["dimer"] = "''"
        self.dlfind["crdfile"] = "''"
        if "wtmask" not in self.dlfind:
            self.dlfind["wtmask"] = "''"
        self.dlfind["hessout"] = "''"
        self.dlfind["ihess"] = -1
        self.dlfind["nhess"] = -1
        self.dlfind["tol"] = 5.e-6
        self.dlfind["tole"] = 1.e-4
        self.dlfind["maxcycle"] = 1200
        self.dlfind["maxene"] = 2400
        self.dlfind["maxstep"] = -1.
        self.dlfind["minstep"] = 1.e-8
        self.dlfind["scalestep"] = 0.5
        self.dlfind["lbfgsmem"] = -1
        self.dlfind["maxupdate"] = -1
        self.dlfind["hessdel"] = 0.005
        self.dlfind["dimer_delta"] = 0.01
        self.dlfind["dimer_maxrot"] = 25
        self.dlfind["dimer_tolrot"] = 2.
        self.dlfind["dimer_mmpercent"] = 0.0

        self.dlfind.pop("crdfile",None)

        
    def Set_DLFIND_TransitionState(self,displaced_rst7):
        self.Set_DLFIND_Minimize()
        self.dlfind["trustrad"] = "'CONST'"
        self.dlfind["hessupd"] = "'BOFILL'"
        self.dlfind["dimer"] = "'LN_NO_EXTRA'"
        self.dlfind["crdfile"] = "'%s'"%(displaced_rst7) 
        
            
    def Set_DLFIND(self,active = None):
        self.cntrl = ddict( str )
        self.qmmm = ddict( str )
        self.ewald = ddict( str )
        self.dlfind = ddict( str )
        self.cntrl["imin"]   = 1  # minimize
        self.cntrl["ntmin"]  = 5  # read &dlfind
        self.cntrl["irest"]  = 0  # no need to read forces from rst
        self.cntrl["ntx"]    = 1
        self.cntrl["ntwx"]   = 50 # traj freq
        self.cntrl["ioutfm"] = 1
        self.cntrl["ntb"]    = 1  # periodic
        self.cntrl["iwrap"]  = 1  # wrap
        self.cntrl["cut"]    = 12.0 # nonbond
        self.cntrl["ig"]     = -1   # random number seed
        self.cntrl["ifqnt"]  = 1    # read &qmmm
        self.cntrl["nmropt"] = 1    # read DISANG,DUMPAVE
        self.cntrl["ntf"]    = 2    # shake
        self.cntrl["ntc"]    = 2    # shake
        if "qmmask" not in self.qmmm:
            self.qmmm["qmmask"] = "''"
        if "qmcharge" not in self.qmmm:
            self.qmmm["qmcharge"]  = 0
        self.qmmm["writepdb"]  = 0
        self.qmmm["qmshake"]   = 0
        self.qmmm["verbosity"] = 0
        self.Set_QMMM_AM1()
        self.dlfind["active"] = "''"
        self.dlfind["prefix"] = ""
        self.Set_DLFIND_Minimize()
        if active is not None:
            self.dlfind["active"] = active

    def WriteMdin(self):
        from collections import defaultdict as ddict
        fh = file(self.MDIN,"w")
        fh.write("%s\n"%(self.title))
        fh.write("&cntrl\n")
        iokeys = ddict(str)
        iokeys["irest"]  = "! 0 = start, 1 = restart"
        iokeys["ntx"]    = "! 1 = start, 5 = restart"
        iokeys["ntxo"]   = "! read/write rst as formatted file"
        iokeys["iwrap"]  = "! wrap crds to unit cell"
        iokeys["ioutfm"] = "! write mdcrd as netcdf"
        iokeys["ntpr"]   = "! mdout print freq"
        iokeys["ntwx"]   = "! mdcrd print freq"
        iokeys["ntwr"]   = "! rst print freq"

        runkeys = ddict(str)
        runkeys["imin"]   = "! 0=dynamics, 1=internal minimizer"
        runkeys["ntmin"]  = "! minimization algo (5=dlfind)" 
        runkeys["nstlim"] = "! number of time steps"
        runkeys["dt"]     = "! ps/step"
        runkeys["numexch"]= "! number of replica exchange attempts"
        runkeys["ntb"]    = "! 1=NVT, 2=NPT, 0=no box"

        tempkeys = ddict(str)
        tempkeys["ntt"]      = "! thermostat (3=Langevin)"
        tempkeys["tempi"]    = "! initial temp for generating velocities"
        tempkeys["temp0"]    = "! target temp"
        tempkeys["gamma_ln"] = "! Langevin collision freq"

        preskeys = ddict(str)
        preskeys["ntp"]      = "! 0=no scaling, 1=isotropic, 2=anisotropic"
        preskeys["barostat"] = "! barostat (1=Berendsen, 2=MC)"
        preskeys["pres0"]    = "! pressure (bar), 1.013 bar/atm"
        preskeys["taup"]     = "! pressure relaxation time"

        shakekeys = ddict(str)
        shakekeys["ntf"] = "! 1=cpt all bond E, 2=ignore HX bond E, 3=ignore all bond E"
        shakekeys["ntc"] = "! 1=no shake, 2=HX constrained, 3=all constrained"
        shakekeys["noshakemask"] = "! do not shake these"

        namedkeys = []
        for ktypes in [ iokeys, runkeys, tempkeys, preskeys, shakekeys ]:
            for key in ktypes:
                namedkeys.append( key )
        
        fh.write("! IO =======================================\n")
        for key in ["irest","ntx","ntxo","iwrap","ioutfm","ntpr","ntwx","ntwr"]:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),iokeys[key]))
        fh.write("! DYNAMICS =================================\n")
        for key in runkeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),runkeys[key]))
        fh.write("! TEMPERATURE ==============================\n")
        for key in tempkeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),tempkeys[key]))
        fh.write("! PRESSURE  ================================\n")
        for key in preskeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),preskeys[key]))
        fh.write("! SHAKE ====================================\n")
        for key in shakekeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),shakekeys[key]))
        is_default_not_shaked = False
        if "ntf" in self.cntrl:
            if self.cntrl["ntf"] == 1:
                is_default_not_shaked = True
        if not self.shake:
            if ("noshakemask" not in self.cntrl) and (not is_default_not_shaked):
                if "active" in self.dlfind:
                    fh.write("%11s = %s\n"%("noshakemask",self.dlfind["active"]))
        fh.write("! MISC =====================================\n")
        for key in sorted(self.cntrl):
            if key not in namedkeys:
                fh.write("%11s = %s\n"%(key,str(self.cntrl[key])))

        fh.write("/\n\n")

        if "nmropt" in sorted(self.cntrl):
            if self.cntrl["nmropt"] == 1:
                fh.write("&wt type='DUMPFREQ', istep1 = %i, /\n"%(self.DUMPFREQ))
                fh.write("&wt type='END', /\n")
                fh.write("LISTOUT=POUT\n")
                fh.write("DISANG=%s\n"%(self.DISANG))
                fh.write("DUMPAVE=%s\n\n"%(self.DUMPAVE))
                if self.restraints is not None:
                    self.restraints.WriteDisang( self.DISANG )
                elif not os.path.isfile( self.DISANG ):
                    raise Exception("%s not found!\n"%(self.DISANG))

        ewald_keys = [ key for key in self.ewald ]
        if len(ewald_keys) > 0:
            fh.write("&ewald\n")
            for key in sorted(self.ewald):
                fh.write("%14s = %s\n"%(key,str(self.ewald[key])))
            fh.write("/\n\n")
                
        if "ifqnt" in self.cntrl:
            if self.cntrl["ifqnt"] > 0:
                fh.write("&qmmm\n")
                for key in sorted(self.qmmm):
                    fh.write("%14s = %s\n"%(key,str(self.qmmm[key])))
                fh.write("/\n\n")
        if "ntmin" in self.cntrl:
            if self.cntrl["ntmin"] == 5:
                self.dlfind["prefix"] = "'%s'"%(self.BASE)
                fh.write("&dlfind\n")
                for key in sorted(self.dlfind):
                    fh.write("%14s = %s\n"%(key,str(self.dlfind[key])))
                fh.write("/\n\n")
                
        fh.close()




##############################################################
##############################################################
##############################################################
##############################################################
##############################################################

def WriteFrcmod(param,native_frcmod,uniqueparams=False):
    self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
    WriteFrcmodObj(self,native_frcmod,angfact=0.9999995714245039,uniqueparams=uniqueparams)

def WriteFrcmodObj(self,native_frcmod,angfact=1.0,uniqueparams=False,selected_names=None,changed_frcmod=None):
    
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from collections import defaultdict as ddict
    import re

    if uniqueparams and selected_names is None:
        selected_names = {}
        for atom, typ in iteritems(self.atom_types):
            selected_names[atom]=1
    elif selected_names is None:
        selected_names = {}
        

    if True:

#        angfact = 0.9999995714245039

        class combofile(object):
            def __init__(self,fh1,fh2):
                self.fh1 = fh1
                self.fh2 = fh2
            def write(self,s):
                self.fh1.write(s)
                self.fh2.write(s)
        
        nfile = file(native_frcmod,"w")
        if changed_frcmod is None:
            cfile = nfile
            outfile = nfile
        else:
            cfile = file(changed_frcmod,"w")
            outfile = combofile( nfile, cfile )

#        self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
        outfile.write("modified parameters")
        outfile.write('\n')
        # Write the atom mass
        outfile.write('MASS\n')
        for atom, typ in iteritems(self.atom_types):
            fh=nfile
            if atom in selected_names:
                fh = cfile
            fh.write('%s%11.8f\n' % (atom.ljust(6), typ.mass))
                
        outfile.write('\n')
        # Write the bonds
        outfile.write('BOND\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict(float) )
        for (a1, a2), typ in iteritems(self.bond_types):
            typ.k = float("%.8f"%(typ.k))
            fh=nfile
            delta = 0
            if a1 in selected_names or a2 in selected_names:
                fh=cfile
                qq = (a1,a2)
                if qq in cdone: continue
                qq = (a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                deltas[typ.k][typ.req] += 1.e-13
                delta = deltas[typ.k][typ.req]
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            fh.write('%s-%s   %19.14f  %11.8f\n' %
                     (a1.ljust(2), a2.ljust(2), typ.k+delta, typ.req))
        outfile.write('\n')
        # Write the angles
        outfile.write('ANGLE\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict(float) )
        for (a1, a2, a3), typ in iteritems(self.angle_types):
            typ.k = float("%.8f"%(typ.k))
            delta = 0.
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names:
                fh=cfile
                qq = (a1,a2,a3)
                if qq in cdone: continue
                qq = (a3,a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                deltas[typ.k][typ.theteq] += 1.e-13
                delta = deltas[typ.k][typ.theteq]
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            fh.write('%s-%s-%s   %19.14f  %17.3f\n' %
                     (a1.ljust(2), a2.ljust(2), a3.ljust(2), typ.k+delta,
                      typ.theteq * angfact))
        outfile.write('\n')
        # Write the dihedrals
        outfile.write('DIHE\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
        for (a1, a2, a3, a4), typ in iteritems(self.dihedral_types):
            isnew = False
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names or a4 in selected_names:
                fh=cfile
                qq = (a1,a2,a3,a4)
                if qq in cdone: continue
                qq = (a4,a3,a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                isnew = True
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            if isinstance(typ, DihedralType) or len(typ) == 1:
                if not isinstance(typ, DihedralType):
                    typ = typ[0]
                    typ.phi_k = float("%.8f"%(typ.phi_k))
                    delta = 0
                    if isnew:
                        deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
                        delta = deltas[typ.phi_k][typ.phase][typ.per]
                if abs(typ.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
            else:
                typ = sorted( typ, key=lambda x: x.per, reverse=False )
                for dtyp in typ[:-1]:
                    dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                    delta = 0
                    if isnew:
                        deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
                        delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
                    if abs(dtyp.phase-180) < 0.0001:
                        #print "%20.16f"%(180.0/dtyp.phase)
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
                                                      a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                    else:
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
                                                      a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                dtyp = typ[-1]
                dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                delta = 0
                if isnew:
                    deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
                    delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
                if abs(dtyp.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                    
        outfile.write('\n')
        # Write the impropers
        deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
        outfile.write('IMPROPER\n')
        for (a1, a2, a3, a4), typ in iteritems(self.improper_periodic_types):
            # Make sure wild-cards come at the beginning
            if a2 == 'X':
                assert a4 == 'X', 'Malformed generic improper!'
                a1, a2, a3, a4 = a2, a4, a3, a1
            elif a4 == 'X':
                a1, a2, a3, a4 = a4, a1, a3, a2

            typ.phi_k = float("%.8f"%(typ.phi_k))
            delta = 0
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names or a4 in selected_names:
                fh=cfile
                deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
                delta = deltas[typ.phi_k][typ.phase][typ.per]
            else:
                fh=nfile
            if abs(typ.phase-180) < 0.0001:
                fh.write('%s-%s-%s-%s %20.14f %13.3f %5.1f\n' %
                         (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
                          typ.phi_k+delta, typ.phase * angfact, typ.per))
            else:
                fh.write('%s-%s-%s-%s %20.14f %13.8f %5.1f\n' %
                         (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
                          typ.phi_k+delta, typ.phase * angfact, typ.per))

                
        outfile.write('\n')
        # Write the LJ terms

        deltas = ddict( lambda: ddict( float ) )

        outfile.write('NONB\n')
        for atom, typ in iteritems(self.atom_types):
            #typ.rmin = float("%.8f"%(typ.rmin))
            typ.epsilon = float("%.9f"%(typ.epsilon))
            delta = 0.
            if atom in selected_names:
                fh=cfile
                deltas[typ.rmin][typ.epsilon] += 1.e-13
                delta = deltas[typ.rmin][typ.epsilon]
            else:
                fh=nfile
            if delta == 0.:
                fh.write('%-3s  %12.8f %18.9f\n' %
                         (atom.ljust(2), typ.rmin, typ.epsilon))
            else:
                fh.write('%-3s  %12.8f %18.14f\n' %
                         (atom.ljust(2), typ.rmin, typ.epsilon+delta))
        outfile.write('\n')
        # Write the NBFIX terms
        if self.nbfix_types:
            outfile.write('LJEDIT\n')
            for (a1, a2), (eps, rmin) in iteritems(self.nbfix_types):
                if a1 in selected_names or a2 in selected_names:
                    fh=cfile
                else:
                    fh=nfile
                fh.write('%s %s %13.8f %13.8f %13.8f %13.8f\n' %
                         (a1.ljust(2), a2.ljust(2), eps, rmin/2,
                          eps, rmin/2))
        cfile.close()
        nfile.close()
        


def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]



def ticopy( parm, molmask, base="ticopy" ):
    import copy
    import os

    fh = file("%s.sh"%(base),"w")
    fh.write("""#!/bin/bash

cat << 'EOF' > %s.sh.cmds

"""%(base))
    
    parmed.tools.writeOFF( parm, "%s.lib"%(base) ).execute()
    WriteFrcmod( parm, "%s.frcmod"%(base), uniqueparams=False )

    fh.write("""

# source leaprc.protein.ff14SB
# source leaprc.water.tip4pew
# source leaprc.constph
# set default PBradii mbondi3

""")
    fh.write("loadOff %s.lib\n"%(base))
    fh.write("loadAmberParams %s.frcmod\n\n\n"%(base))

    
    p = CopyParm( parm )
    p.strip( "!(%s)"%( molmask ) )
    parmed.tools.writeCoordinates(p, "%s.mol0.pdb"%(base)).execute()
    parmed.tools.writeCoordinates(p, "%s.mol1.pdb"%(base)).execute()

    seq = GetResSeq( p )
    fh.write("mol0 = loadPdbUsingSeq %s.mol0.pdb { %s }\n\n"%(base," ".join(seq)))
    fh.write("mol1 = loadPdbUsingSeq %s.mol1.pdb { %s }\n\n\n"%(base," ".join(seq)))

    p = CopyParm( parm )
    p.strip( "(%s)|:WAT"%( molmask ) )
    nonwat_len = len(p.atoms)
    if nonwat_len > 0:
        parmed.tools.writeCoordinates(p, "%s.nonwater.pdb"%(base)).execute()
        seq = GetResSeq( p )
        fh.write("nonwater = loadPdbUsingSeq %s.nonwater.pdb { %s }\n\n"%(base," ".join(seq)))
        
    p = CopyParm( parm )
    p.strip( "(%s)|(!:WAT)"%( molmask ) )
    wat_len = len(p.atoms)
    if wat_len > 0:
        parmed.tools.writeCoordinates(p, "%s.solvent.pdb"%(base)).execute()
        fh.write("solvent = loadPdb %s.solvent.pdb\n\n"%(base))
        
    fh.write("\n\nx = combine { mol0 mol1")
    if nonwat_len > 0:
        fh.write(" nonwater")
    if wat_len > 0:
        fh.write(" solvent")
    fh.write(" }")
    fh.write("""
setbox x centers
saveAmberParm x %s.parm7 %s.rst7
quit
EOF

tleap -s -f %s.sh.cmds | grep -v "+---" | grep -v "+Currently" > %s.sh.out
cat %s.sh.out

"""%(base,base,base,base,base))

    if parm.box is not None:
        fh.write("""
# Set the box dimensions
ChBox -X %.12f -Y %.12f -Z %.12f -al %.12f -bt %.12f -gm %.12f -c %s.rst7 -o %s.rst7.tmp; mv %s.rst7.tmp %s.rst7
"""%(parm.box[0],parm.box[1],parm.box[2],parm.box[3],parm.box[4],parm.box[5],base,base,base,base))

        if abs(parm.box[-1]-109.471219000000) < 1.e-4:
            fh.write("""
# Reset the ifbox flag
sed -i 's|0       0       1|0       0       2|' %s.parm7

"""%(base))



##############################################################
##############################################################
##############################################################
##############################################################
##############################################################


def write_makelatex():
    fh=file("makelatex.py","w")
    fh.write("""#!/usr/bin/env python2.7

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    from tianalysis import DataLoc as dloc
    from tianalysis import UsePkaUnits
    from tianalysis import make_latex_document
    from collections import defaultdict as ddict
    import os
    import glob
    
    #UsePkaUnits()
    
    tequil=0
    tmax=1.e+10
    odir="latex"
    
    D = ddict( lambda: ddict( lambda: ddict( list ) ) )
    for i in ["LIG"]:

        deqs=[]
        scs=[]
        reqs=[]
        for d in glob.glob("stddeq/*/results/data"):
            trial=d.split("/")[1]
            if os.path.exists(d):
                if len(glob.glob(os.path.join(d,"dvdl_*.dat"))) > 0:
                    deqs.append( dloc(trial,d,"",tequil,tmax) )
        for d in glob.glob("stdsc/*/results/data"):
            trial=d.split("/")[1]
            if os.path.exists(d):
                if len(glob.glob(os.path.join(d,"dvdl_*.dat"))) > 0:
                    scs.append( dloc(trial,d,"",tequil,tmax) )
        for d in glob.glob("stdreq/*/results/data"):
            trial=d.split("/")[1]
            if os.path.exists(d):
                if len(glob.glob(os.path.join(d,"dvdl_*.dat"))) > 0:
                    reqs.append( dloc(trial,d,"",tequil,tmax) )

        
        if len(deqs) > 0:
            D[i]["bio"]["deq"] = deqs
        if len(scs) > 0:
            D[i]["bio"]["sc"]  = scs
        if len(reqs) > 0:
            D[i]["bio"]["req"]  = reqs

        deqs=[]
        scs=[]
        reqs=[]
        if len(deqs) > 0:
            D[i]["ref"]["deq"] = deqs
        if len(scs) > 0:
            D[i]["ref"]["sc"]  = scs
        if len(reqs) > 0:
            D[i]["ref"]["req"]  = reqs


    import traceback
    try:
        make_latex_document( odir, D, methods=["TI","TI3","BAR","MBAR"] )
    except Exception, err:
        try:
            print "\n\n\n\n\n\n"
            print traceback.format_exc()
            print "***************************************************"
            print "Could not use methods TI,TI3,BAR,MBAR"
            print "Retrying with methods TI,TI3,BAR"
            print "***************************************************"
            make_latex_document( odir, D, methods=["TI","TI3","BAR"] )
        except:
            print "\n\n\n\n\n\n"
            print traceback.format_exc()
            print "***************************************************"
            print "Could not use methods TI,TI3,BAR"
            print "Retrying with methods TI,TI3"
            print "***************************************************"
            make_latex_document( odir, D, methods=["TI","TI3"] )

""")
    
    import os
    import stat
    st = os.stat("makelatex.py")
    os.chmod("makelatex.py", st.st_mode | stat.S_IEXEC)

    

##############################################################
##############################################################
##############################################################
##############################################################
##############################################################





def RefitCharges( p, timask, scmask, equivH=False, equivQ=False ):
    import subprocess
    
    atoms = [ p.atoms[i].idx for i in parmed.amber.mask.AmberMask( p, timask ).Selected() ]
    disappearing = [ p.atoms[i].idx for i in parmed.amber.mask.AmberMask( p, scmask ).Selected() ]

    atoms.sort()
    disappearing.sort()
    atomq = [ p.atoms[i].charge for i in atoms ]
    atns  = [ p.atoms[i].atomic_number for i in atoms ]

    umask = [0]*len(atoms)
    for i in range(len(atoms)):
        if atoms[i] in disappearing:
            umask[i] = -1
    for i in range(len(atoms)):
        if atoms[i] in disappearing:
            continue
        if umask[i] > 0:
            continue
        qi = atomq[i]
        for j in range(i+1,len(atoms)):
            qj = atomq[j]
            if equivQ:
                if abs(qj-qi) < 1.e-5:
                    umask[j]=i+1
            elif equivH:
                if abs(qj-qi) < 1.e-5 and \
                   atns[i]==1 and \
                   atns[j]==1:
                    umask[j]=i+1
        
    q0=0
    fh = file("mmresp.inp","w")
    fh.write("%i\n\n"%(len(atoms)))
    for i,msk in zip(atoms,umask):
        a = p.atoms[i]
        n = a.atomic_number
        if i in disappearing:
            n = -n
        fh.write("%3i %14.8f %14.8f %14.8f   %14.8f  %3i\n"%(n,a.xx,a.xy,a.xz,a.charge,msk))
        q0 += a.charge
    
    fh.close()
    subprocess.call("mmresp mmresp.inp > mmresp.out",shell=True)
    fh = file("mmresp.out","r")
    q=[]
    for line in fh:
        line = line.strip()
        if len(line) > 0:
            q.append( float(line) )
    fh.close()
    

    print "qin "," ".join(["%10.6f"%(qa) for qa in atomq])
    print "qout"," ".join(["%10.6f"%(qa) for qa in q])

    if len(q) != len(atomq):
        raise Exception("Failed to refit charges. (size mismatch)")
    
    
    return q

    
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################




class computer_info(object):
    def __init__(self):
        self.num_cpu_per_node = 24
        self.num_gpu_per_node =  4
        self.min_cpu_nodes = -1
        self.max_cpu_nodes = -1
        self.num_days  =  2
        self.nstlim = 1000000
        self.ntpr = 20000
        self.ntwx = 20000
        self.numexchg = 0
        self.xstream = False
        self.cpu_partition = "main"
        self.gpu_partition = "gpu"
        self.launch = "srun --mpi=pmi2"
        self.steps = 100
        self.steps_per_slurm = 10
        self.qos=""
        self.account=""
        self.exclude=""
        self.constraint=""
        self.parallel_cpu=True
        self.parallel_gpu=False
        self.serial_cpu=False
        self.serial_gpu=False

def find_num_nodes_and_cores_gas( njobs, cpu_per_node, min_num_nodes, max_num_nodes ):
    max_nodes=max_num_nodes
    if max_nodes < 1:
        max_nodes=njobs
    for ii in range(1,max_nodes):
        if 2*njobs < ii*cpu_per_node:
            num_nodes=ii
            if (2*njobs) % ii == 0:
                num_cores_per_node = (2*njobs)/ii
                break
    return num_nodes,num_cores_per_node

def find_num_nodes_and_cores( njobs, cpu_per_node, min_num_nodes, max_num_nodes ):
    if min_num_nodes < 1:
        min_num_nodes = njobs
    if max_num_nodes < 1:
        max_num_nodes = njobs
    orig_cpu_per_node = cpu_per_node
    while True:
        num_nodes = min_num_nodes
        while True:
            if ( num_nodes*cpu_per_node ) % njobs == 0:
                break
            else:
                num_nodes += 1
        if num_nodes > max_num_nodes:
            cpu_per_node -= 1
        else:
            break
        if cpu_per_node < orig_cpu_per_node * 0.75:
            cpu_per_node = orig_cpu_per_node
            num_nodes = njobs
            break
    return num_nodes,cpu_per_node


def min_num_nodes_and_cores( njobs, cpu_per_node ):
    orig_cpu_per_node = cpu_per_node
    num_nodes = 1
    cpu_per_node = 1
    while True:
        if num_nodes*cpu_per_node == njobs:
            break
        elif cpu_per_node < orig_cpu_per_node:
            cpu_per_node += 1
        else:
            num_nodes += 1
            cpu_per_node = 1
    return num_nodes,cpu_per_node
            

def write_stdti_mdin(fname,basename,clambda,lams,crgmask,noshakemask,timask1,timask2,scmask1,scmask2,compinfo,is_piti=False,is_gas=False):
    fh=open(fname,"w")

    icfe=1
    ifsc=0
    chngmask=1
    acyclic_remd=""
    if len(scmask1) != 0 or len(scmask2) != 0:
        ifsc=1

    is_piscti=False
    if is_piti and ifsc==1:
        is_piscti=True
        
    if is_piscti:
        icfe=0
        ifsc=0
        chngmask=0
        acyclic_remd="acyclic_remd= 0"

    lam="CLAMBDA"
    try:
        lam = "%.8f"%(clambda)
    except:
        pass

        
        
    fh.write("""Title
&cntrl
! STEPS ===========================
imin        = 0
maxcyc      = 0
ntmin       = 2
nstlim      = %i
numexchg    = %i
dt          = 0.001
! START/RESTART ===================
irest       = 1
ntx         = 5
ntxo        = 1
! ENERGY
noshakemask = '%s'
ntc         = 2
ntf         = 1
tol         = 0.00001
! PRINTING ========================
ntwx        = %i
ntwr        = %i
ntpr        = %i
ntwe        = 0
"""%(compinfo.nstlim,
         compinfo.numexchg,
         noshakemask,
         compinfo.ntwx,
         compinfo.ntpr,
         compinfo.ntpr))
    if not is_gas:
        fh.write("""! PERIODIC ========================
cut         = %.1f
ntb         = 1
iwrap       = 1
nscm        = 0
"""%(compinfo.cut))
    else:
        fh.write("""! PERIODIC ========================
cut         = %.1f
ntb         = 0
iwrap       = 0
nscm        = 1000
igb         = 6
fixcom      = 1
vrand       = 1000
"""%(compinfo.cut))
    fh.write("""! TEMPERATURE =====================
tempi       = 298.
temp0       = 298.
ntt         = 3
gamma_ln    = 1.
! PRESSURE ========================
ntp         = 0
! RESTRAINTS ======================
ntr         = 0
restraintmask = ''
restraint_wt  = 0.
nmropt      = %i
! MISC ============================
klambda     = 1
ig          = -1"""%(compinfo.nmropt))
    
    if is_piti and not is_piscti:
        fh.write("""
/
""")
    else:
        if compinfo.gti:
            fh.write("\ngti_add_sc   = 2")
            fh.write("\ngti_cut      = 1")
            fh.write("\ngti_chg_keep = 1")
            fh.write("\ngti_output   = 1")
            fh.write("\ngti_lam_sch  = 0")
            fh.write("\ngti_vdw_sc   = 1")
            fh.write("\ngti_ele_sc   = 1")
            fh.write("\nscalpha      = 0.2")
            fh.write("\nscbeta       = 50")
            if compinfo.numexchg > 0:
                if len(lams) % 2 == 1:
                    fh.write("\ngremd_acyc = 1")
                else:
                    fh.write("\ngremd_acyc = 0")
        fh.write("""
       icfe = %i
       ifsc = %i
    clambda = %s
"""%(icfe,ifsc,clambda))
        if not is_piti:
            fh.write("ifmbar=1, bar_intervall=1, mbar_states=%i\n"%(len(lams)))
            for ilam,lam in enumerate(lams):
                fh.write("mbar_lambda(%i)=%.8f\n"%(ilam+1,lams[ilam]))
            #fh.write("mbar_lambda=%s\n"%(",".join(["%.8f"%(lam) for lam in lams])))
        fh.write("""crgmask = '%s'
    timask1 = '%s'
    timask2 = '%s'
    scmask1 = '%s'
    scmask2 = '%s'
    %s !
/

&ewald
   chngmask = %i ! default 1; set to 0 when using piscti softcore
/

"""%(crgmask,timask1,timask2,scmask1,scmask2,acyclic_remd,chngmask))
    if compinfo.nmropt == 1:
        fh.write("""
&wt
type='DUMPFREQ', istep1=%i
&end
&wt
  type='END',
&end
DISANG=%s.disang
DUMPAVE=current/%s_%s.dumpave

"""%(compinfo.ntpr,basename,basename,clambda))
        
    fh.close()


def write_stdti_mdin_script(lams,basename,compinfo,
                            crgmask,
                            noshakemask,
                            timask1,timask2,
                            scmask1,scmask2,
                            is_piti=False,
                            is_gas=False):
    import os.path

    ntmin=1
    ifsc=0
    if len(scmask1) != 0 or len(scmask2) != 0:
        ifsc=1
        ntmin=2

    piti=0
    if is_piti:
        piti=1
        ntmin=1
        
    is_piscti=False
    if is_piti and ifsc==1:
        is_piscti=True
        
    # if len(lams) > 1:
    #     second_lam = 1.0 / ( len(lams) - 1 )
    #     #mbarstr = "ifmbar=1, bar_intervall=1, bar_l_min=0.0, bar_l_max=1.0, bar_l_incr=%.15f"%(second_lam)
    #     #mbarstr = "ifmbar=1, bar_intervall=1, mbar_states=%i, mbar_lambda=%s"%(len(lams),",".join(["%.8f"%(lam) for lam in lams]))
    #     mbarstr=""

    # else:
    #     mbarstr = "ifmbar=0"

    # if is_piti and not is_piscti:
    #     mbarstr = "ifmbar=0"

    dirname = os.path.join(basename,"template") 
    fname = os.path.join(dirname,"make.mdin.sh")
    template = os.path.join(dirname,"TEMPLATE.mdin")
    write_stdti_mdin(template,basename,"CLAMBDA",lams,crgmask,noshakemask,timask1,timask2,scmask1,scmask2,compinfo,is_piti=is_piti,is_gas=is_gas)
    fh = open(fname,"w")
    fh.write("#!/bin/bash\n\n")
    fh.write("lams=( %s )\n\n"%( " ".join( ["%.8f"%(lam) for lam in lams] ) ))
    fh.write("prefix=\"%s\"\n\n"%(basename))
    fh.write("mdin=TEMPLATE.mdin\n\n")
    piscti=0
    if is_piscti:
        piscti=1
    fh.write("""

#
# if you set cpu_equil=1, then gti options are removed from equilibration mdin files
#
cpu_equil=0

#
# is this softcore?
#
ifsc="%i"

#
# is this parameter interpolated?
#
piti="%i"

#
# is this a parameter-interpolated softcore calculation?
#

piscti="%i"
if [ "${piscti}" == "1" ]; then
    comment="\!"
fi

noshake=$(grep noshakemask TEMPLATE.mdin | sed -e "s/^.*= *//" -e "s/'//g")


#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg ${mdin} | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if this is replica exchange AND we are doing piscti, then do acyclic exchange loops
#

acyclic=0
if [ "${piscti}" == "1" -a "${rem}" == "3" ]; then
   acyclic=1
fi

    """%(ifsc,piti,piscti))

    
    fh.write("""

if [ ! -d inputs ]; then
   mkdir inputs
fi

for lam in ${lams[@]}; do   
    base=${prefix}_${lam}
    init=inputs/${base}_initial
    rest=inputs/${base}_restart
    anal=inputs/${base}_analyze
    mini=inputs/${base}_minimiz
    heat=inputs/${base}_heating
    pre1=inputs/${base}_press01
    pre2=inputs/${base}_press02
    pre3=inputs/${base}_press03
    pre4=inputs/${base}_press04
    pre5=inputs/${base}_press05
    pre6=inputs/${base}_press06
    pre7=inputs/${base}_press07
    pre8=inputs/${base}_press08


    sed -e "s/CLAMBDA/${lam}/" \\
        -e "s/irest *= *[0-9]*/irest = 1/" \\
        -e "s/ntx *= *[0-9]*/ntx = 5/" \\
        -e "s/acyclic_remd *= *[0-9]*/acyclic_remd = ${acyclic}/" \\
        -e "s/clambda/${comment}clambda/" \\
        -e "s/ifsc/${comment}ifsc/" \\
        -e "s/icfe/${comment}icfe/" \\
        -e "s/timask/${comment}timask/" \\
        -e "s/scmask/${comment}scmask/" \\
        -e "s/crgmask/${comment}crgmask/" \\
         ${mdin} > ${rest}.mdin

    sed -e "s/irest *= *[0-9]*/irest = 0/" \\
        -e "s/ntx *= *[0-9]*/ntx = 1/" \\
         ${rest}.mdin > ${init}.mdin

    sed -e "s/imin *= *[0-9]*/imin = 1/" \\
	-e "s/maxcyc *= *[0-9]*/maxcyc = 1000/" \\
	-e "s/ntc *= *[0-9]*/ntc = 1/" \\
	-e "s/iwrap *= *[0-9]*/iwrap = 0/" \\
        -e "s/ntwx *= *[0-9]*/ntwx = 500/" \\
	-e "s/ntwr *= *[0-9]*/ntwr = 500/" \\
	-e "s/ntpr *= *[0-9]*/ntpr = 500/" \\
	-e "s/ntr *= *[0-9]*/ntr = 1/" \\
	-e "s/restraintmask *=.*/restraintmask = '!:WAT \\& !@H= \& !${noshake}'/" \\
	-e "s/restraint_wt *= *[0-9\\.]*/restraint_wt = 100./" \\
        -e "s/acyclic_remd *= *[0-9]*/acyclic_remd = 0/" \\
        -e 's/gremd_acyc/!gremd_acyc/' \\
	-e "s/numexchg *= *[0-9]*/numexchg = 0/" \\
        -e "s/ifmbar *= *[0-9]*/ifmbar = 0/" \\
        ${init}.mdin > ${mini}.mdin

    if [ "${cpu_equil}" == "1" ]; then
       sed -i -e 's/gti/!gti/' ${mini}.mdin
    fi

    sed -e "s/imin *= *[0-9]*/imin = 0/" \\
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \\
	-e "s/nstlim *= *[0-9]*/nstlim = 100000/" \\
	-e "s/tempi *= *[0-9\\.]*/tempi = 100./" \\
	-e "s/temp0 *= *[0-9\\.]*/temp0 = 298./" \\
	-e "s/ntc *= *[0-9]*/ntc = 2/" \\
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \\
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \\
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \\
	-e "s/nmropt *= *[0-9]*/nmropt = 1/" \\
	${mini}.mdin > ${heat}.mdin

    
    cat <<EOF >> ${heat}.mdin

 &wt
     TYPE="TEMP0", istep1=0, istep2=100000,
     value1=100., value2=298.,
 /
 &wt
     TYPE="END",
 /

EOF
    
    if [ "${rem}" == "0" ]; then

    sed -e "s/imin *= *[0-9]*/imin = 0/" \\
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \\
	-e "s/nstlim *= *[0-9]*/nstlim = 2000/" \\
	-e "s/irest *= *[0-9]*/irest = 1/" \\
	-e "s/ntx *= *[0-9]*/ntx = 5/" \\
	-e "s/tempi *= *[0-9\\.]*/tempi = 298./" \\
	-e "s/temp0 *= *[0-9\\.]*/temp0 = 298./" \\
	-e "s/ntc *= *[0-9]*/ntc = 2/" \\
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \\
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \\
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \\
	-e "s/ntb *= *[0-9]*/ntb = 2/" \\
	-e "s/ntp *= *[0-9]*/ntp = 1/" \\
	${mini}.mdin > ${pre1}.mdin

    else


    sed -e "s/imin *= *[0-9]*/imin = 0/" \\
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \\
	-e "s/nstlim *= *[0-9]*/nstlim = 2000/" \\
	-e "s/irest *= *[0-9]*/irest = 1/" \\
	-e "s/ntx *= *[0-9]*/ntx = 5/" \\
	-e "s/tempi *= *[0-9\\.]*/tempi = 298./" \\
	-e "s/temp0 *= *[0-9\\.]*/temp0 = 298./" \\
	-e "s/ntc *= *[0-9]*/ntc = 2/" \\
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \\
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \\
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \\
	${mini}.mdin > ${pre1}.mdin

    fi

    sed -e "s/nstlim *= *[0-9\\.]*/nstlim = 18000/" \\
	${pre1}.mdin > ${pre2}.mdin

    sed -e "s/nstlim *= *[0-9\\.]*/nstlim = 25000/" \\
	${pre1}.mdin > ${pre3}.mdin

    sed -e "s/nstlim *= *[0-9\\.]*/nstlim = 55000/" \\
	${pre1}.mdin > ${pre4}.mdin

    sed -e "s/nstlim *= *[0-9\\.]*/nstlim = 100000/" \\
        -e "s/restraint_wt *= *[0-9\\.]*/restraint_wt = 10./" \\
	-e "s/irest *= *[0-9]*/irest = 1/" \\
	-e "s/ntx *= *[0-9]*/ntx = 5/" \\
	${pre4}.mdin > ${pre5}.mdin

    sed -e "s/restraintmask *=.*/restraintmask = '@CA,N,C \& !${noshake}'/" \\
	${pre5}.mdin > ${pre6}.mdin

    sed -e "s/restraint_wt *= *[0-9\\.]*/restraint_wt = 1./" \\
	${pre6}.mdin > ${pre7}.mdin

    sed -e "s/restraint_wt *= *[0-9\\.]*/restraint_wt = 0.1/" \\
	${pre7}.mdin > ${pre8}.mdin

    sed -e "s/${comment}clambda/clambda/" \\
        -e "s/${comment}ifsc/ifsc/" \\
        -e "s/${comment}icfe/icfe/" \\
        -e "s/${comment}timask/timask/" \\
        -e "s/${comment}scmask/scmask/" \\
        -e "s/${comment}crgmask/crgmask/" \\
        -e "s/imin *= *[0-9]*/imin = 6/" \\
        -e "s/ntpr *= *[0-9]*/ntpr = 1/" \\
        -e "s/ntwx *= *[0-9]*/ntwx = 0/" \\
        -e "s/ntwr *= *[0-9]*/ntwr = 0/" \\
        -e "s/icfe *= *[0-9]*/icfe = 1/" \\
        -e "s/ifsc *= *[0-9]*/ifsc = ${ifsc}/" \\
        -e "s/numexchg *= *[0-9]*/numexchg = 0/" \\
        -e "s/acyclic_remd/\\!acyclic_remd/" \\
        -e "s/chngmask *= *[0-9]*/chngmask = 1/" \\
         ${rest}.mdin > ${anal}.mdin
done

""")
    
    fh.write("""
truncate -s0 inputs/${prefix}_initial.groupfile
truncate -s0 inputs/${prefix}_restart.groupfile
truncate -s0 inputs/${prefix}_analyze.groupfile
truncate -s0 inputs/${prefix}_minimiz.groupfile
truncate -s0 inputs/${prefix}_heating.groupfile
truncate -s0 inputs/${prefix}_press01.groupfile
truncate -s0 inputs/${prefix}_press02.groupfile
truncate -s0 inputs/${prefix}_press03.groupfile
truncate -s0 inputs/${prefix}_press04.groupfile
truncate -s0 inputs/${prefix}_press05.groupfile
truncate -s0 inputs/${prefix}_press06.groupfile
truncate -s0 inputs/${prefix}_press07.groupfile
truncate -s0 inputs/${prefix}_press08.groupfile

for lam in ${lams[@]}; do
    base=${prefix}_${lam}
    init=${base}_initial
    rest=${base}_restart
    anal=${base}_analyze
    mini=${base}_minimiz
    heat=${base}_heating
    pre1=${base}_press01
    pre2=${base}_press02
    pre3=${base}_press03
    pre4=${base}_press04
    pre5=${base}_press05
    pre6=${base}_press06
    pre7=${base}_press07
    pre8=${base}_press08

    parm=${prefix}.parm7
    if [ "${piti}" == "1" ]; then
       parm=${base}.parm7
    fi

    echo "-O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o equilib/${mini}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_minimiz.groupfile
    echo "-O -c equilib/${mini}.rst7 -p ${parm} -i inputs/${heat}.mdin -o equilib/${heat}.mdout -r equilib/${heat}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_heating.groupfile
    echo "-O -c equilib/${heat}.rst7 -p ${parm} -i inputs/${pre1}.mdin -o equilib/${pre1}.mdout -r equilib/${pre1}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_press01.groupfile
    echo "-O -c equilib/${pre1}.rst7 -p ${parm} -i inputs/${pre2}.mdin -o equilib/${pre2}.mdout -r equilib/${pre2}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre1}.rst7" >> inputs/${prefix}_press02.groupfile
    echo "-O -c equilib/${pre2}.rst7 -p ${parm} -i inputs/${pre3}.mdin -o equilib/${pre3}.mdout -r equilib/${pre3}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre2}.rst7" >> inputs/${prefix}_press03.groupfile
    echo "-O -c equilib/${pre3}.rst7 -p ${parm} -i inputs/${pre4}.mdin -o equilib/${pre4}.mdout -r equilib/${pre4}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre3}.rst7" >> inputs/${prefix}_press04.groupfile
    echo "-O -c equilib/${pre4}.rst7 -p ${parm} -i inputs/${pre5}.mdin -o equilib/${pre5}.mdout -r equilib/${pre5}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press05.groupfile
    echo "-O -c equilib/${pre5}.rst7 -p ${parm} -i inputs/${pre6}.mdin -o equilib/${pre6}.mdout -r equilib/${pre6}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press06.groupfile
    echo "-O -c equilib/${pre6}.rst7 -p ${parm} -i inputs/${pre7}.mdin -o equilib/${pre7}.mdout -r equilib/${pre7}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press07.groupfile
    echo "-O -c equilib/${pre7}.rst7 -p ${parm} -i inputs/${pre8}.mdin -o equilib/${pre8}.mdout -r equilib/${pre8}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press08.groupfile


    
    echo "-O -c inputs/${init}.rst7 -p ${parm} -i inputs/${init}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile" >> inputs/${prefix}_initial.groupfile
    echo "-O -c current/${rest}.rst7 -p ${parm} -i inputs/${rest}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile" >> inputs/${prefix}_restart.groupfile
    echo "-O -c ${base}.rst7 -p ../../${parm} -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -l ${anal}.logfile -y ${base}.nc " >> inputs/${prefix}_analyze.groupfile

    if [ ! -e "inputs/${init}.rst7" ]; then
       if [ -e "${base}.rst7" ]; then
          cd inputs
          ln -s ../${base}.rst7 ${init}.rst7
          cd ../
       elif [ -e "${prefix}.rst7" ]; then
          cd inputs
          ln -s ../${prefix}.rst7 ${init}.rst7
          cd ../
       fi
    fi
done
""")


    fh.close()

    
    # make script executable by user
    import os
    import stat
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)

    mergescript = os.path.join(dirname,"analysis2results.py")
    fh=file(mergescript,"w")
    fh.write("""#!/usr/bin/env python2.7
import os
from collections import defaultdict as ddict
prefix="%s"
lams = [ %s ]
merge_gaps=False
"""%(basename,",".join( ["\"%.8f\""%(lam) for lam in lams ] )))
    fh.write("""
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
            fh.write("%12.1f %s\\n"%(t,dvdl_data[t][lam]))
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
                fh.write("%12.1f %s\\n"%(t,efep_data[t][tlam][plam]))
            fh.close()

""")
    fh.close()
    st = os.stat(mergescript)
    os.chmod(mergescript, st.st_mode | stat.S_IEXEC)

    write_stdti_production_slurm_script(lams,basename,compinfo,is_piti=is_piti,is_piscti=is_piscti,is_gas=is_gas)
    write_stdti_heating_slurm_script(lams,basename,compinfo,is_piti=is_piti,is_piscti=is_piscti,is_gas=is_gas)
    if (not is_piti) and (not is_piscti):
        write_stdti_analysis_slurm_script(lams,basename,compinfo,is_piti=is_piti,is_piscti=is_piscti,is_gas=is_gas)
    else:
        write_piti_analysis_slurm_script(lams,basename,compinfo)




def write_stdti_production_slurm_script(lams,basename,compinfo,is_piti=False,is_piscti=False,is_gas=False):
    import os.path
    fname = os.path.join( os.path.join(basename,"template"),"production.slurm.sh")
    lamsstr="lams=( %s )\n"%( " ".join(["%.8f"%(lam) for lam in lams]) )
    lamsstr+="REVERSE=false\n"
    lamsstr+="if [ \"${REVERSE}\" = true ]; then lams=( %s ); fi\n"%( " ".join(["%.8f"%(lam) for lam in reversed(lams)]) )

    pmemd="pmemd"
    pmemdCUDA="pmemd.cuda"
    if is_piscti:
        pmemd="pmemdPI"
        pmemdCUDA="pmemdPI.cuda"
    elif is_piti:
        pmemd="pmemd"
        pmemdCUDA="pmemd.cuda"
    
    fh = open(fname,"w")
    fh.write("""#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi


    ###################################################################
""")
    
    fh.write("""
    local first_job_gets_to_step=%i
    local final_job_gets_to_step=%i
    local steps_per_job=%i

    CHECKNAN="F"
"""%(compinfo.steps_per_slurm,compinfo.steps,compinfo.steps_per_slurm))


    pc="#"
    pg="#"
    sc="#"
    sg="#"
    if compinfo.parallel_cpu:
        pc=""
    elif compinfo.parallel_gpu:
        pg=""
    elif compinfo.serial_cpu:
        sc=""
    elif compinfo.serial_gpu:
        sg=""
        
    fh.write("""
    %slocal write_template=write_parallel_cpu_template
    %slocal write_template=write_parallel_gpu_template
    # THE SERIAL-CASES CAN ONLY BE USED FOR NON-REMD SIMULATIONS
    %slocal write_template=write_serial_cpu_template
    %slocal write_template=write_serial_gpu_template

    ###################################################################


"""%(pc,pg,sc,sg))

    fh.write("    local template=production.slurm")

    fh.write("""
    
    local JARR=($( seq ${first_job_gets_to_step} ${steps_per_job} ${final_job_gets_to_step} ))
    if [ "${JARR[${#JARR[@]}-1]}" != "${final_job_gets_to_step}" ]; then JARR+=(${final_job_gets_to_step}); fi

    ${write_template} "${template}" "${CHECKNAN}" "${JARR[@]}"

    local lastid=$(squeue -h -o "%F %o" | grep ${PWD} | sed 's/ .*//' | sort -rn | head -n 1)
    if [ "${lastid}" != "" ]; then
       cmdstr="sbatch --dependency=afterany:${lastid} ${template}"
       line=$( sbatch --dependency=afterany:${lastid} ${template} )
    else
       cmdstr="sbatch ${template}"
       line=$( sbatch ${template} )
    fi
    echo "${cmdstr}"
    echo "${line}"
}


""")

    num_nodes,num_cores_per_node = find_num_nodes_and_cores( len(lams), compinfo.num_cpu_per_node, compinfo.min_cpu_nodes, compinfo.max_cpu_nodes )
    if is_gas:
        num_nodes,num_cores_per_node = find_num_nodes_and_cores_gas( len(lams), compinfo.num_cpu_per_node, compinfo.min_cpu_nodes, compinfo.max_cpu_nodes )


    fh.write("""
##############################################################################

write_parallel_cpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=%s
#SBATCH --nodes=%i
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
##SBATCH --exclude=cuda[001-008]
"""%(compinfo.cpu_partition,
     num_nodes,
     num_cores_per_node,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))


    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"

export EXE=${AMBERHOME}/bin/%s.MPI
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_parallel_body "${fname}"
}
"""%(compinfo.launch,
     pmemd,
     basename,
     lamsstr))

    pmemdMPI="%s.MPI"%(pmemd)
    ncpu = compinfo.num_cpu_per_node
    if is_gas:
        ncpu=2
        if is_piti:
            ncpu=1
            pmemdMPI=pmemd
    fh.write("""
##############################################################################
write_serial_cpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
##SBATCH --exclude=cuda[001-008]
"""%(compinfo.cpu_partition,
     ncpu,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        
    

    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_serial_body "${fname}"
}
"""%(compinfo.launch,
     pmemdMPI,
     basename,
     lamsstr))

    if compinfo.oversub:
        num_nodes = 1
        num_cores_per_node = len(lams)
    else:
        num_nodes,num_cores_per_node = min_num_nodes_and_cores( len(lams), compinfo.num_gpu_per_node )

        
    if compinfo.xstream:
        num_nodes,num_cores_per_node = min_num_nodes_and_cores( len(lams), 16 )
        fh.write("""
##############################################################################

write_parallel_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=normal
#
#SBATCH --nodes=%i
#SBATCH --ntasks=%i
#SBATCH --mincpus=%i
#SBATCH --gres=gpu:%i
#SBATCH --cpus-per-task=1
#
##SBATCH --share
#SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s.MPI
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_parallel_body "${fname}"
}
"""%(num_nodes,
     len(lams),
     num_cores_per_node,
     num_cores_per_node,
     compinfo.num_days,
     compinfo.launch,
     pmemdCUDA,
     basename,
     lamsstr))


        fh.write("""
##############################################################################
write_serial_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=normal
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mincpus=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#
##SBATCH --share
#SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
"""%(compinfo.num_days))
        
        fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"

export EXE=${AMBERHOME}/bin/%s
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_serial_body "${fname}"
}
"""%(compinfo.launch,
     pmemdCUDA,
     basename,
     lamsstr))

        
    else:
        fh.write("""
##############################################################################

write_parallel_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=%s
#
#SBATCH --nodes=%i
#SBATCH --ntasks-per-node=%i
##SBATCH --mincpus=%i
#SBATCH --gres=gpu:%i
#SBATCH --cpus-per-task=1
#
#SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
"""%(compinfo.gpu_partition,
     num_nodes,
     num_cores_per_node,
     num_cores_per_node,
     compinfo.num_gpu_per_node,
     compinfo.num_days))
    
        if len(compinfo.qos) > 0:
            fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
        if len(compinfo.account) > 0:
            fh.write("#SBATCH --account=%s\n"%(compinfo.account))
        if len(compinfo.exclude) > 0:
            fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
        if len(compinfo.constraint) > 0:
            fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        
        
        fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s.MPI
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_parallel_body "${fname}"
}
"""%(compinfo.launch,
     pmemdCUDA,
     basename,
     lamsstr))


        fh.write("""
##############################################################################
write_serial_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=%s
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mincpus=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#
#SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
"""%(compinfo.gpu_partition,
     compinfo.num_days))
    
        if len(compinfo.qos) > 0:
            fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
        if len(compinfo.account) > 0:
            fh.write("#SBATCH --account=%s\n"%(compinfo.account))
        if len(compinfo.exclude) > 0:
            fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
        if len(compinfo.constraint) > 0:
            fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        
        fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"

export EXE=${AMBERHOME}/bin/%s
prefix="%s"
%s
EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF

write_serial_body "${fname}"
}
"""%(compinfo.launch,
     pmemdCUDA,
     basename,
     lamsstr))

        



    fh.write("""
##############################################################################
write_parallel_body() {
   local fname="$1"
   cat << 'EOF' >> ${fname}
#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg inputs/${prefix}_${lams[0]}_initial.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#


if [ ! -d current ]; then
   mkdir current
fi


#
# each step consists of running all windows for some length of time
#

for istep in $(seq ${laststep}); do

   istep=$(printf "%06i" ${istep})

   #
   # if we've already run this step, then skip to the next step
   #

   if [ -d production/${istep} ]; then
      continue
   fi

   #
   # this infinite-loop will exit if all simulations run correctly
   #

   while true; do

      #
      # run all simulations at the same time using a groupfile.
      # there are different mdin and rst7 files depending on if we are
      # starting or restarting the simulation
      #

      if [ -e "current/${prefix}_${lams[0]}_restart.rst7" ]; then

         ${LAUNCH} ${EXE} -rem ${rem} -ng ${#lams[@]} -groupfile inputs/${prefix}_restart.groupfile

      else 

         echo "Missing current/${prefix}_${lams[0]}_restart.rst7"
         echo "You should first run equilibration.slurm.sh"
         exit 1

#         if [ ! -d equilib ]; then
#            mkdir equilib
#         fi
#         ${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_minimiz.groupfile
#         for lam in ${lams[@]}; do
#             initial=${prefix}_${lam}_initial
#             minimiz=${prefix}_${lam}_minimiz
#             if [ -L inputs/${initial}.rst7 ]; then
#                  rm -f inputs/${initial}.rst7
#             fi
#             cp equilib/${minimiz}.rst7 inputs/${initial}.rst7
#         done
#         ${LAUNCH} ${EXE} -rem ${rem} -ng ${#lams[@]} -groupfile inputs/${prefix}_initial.groupfile

      fi

      ok=1
      for lambda in ${lams[@]}; do
           base="${prefix}_${lambda}"
           if [ ! -e current/${base}.mdout ]; then
              ok=0
           fi
      done
      if [ "${ok}" == "0" ]; then
         echo "There was an unrecoverable error. A mdout file was not produced"
         exit 1
      else
	  incomplete=($(grep -L '5.  TIMINGS' current/${prefix}_*.mdout))
	  if [ "${#incomplete[@]}" -gt 0 ]; then
	      echo "The jobs did not complete for an unknown reason"
	      exit 1
	  fi
      fi

      #
      # if there aren't any not-a-numbers in the mdout files,
      # then the simulations must have run correctly,
      # so copy the output files to a subdirectory and exit
      # the infinite loop
      #

      if [ "$(grep -l NaN current/*.mdout)" == "" -o "${CHECKNAN}" == "F" ]; then

         mkdir -p production/${istep}
         for lambda in ${lams[@]}; do
            base="${prefix}_${lambda}"

            #
            # save the coordinates to restart the next simulation
            #

            cp current/${base}.rst7 current/${base}_restart.rst7

            # 
            # now copy all the outputs to a subdirectory
            #

            mv current/${base}.nc current/${base}.rst7 current/${base}.mdout current/${base}.mdinfo production/${istep}/
            for dumpave in current/*.dumpave; do
               if [ -e "${dumpave}" ]; then
                  mv ${dumpave} production/${istep}/$(basename ${dumpave})
               fi
            done
            if [ -e "current/${base}.logfile" ]; then
               mv current/${base}.logfile production/${istep}/
            fi
         done
         if [ -e "rem.log" ]; then
            mv rem.log production/${istep}/
         fi

         # cd production/${istep}
         # ln -s ../../*.parm7 .
         # cd ../../

         break # exit the infinite loop

     fi
  done
done

EOF
    
}
""")
    fh.write("""
##############################################################################
write_serial_body() {
local fname="$1"
cat << 'EOF' >> ${fname}

#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg inputs/${prefix}_${lams[0]}_initial.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN SERIAL-MODE"
   exit 1
fi


if [ ! -d current ]; then
    mkdir current
fi


#
# each step consists of running all windows for some length of time
#

for istep in $(seq ${laststep}); do

   istep=$(printf "%06i" ${istep})

   #
   # if we've already run this step, then skip to the next step
   #

   if [ -d production/${istep} ]; then
      continue
   fi


   #
   # run a simulation for each window
   #

   for lambda in ${lams[@]}; do
       base="${prefix}_${lambda}"
       rest="${base}_restart"
       init="${base}_initial"
       mini="${base}_minimiz"
""")
    if not is_piti:
        fh.write("       parm=\"${prefix}.parm7\"\n")
    else:
        fh.write("       parm=\"${base}.parm7\"\n")        
    fh.write("""

       #
       # do we need to run this window?
       #

       needed=1
       if [ -e "current/${base}.mdout" ]; then
          if [ "$(grep -LE 'Run   done at|Final Performance Info' current/${base}.mdout)" == "" ]; then
             needed=0
          fi
       fi

       #
       # if we need this window, then run it until it until it terminates normally
       #

       if [ "${needed}" == "1" ]; then
          while true; do

             #
             # are we starting or restarting a simulation?
             # we need to use different mdin and rst7 files depending on the answer
             #

             if [ -e "current/${rest}.rst7" ]; then

                ${LAUNCH} ${EXE} -O -c current/${rest}.rst7 -p ${parm} -i inputs/${rest}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo

             else 

                echo "Missing current/${rest}.rst7"
                echo "You should first run equilibration.slurm.sh"
                exit 1

#                if [ ! -d equilib ]; then
#                    mkdir equilib
#                fi
#                ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o current/${base}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile
#                if [ -L inputs/${init}.rst7 ]; then
#                    rm -f inputs/${init}.rst7
#                fi
#                cp equilib/${mini}.rst7 inputs/${init}.rst7
#                ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${init}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo

             fi

             #
             # if we don't get not-a-numbers, then it must've worked, so leave this loop
             #

             if [ ! -e current/${base}.mdout ]; then
                echo "There was an unrecoverable error. A mdout file was not produced"
                exit 1
             else
                incomplete=($(grep -L '5.  TIMINGS' current/${base}.mdout))
                if [ "${#incomplete[@]}" -gt 0 ]; then
                   echo "The jobs did not complete for an unknown reason"
                   exit 1
                fi
             fi 

             if [ "$(grep -l NaN current/${base}.mdout)" == "" -o "${CHECKNAN}" == "F" ]; then
                break
             fi
          done
       fi
   done

   #
   # we've run all the windows for this step, so store the outputs in a subdirectory
   #

   mkdir -p production/${istep}
   for lambda in ${lams[@]}; do
      base="${prefix}_${lambda}"

      #
      # save these coordinates to restart the next step
      # 

      cp current/${base}.rst7 current/${base}_restart.rst7

      # 
      # copy the outputs
      #

      mv current/${base}.nc current/${base}.rst7 current/${base}.mdout current/${base}.mdinfo production/${istep}/
      for dumpave in current/*.dumpave; do
          if [ -e "${dumpave}" ]; then
              mv ${dumpave} production/${istep}/$(basename ${dumpave})
          fi
      done
      if [ -e "current/${base}.logfile" ]; then
         mv current/${base}.logfile production/${istep}/
      fi

   done
   if [ -e "rem.log" ]; then
      mv rem.log production/${istep}/
   fi

   if [ -e logfile ]; then
      rm -f logfile
   fi

   # cd production/${istep}
   # ln -s ../../*.parm7 .
   # cd ../../

done

EOF
    
}


# ---------------------------
# call to the main function
# ---------------------------

main

""")

    # make script executable by user
    import os
    import stat
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)
    




def write_stdti_analysis_slurm_script(lams,basename,compinfo,is_piti=False,is_piscti=False,is_gas=False):
    import os.path
    fname = os.path.join( os.path.join(basename,"template"),"analyze.slurm.sh")
    fh = open(fname,"w")
    fh.write("""#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi

    local base="%s" 
    local lams=( %s )"""%(basename," ".join("%.8f"%(lam) for lam in lams)))
    fh.write("""

    #
    # if we don't have any completed production, then there's nothing to analyze
    #

    if [ ! -d production ]; then 
       echo "Nothing to analyze because you don't have a production directory yet"
       exit
    fi

    #
    # the analysis will use a python utility to transform the data in the
    # re-analyzed mdout files to a series of .dat files
    #

    if [ ! -e "production/stdti_step2dats.py" ]; then
       write_python_script
    fi

    #
    # the production dir contains many subdirs that are 0-padded integers.
    # what is the largest number that we can find?
    #

    local last_step=$(for f in $(ls production | grep -v '\.' | tail -n 1); do bc -l <<< $f; done)

    #
    # if we couldn't find a valid subdirectory, then exit now
    #

    if [ "${last_step}" == "" ]; then 
       echo "Nothing to analyze because I couldn't find a valid subdirectory name in production/"
       exit
    fi


    #
    # collect the existing results, if any
    #

    echo "Collecting existing results before submitting new analysis..."
    python2.7 analysis2results.py

    echo ""
    echo ""
    echo "Searching for un-analyzed production directories..."

    local script="analysis.slurm"
    local jarr=()
    for step in $(seq ${last_step}); do

       #
       # this is the zero-padded name
       # 

       local step_name=$(printf "%06i" ${step})

       #
       # if the subdir does not exist, then skip it
       #

       if [ ! -d "production/${step_name}" ]; then
          echo "skipping ${step_name} because dir does not exist"
          continue
       fi


       #
       # do we actually have to analyze this subdir?
       #


       #
       # if it doesn't have mdout files, then I think I already deleted it to save disk space
       #

       local ok=1
       local cnt=0
       for lam in ${lams[@]}; do
           if [ -e "production/${step_name}/${base}_${lam}.mdout" ]; then
               cnt=$(( ${cnt} + 1 ))
           else
               ok=0
           fi
       done
       if [ "${cnt}" == "0" ]; then
          echo "skipping ${step_name} because it was probably already deleted from this machine"
          continue
       elif [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it has some, but not all, of the mdout files (something wrong here?)"
          continue
       fi


       ok=1
       for lam in ${lams[@]}; do
           if [ ! -e "production/${step_name}/${base}_${lam}.nc" ]; then
               ok=0
               echo "Missing production/${step_name}/${base}_${lam}.nc"
           fi
       done
       if [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it is missing 1-or-more nc files"
          continue
       fi


       #
       # do we have the mbar trace file generated from the analysis
       #

       local testfile=$(printf "production/${step_name}/efep_%.8f_%.8f.dat" 1. 1.)
       if [ -e ${testfile} ]; then
          echo "skipping ${step_name} because it has already been analyzed"
          continue
       fi

""")

    if not is_piti:
        fh.write("""
       cd production/${step_name}
       python2.7 ../stdti_step2dats.py ${base}_*.mdout
       cd ../../
""")
    
    if is_piti:
        fh.write("""

       local queued=$(squeue -ro "%i %o" | grep ${PWD} | grep ${script} | sed -e 's/ .*/ /' -e 's/.*_/_/' -e 's/ /_/' | grep "_${step}_")
       if [ -e "${queued}" ]; then
          echo "skipping ${step_name} because it is already queued"
          continue
       fi 

       jarr+=(${step})

""")

    fh.write("""

    done

""")

    if not is_piti:
        fh.write("""
    python2.7 analysis2results.py

""")

    if is_piti:
        fh.write("""

    if [ "${#jarr[@]}" -gt "0" ]; then
       write_template "${script}" "${jarr[@]}"
       echo "submitting ${jarr[@]}"
       cd production
       sbatch ${script}
       cd ..
    fi

""")

    fh.write("""
}



##############################################################################
""")
    pmemdMPI="%s.MPI"%("pmemd")
    ncpu=compinfo.num_cpu_per_node
    if is_gas:
        ncpu=2
        if is_piti:
            ncpu=1
            pmemdMPI="pmemd"
    fh.write("""

write_template() {
    local fname="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > production/${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}
EOF
    cat << 'EOF' >> production/${fname}
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
##SBATCH --exclude=cuda[001-008]
"""%(compinfo.cpu_partition,
     ncpu,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))

    
    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE="${AMBERHOME}/bin/%s"

prefix="%s"
"""%(compinfo.launch,
     pmemdMPI,
     basename))

    fh.write("lams=( %s )\n"%( " ".join(["%.8f"%(lam) for lam in lams]) ))
    fh.write("""EOF
    cat << 'EOF' >> production/${fname}
istep=${SLURM_ARRAY_TASK_ID}
""")
    fh.write("""

#
# is the input file asking for a replica-exchange simulation?
#
numexchg=$( grep numexchg ../inputs/${prefix}_${lams[0]}_analyze.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN ANALYSIS-MODE"
   exit 1
fi

#
# each step consists of running all windows for some length of time
#

istep=$(printf "%06i" ${istep})

#
# if don't have the subdirectory, then something's gone wrong
#

if [ ! -d ${istep} ]; then
   exit
fi

cd ${istep}

#
# run analysis for each window
#

for lambda in ${lams[@]}; do
    base="${prefix}_${lambda}"
    rest="${base}_restart"
    init="${base}_initial"
    anal="${base}_analyze"
""")
    if not is_piscti:
        fh.write("""
    ${LAUNCH} ${EXE} -O -c ${base}.rst7 -p ../../${prefix}.parm7 -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -y ${base}.nc
""")
    else:
        fh.write("""
    ${LAUNCH} ${EXE} -O -c ${base}.rst7 -p ../../${base}.parm7 -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -y ${base}.nc
""")
    fh.write("""
    for tmpfile in ${anal}.rst7 ${anal}.nc ${anal}.mdinfo logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

python2.7 ../stdti_step2dats.py *_analyze.mdout

ok=1
for lambda in ${lams[@]}; do
    if [ ! -e "dvdl_${lambda}.dat" ]; then
      ok=0
    fi
done
if [ "${ok}" == "1" ]; then
    for lambda in ${lams[@]}; do
       base="${prefix}_${lambda}"
       anal="${base}_analyze"
       if [ -e "${anal}.mdout" ]; then
          rm -f "${anal}.mdout"
       fi
    done
fi

EOF
    
}




write_python_script() {

cat << 'EOF' > production/stdti_step2dats.py
#!/usr/bin/env python2.7
import sys,os

def extract_traditional_ti( fname, write=False ):
    import os
    from collections import defaultdict as ddict

    fh = open(fname,"r")
    if not fh:
        raise Exception("Could not open %s\\n"%(fname))



    numexchg=0
    nstlim=None
    ntpr=None
    dt=None
    irest=0
    for line in fh:
        cmdstr,sepstr,comstr = line.partition("!")
        if "ntpr" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
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

    dt = dt
    nstep_per_sim = nstlim * numexchg
    nframe_per_sim = nstep_per_sim / ntpr

    if nstep_per_sim % ntpr != 0:
        print "num md steps per simulation is not a multiple of ntpr. Unclear how the simulation time works"

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
        lams = [ float(i) / ( nlam-1. ) for i in range(nlam) ]
        for l in lams:
            if abs(l-lam) < 0.001:
                lam = l
                break
        head, tail = os.path.split(fname)
        dvdl_fname = os.path.join( head, "dvdl_%.8f.dat"%( lam ) )

        if irest == 0:
           dvdls=dvdls[1:]

        fh = file(dvdl_fname,"w")
        for i in range(len(dvdls)):
            fh.write("%.4f %18.6f\\n"%((i+1)*t_per_frame,dvdls[i]))
        fh.close()
        for ilam,plam in enumerate(lams):
            efep_fname = os.path.join( head, "efep_%.8f_%.8f.dat"%( lam, plam ) )
            fh = file(efep_fname,"w")
            for i in range(len(efeps)):
                fh.write("%.4f %18.6f\\n"%((i+1)*t_per_frame,efeps[i][ilam]))
            fh.close()

    return dvdls,efeps


for arg in sys.argv[1:]:
    if os.path.isfile( arg ):
        if ".mdout" in arg:
            extract_traditional_ti( arg, write=True )
        else:
            print "File does not end in .mdout: %s"%(arg)
    else:
        print "File not found: %s"%(arg)

EOF

chmod u+x production/stdti_step2dats.py

}





# ---------------------------
# call to the main function
# ---------------------------

main

""")

    # make script executable by user
    import os
    import stat
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)
    



def write_piti_analysis_slurm_script(lams,basename,compinfo):
    import os.path
    fname = os.path.join( os.path.join(basename,"template"),"analyze.slurm.sh")
    fh = open(fname,"w")
    fh.write("""#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    else
       source ${AMBERHOME}/amber.sh
    fi

    has_parmed=$( python2.7 << 'EOF'
try:
    import parmed
    print "1"
except:
    print "0"
EOF
)

    if [ "${has_parmed}" != "1" ]; then
       echo "python2.7 could not import parmed."
       echo "Your PYTHONPATH variable may be incorrect."
       pathdir=$( grep 'export PYTHONPATH' ${AMBERHOME}/amber.sh | head -n 1 | sed -e 's/ *export PYTHONPATH="//' -e 'sA/site-packages"AA' )
       echo "${AMBERHOME}/amber.sh sets PYTHONPATH to ${pathdir}"
       if [ ! -d ${pathdir} ]; then
          echo "but that directory does not exist."
          exit 1
       else
          exit 1
       fi
    fi


    #######################################################
    local write_template=write_parallel_template
    #local write_template=write_serial_template
    #######################################################

    local base="%s" 
    local lams=( %s )"""%(basename," ".join("%.8f"%(lam) for lam in lams)))

    fh.write("""

    #
    # if we don't have any completed production, then there's nothing to analyze
    #

    if [ ! -d production ]; then 
       echo "Nothing to analyze because you don't have a production directory yet"
       exit
    fi

    #
    # the analysis will use a python utility to transform the data in the
    # re-analyzed mdout files to a series of .dat files
    #

    if [ ! -e "production/stdti_step2dats.py" ]; then
       write_python_script
    fi

    #
    # the production dir contains many subdirs that are 0-padded integers.
    # what is the largest number that we can find?
    #

    local last_step=$(for f in $(ls production | grep -v '\.' | tail -n 1); do bc -l <<< $f; done)

    #
    # if we couldn't find a valid subdirectory, then exit now
    #

    if [ "${last_step}" == "" ]; then 
       echo "Nothing to analyze because I couldn't find a valid subdirectory name in production/"
       exit
    fi



    #
    # collect the existing results, if any
    #

    echo "Collecting existing results before submitting new analysis..."
    python2.7 analysis2results.py

    echo ""
    echo ""
    echo "Searching for un-analyzed production directories..."


    local script="analysis.slurm"

    local jarr=()
    for step in $(seq ${last_step}); do

       #
       # this is the zero-padded name
       # 

       local step_name=$(printf "%06i" ${step})

       #
       # if the subdir does not exist, then skip it
       #

       if [ ! -d "production/${step_name}" ]; then
          echo "skipping ${step_name} because dir does not exist"
          continue
       fi


       #
       # do we actually have to analyze this subdir?
       #

       #
       # if it doesn't have mdout files, then I think I already deleted it to save disk space
       #

       local ok=1
       local cnt=0
       for lam in ${lams[@]}; do
           if [ -e "production/${step_name}/${base}_${lam}.mdout" ]; then
               cnt=$(( ${cnt} + 1 ))
           else
               ok=0
           fi
       done
       if [ "${cnt}" == "0" ]; then
          echo "skipping ${step_name} because it was probably already deleted from this machine"
          continue
       elif [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it has some, but not all, of the mdout files (something wrong here?)"
          continue
       fi


       ok=1
       for lam in ${lams[@]}; do
           if [ ! -e "production/${step_name}/${base}_${lam}.nc" ]; then
               ok=0
               echo "Missing production/${step_name}/${base}_${lam}.nc"
           fi
       done
       if [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it is missing 1-or-more nc files"
          continue
       fi


       #
       # do we have the mbar trace file generated from the analysis
       #

       local testfile=$(printf "production/${step_name}/efep_%.8f_%.8f.dat" 1. 1.)
       if [ -e ${testfile} ]; then
          echo "skipping ${step_name} because it has already been analyzed"
          continue
       fi


       local queued=$(squeue -ro "%i %o" | grep ${PWD} | grep ${script} | sed -e 's/ .*/ /' -e 's/.*_/_/' -e 's/ /_/' | grep "_${step}_")
       if [ -e "${queued}" ]; then
          echo "skipping ${step_name} because it is already queued"
          continue
       fi 
       jarr+=(${step})
    done

    if [ "${#jarr[@]}" -gt "0" ]; then
       ${write_template} "${script}" "${jarr[@]}"
       echo "submitting ${jarr[@]}"
       cd production
       sbatch ${script}
       cd ..
    fi
}



##############################################################################
""")
    
    fh.write("""

write_serial_template() {
    local fname="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > production/${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}
EOF
    cat << 'EOF' >> production/${fname}
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=0-06:00:00
##SBATCH --exclude=cuda[001-008]
"""%(compinfo.cpu_partition,
     compinfo.num_cpu_per_node))

    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        
    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE="${AMBERHOME}/bin/sander.MPI"

prefix="%s"
"""%(compinfo.launch,
     basename))

    fh.write("lams=( %s )\n"%( " ".join(["%.8f"%(lam) for lam in lams]) ))

    fh.write("""EOF
    cat << 'EOF' >> production/${fname}
istep=${SLURM_ARRAY_TASK_ID}
""")
    
    fh.write("""

#
# is the input file asking for a replica-exchange simulation?
#
numexchg=$( grep numexchg ../inputs/${prefix}_${lams[0]}_analyze.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN ANALYSIS-MODE"
   exit 1
fi

#
# each step consists of running all windows for some length of time
#

istep=$(printf "%06i" ${istep})

#
# if don't have the subdirectory, then something's gone wrong
#

if [ ! -d ${istep} ]; then
   echo "DIRECTORY NOT FOUND ${istep}"
   exit 1
fi

cd ${istep}

#############################################################
#
# THERMODYNAMIC INTEGRATION ANALYSIS
#

for lambda in ${lams[@]}; do
    base="${prefix}_${lambda}"
    rest="${base}_restart"
    init="${base}_initial"
    anal="${base}_analyze"

    ${LAUNCH} ${EXE} -O -c ${base}.rst7 -p ../../${base}.parm7 -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -y ${base}.nc

    for tmpfile in ${anal}.rst7 ${anal}.nc ${anal}.mdinfo ${anal}.mdout logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

#
# dU/dlam = dU/dParam * dParam/dlam
# The values are stored in a series of file: dvdl_%{lam}.dat
#

lam0=${lams[0]}
lam1=${lams[ $(( ${#lams[@]} - 1 )) ]}
python2.7 ../piti_chainrule.py ${prefix} ${#lams[@]} ../../${prefix}_${lam0}.parm7 ../../${prefix}_${lam1}.parm7 




#############################################################
#
# MBAR ANALYSIS
#

#
# The mbar analysis requires the calculation of the matrix
#   U( R_{tlam}; Param_{plam} )
# where R_{tlam} are the coordinates corresponding to trajectory tlam
# and Param_{plam} are the MM parameters corresponding to parameter file plam
# The values are stored in a series of files: efep_${tlam}_${plam}.dat
#

#
# imin=6 with ntpr=1 writes the .nc.dedparm file for TI
# however, for mbar, we don't need the .nc.dedparm file
# so we set it to 2
#

sed "s/ntpr *= *[0-9]*/ntpr = 2/" ../../inputs/${prefix}_${lam0}_analyze.mdin > mbar.mdin

for tlam in ${lams[@]}; do
    for plam in ${lams[@]}; do
       base="${prefix}_${tlam}_${plam}"
       traj="${prefix}_${tlam}"
       parm="${prefix}_${plam}"

       ${LAUNCH} ${EXE} -O -p ../../${parm}.parm7 -i mbar.mdin -o ${base}.mdout -c ${traj}.rst7 -y ${traj}.nc -r reanal.rst7 -inf ${base}.mdinfo

       grep 'minimization completed, ENE=' ${base}.mdout | sed -e 's/minimization completed, ENE= *//' -e 's/ RMS= .*//' | nl > efep_${tlam}_${plam}.dat

       for tmpfile in ${base}.mdinfo ${base}.mdout reanal.rst7 logfile; do
          if [ -e "${tmpfile}" ]; then
             rm -f "${tmpfile}"
          fi
       done
    done
done

if [ -e mbar.mdin ]; then
    rm -f mbar.mdin
fi


EOF

   ############################################################################## 
}

""")


    num_nodes,cpu_per_node = find_num_nodes_and_cores( len(lams), compinfo.num_cpu_per_node, 1, -1 )
    if is_gas:
        num_nodes,cpu_per_node = find_num_nodes_and_cores_gas( len(lams), compinfo.num_cpu_per_node, 1, -1 )
    fh.write("""

write_parallel_template() {
    local fname="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > production/${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}
EOF
    cat << 'EOF' >> production/${fname}
#SBATCH --partition=%s
#SBATCH --nodes=%s
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=0-06:00:00
##SBATCH --exclude=cuda[001-008]
"""%(compinfo.cpu_partition,
     num_nodes,
     cpu_per_node,
     compinfo.lau))

    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        
    
    fh.write("""

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE="${AMBERHOME}/bin/sander.MPI"

prefix="%s"
"""%(compinfo.launch,
     basename))

    fh.write("lams=( %s )\n"%( " ".join(["%.8f"%(lam) for lam in lams]) ))

    
    fh.write("""EOF
    cat << 'EOF' >> production/${fname}
istep=${SLURM_ARRAY_TASK_ID}
""")
    fh.write("""

#
# is the input file asking for a replica-exchange simulation?
#
numexchg=$( grep numexchg ../inputs/${prefix}_${lams[0]}_analyze.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN ANALYSIS-MODE"
   exit 1
fi

#
# each step consists of running all windows for some length of time
#

istep=$(printf "%06i" ${istep})

#
# if don't have the subdirectory, then something's gone wrong
#

if [ ! -d ${istep} ]; then
   echo "DIRECTORY NOT FOUND ${istep}"
   exit 1
fi

cd ${istep}

#############################################################
#
# THERMODYNAMIC INTEGRATION ANALYSIS
#

truncate -s0 ${prefix}_analyze.groupfile
for lambda in ${lams[@]}; do
    base="${prefix}_${lambda}"
    rest="${base}_restart"
    init="${base}_initial"
    anal="${base}_analyze"

    echo "-O -c ${base}.rst7 -p ../../${base}.parm7 -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -y ${base}.nc" >> ${prefix}_analyze.groupfile
done

${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile ${prefix}_analyze.groupfile

for lambda in ${lams[@]}; do
    base="${prefix}_${lambda}"
    anal="${base}_analyze"
    for tmpfile in ${anal}.rst7 ${anal}.nc ${anal}.mdinfo ${anal}.mdout logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f ${tmpfile}
       fi
    done
done

#
# dU/dlam = dU/dParam * dParam/dlam
# The values are stored in a series of file: dvdl_%{lam}.dat
#

lam0=${lams[0]}
lam1=${lams[ $(( ${#lams[@]} - 1 )) ]}
python2.7 ../piti_chainrule.py ${prefix} ${#lams[@]} ../../${prefix}_${lam0}.parm7 ../../${prefix}_${lam1}.parm7 




#############################################################
#
# MBAR ANALYSIS
#

#
# The mbar analysis requires the calculation of the matrix
#   U( R_{tlam}; Param_{plam} )
# where R_{tlam} are the coordinates corresponding to trajectory tlam
# and Param_{plam} are the MM parameters corresponding to parameter file plam
# The values are stored in a series of files: efep_${tlam}_${plam}.dat
#

#
# imin=6 with ntpr=1 writes the .nc.dedparm file for TI
# however, for mbar, we don't need the .nc.dedparm file
# so we set it to 2
#

sed "s/ntpr *= *[0-9]*/ntpr = 2/" ../../inputs/${prefix}_${lam0}_analyze.mdin > mbar.mdin

for tlam in ${lams[@]}; do

    truncate -s0 ${prefix}_analyze.groupfile

    for plam in ${lams[@]}; do
       base="${prefix}_${tlam}_${plam}"
       traj="${prefix}_${tlam}"
       parm="${prefix}_${plam}"

       echo "-O -p ../../${parm}.parm7 -i mbar.mdin -o ${base}.mdout -c ${traj}.rst7 -y ${traj}.nc -r ${base}.rst7 -inf ${base}.mdinfo" >> ${prefix}_analyze.groupfile
     done


     ${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile ${prefix}_analyze.groupfile


    for plam in ${lams[@]}; do
       base="${prefix}_${tlam}_${plam}"
       traj="${prefix}_${tlam}"
       parm="${prefix}_${plam}"

       grep 'minimization completed, ENE=' ${base}.mdout | sed -e 's/minimization completed, ENE= *//' -e 's/ RMS= .*//' | nl > efep_${tlam}_${plam}.dat

       for tmpfile in ${base}.mdinfo ${base}.mdout ${base}.rst7 logfile; do
          if [ -e "${tmpfile}" ]; then
            rm -f "${tmpfile}"
          fi
       done
    done

done

for tmpfile in mbar.mdin ${prefix}_analyze.groupfile; do
    if [ -e "${tmpfile}" ]; then
       rm -f "${tmpfile}"
    fi
done


EOF

   ############################################################################## 
}

""")



    fh.write("""


write_python_script() {

cat << 'EOF' > production/piti_chainrule.py
#!/usr/bin/env python2.7
import sys,os

        
def CptdEdLambda( param_0, param_1, lam, dedparm ):
    from math import sqrt

    from scipy.io import netcdf
    import numpy as np

    nc = netcdf.NetCDFFile(dedparm,'r')

    vdedq = nc.variables['q'] # kcal/mol per amber-charge
    nframe, npar = vdedq.shape
    if npar != len(param_0.atoms):
        raise Exception("CptdEdLambda num atom mistmatch %i %i"%(npar,len(param_0.atoms)))
    dedlam = [0.]*nframe
    for ipar in range(npar):
        dpardlam = 18.2223 * ( param_1.atoms[ipar].charge - param_0.atoms[ipar].charge )
        for iframe in range(nframe):
            dedq = vdedq.data[iframe,ipar]
            dedlam[iframe] += dedq * dpardlam
    
    ntypes = param_0.ptr('NTYPES')
    nttyp  = ntypes*(ntypes+1)/2
    pd = param_0.parm_data

    Rc0 = param_0.LJ_radius
    Rc1 = param_1.LJ_radius
    E0 = param_0.LJ_depth
    E1 = param_1.LJ_depth


    vdeda = nc.variables['a']
    vdedb = nc.variables['b']
    nframe, npar = vdeda.shape

    for i in range(ntypes):
        for j in range(0, i+1):
            index = param_0.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            if index < 0:
                continue

            Rm0 = Rc0[i] + Rc0[j]
            Rm1 = Rc1[i] + Rc1[j]
            e0  = sqrt( E0[i]*E0[j] )
            e1  = sqrt( E1[i]*E1[j] )
            Rm  = (1.-lam) * Rm0 + lam * Rm1
            e   = (1.-lam) * e0  + lam * e1
            
            A     = e * Rm**12
            B     = 2 * e * Rm**6
            
            dRdL = Rm1-Rm0
            dedL = e1-e0
            
            dAde  = Rm**12
            dAdR  = e * 12 * Rm**11
            dBde  = 2 * Rm**6
            dBdR  = 2 * e * 6 * Rm**5
            dadl = dAde*dedL + dAdR*dRdL
            dbdl = dBde*dedL + dBdR*dRdL
            for iframe in range(nframe):
                deda = vdeda.data[iframe,index]
                dedb = vdedb.data[iframe,index]
                dedlam[iframe] += deda*dadl + dedb*dbdl

    vdedk = nc.variables['rk']
    vdedq = nc.variables['req']
    for var,par in [ ("rk","BOND_FORCE_CONSTANT"), 
                     ("req","BOND_EQUIL_VALUE"), 
                     ("tk","ANGLE_FORCE_CONSTANT"), 
                     ("teq","ANGLE_EQUIL_VALUE"), 
                     ("pk","DIHEDRAL_FORCE_CONSTANT") ]:
        vdedp = nc.variables[var]
        nframe, npar = vdedp.shape
        for ipar in range(npar):
            dpdl = param_1.parm_data[par][ipar] - param_0.parm_data[par][ipar]
            for iframe in range(nframe):
                dedp = vdedp.data[iframe,ipar]
                dedlam[iframe] += dedp*dpdl
    return dedlam

        


def analyze_interpolated_charge_dedparms( parmfile_lambda0, parmfile_lambda1, prefix, nquad ):
    import parmed
    from collections import defaultdict as ddict
    param_0 = parmed.load_file(parmfile_lambda0)
    param_1 = parmed.load_file(parmfile_lambda1)
    pts = [ float(i) / ( nquad - 1. ) for i in range(nquad) ]
    lam0 = pts[0]
    dEdLambdas = ddict( float )
    for lam in pts:
        dedparm = "%s_%.8f.nc.dedparm"%(prefix,lam)
        dEdLambdas[lam] = CptdEdLambda( param_0, param_1, lam, dedparm )
    nframes = len(dEdLambdas[lam0])
    for lam in dEdLambdas:
        fh = file("dvdl_%.8f.dat"%(lam),"w")
        for frame in range(nframes):
            fh.write("%8i %15.6f\\n"%(frame+1,dEdLambdas[lam][frame]))
        fh.close()
        


prefix = sys.argv[1]
nquad = int(sys.argv[2])
parm0 = sys.argv[3]
parm1 = sys.argv[4]

analyze_interpolated_charge_dedparms(parm0,parm1,prefix,nquad)

EOF

chmod u+x production/piti_chainrule.py

}





# ---------------------------
# call to the main function
# ---------------------------

main

""")

    # make script executable by user
    import os
    import stat
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)
    


    

def write_stdti_heating_slurm_script(lams,basename,compinfo,is_piti,is_piscti,is_gas=False):
    import os.path
    fname = os.path.join( os.path.join(basename,"template"),"equilibration.slurm.sh")
    lamsstr="lams=( %s )\n"%( " ".join(["%.8f"%(lam) for lam in lams]) )

    pmemd="pmemd"
    if is_piscti:
        pmemd="pmemdPI"


    sc="#"
    pc="#"
    sg="#"
    if compinfo.serial_cpu:
        sc=""
    if compinfo.parallel_cpu:
        pc=""
    if compinfo.serial_gpu:
        sg=""
        #pmemd=pmemd+".cuda"
    if compinfo.parallel_gpu:
        sg=""
        #pmemd=pmemd+".cuda"
        
    fh = open(fname,"w")
    fh.write("""#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi

    ###################################################################
    %slocal write_template=write_parallel_cpu_template
    %slocal write_template=write_serial_cpu_template
    %slocal write_template=write_serial_gpu_template
    ###################################################################
"""%(pc,sc,sg))

    fh.write("    local template=equilibration.slurm")

    fh.write("""
    ${write_template} "${template}"

    #
    # do we already have stuff running in this directory?
    # if so, then append the current jobs to the end of the dependency chain
    #

    local lastid=$(squeue -h -o "%i %o" | grep ${PWD} | sed 's/ .*//' | sort -rn | head -n 1)
    if [ "${lastid}" != "" ]; then
       sbatch --dependency=afterany:${lastid} ${template}
    else
       sbatch ${template}
    fi

    
}

""")

    num_nodes,num_cores_per_node = find_num_nodes_and_cores( len(lams), compinfo.num_cpu_per_node, compinfo.min_cpu_nodes, compinfo.max_cpu_nodes )
    if is_gas:
        num_nodes,cpu_per_node = find_num_nodes_and_cores_gas( len(lams), compinfo.num_cpu_per_node, compinfo.min_cpu_nodes, compinfo.max_cpu_nodes )


    fh.write("""
##############################################################################

write_parallel_cpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=%s
#SBATCH --nodes=%i
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
"""%(compinfo.cpu_partition,
     num_nodes,
     num_cores_per_node,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
        

    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s.MPI
prefix="%s"
%s
EOF
write_parallel_template "${fname}"
}
"""%(compinfo.launch,
     pmemd,
     basename,
     lamsstr))
    

    fh.write("""
write_parallel_template() {
    local fname="$1"
cat << 'EOF' >> ${fname}
if [ ! -d equilib ]; then
    mkdir equilib
fi
if [ ! -d current ]; then
    mkdir current
fi

${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_minimiz.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_heating.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press01.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press02.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press03.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press04.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press05.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press06.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press07.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press08.groupfile

for lam in ${lams[@]}; do
    base=${prefix}_${lam}
    rest=${base}_restart
    heat=${base}_heating
    pre1=${base}_press01
    pre2=${base}_press02
    pre3=${base}_press03
    pre4=${base}_press04
    pre5=${base}_press05
    pre6=${base}_press06
    pre7=${base}_press07
    pre8=${base}_press08

    cp equilib/${pre8}.rst7 current/${rest}.rst7

    for tmpfile in current/${base}.mdout current/${base}.nc current/${base}.mdinfo current/${base}.logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

EOF
}
""")


    parmstr="${prefix}.parm7"
    if is_piti:
        parmstr="${base}.parm7"

    pmemdMPI="%s.MPI"%(pmemd)
    ncpu=compinfo.num_cpu_per_node
    if is_gas:
        ncpu=2
        if is_piti:
            ncpu=1
            pmemdMPI=pmemd
        
    fh.write("""
##############################################################################
write_serial_cpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --array=1-%i
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=%i
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
##SBATCH --exclude=cuda[001-008]
"""%(len(lams),
     compinfo.cpu_partition,
     ncpu,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
    
    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s.MPI
prefix="%s"
%s
export parm=%s

EOF
write_serial_template ${fname}
}
"""%(compinfo.launch,
     pmemd,
     basename,
     lamsstr,parmstr))

    
    fh.write("""
##############################################################################
write_serial_gpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --array=1-%i
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mincpus=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#
##SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --export=ALL
#SBATCH --time=%i-00:00:00
"""%(len(lams),
     compinfo.gpu_partition,
     compinfo.num_days))
    
    if len(compinfo.qos) > 0:
        fh.write("#SBATCH --qos=%s\n"%(compinfo.qos))
    if len(compinfo.account) > 0:
        fh.write("#SBATCH --account=%s\n"%(compinfo.account))
    if len(compinfo.exclude) > 0:
        fh.write("#SBATCH --exclude=%s\n"%(compinfo.exclude))
    if len(compinfo.constraint) > 0:
        fh.write("#SBATCH --constraint=%s\n"%(compinfo.constraint))
    
    fh.write("""
export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="%s"


export EXE=${AMBERHOME}/bin/%s.cuda
prefix="%s"
%s
export parm=%s

EOF
write_serial_template ${fname}
}
"""%(compinfo.launch,
     pmemd,
     basename,
     lamsstr,parmstr))

    
    fh.write("""
write_serial_template() {
    local fname="$1"
cat <<'EOF' >> ${fname} 
if [ ! -d equilib ]; then
    mkdir equilib
fi
if [ ! -d current ]; then
    mkdir current
fi


myidx=$(( ${SLURM_ARRAY_TASK_ID} - 1 ))
mylam=${lams[${myidx}]}

for lam in ${mylam}; do
    base=${prefix}_${lam}
    init=${base}_initial
    rest=${base}_restart
    mini=${base}_minimiz
    heat=${base}_heating
    pre1=${base}_press01
    pre2=${base}_press02
    pre3=${base}_press03
    pre4=${base}_press04
    pre5=${base}_press05
    pre6=${base}_press06
    pre7=${base}_press07
    pre8=${base}_press08

    ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o equilib/${mini}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${mini}.rst7 -p ${parm} -i inputs/${heat}.mdin -o equilib/${heat}.mdout -r equilib/${heat}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${heat}.rst7 -p ${parm} -i inputs/${pre1}.mdin -o equilib/${pre1}.mdout -r equilib/${pre1}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre1}.rst7 -p ${parm} -i inputs/${pre2}.mdin -o equilib/${pre2}.mdout -r equilib/${pre2}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre1}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre2}.rst7 -p ${parm} -i inputs/${pre3}.mdin -o equilib/${pre3}.mdout -r equilib/${pre3}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre2}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre3}.rst7 -p ${parm} -i inputs/${pre4}.mdin -o equilib/${pre4}.mdout -r equilib/${pre4}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre3}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre4}.rst7 -p ${parm} -i inputs/${pre5}.mdin -o equilib/${pre5}.mdout -r equilib/${pre5}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre5}.rst7 -p ${parm} -i inputs/${pre6}.mdin -o equilib/${pre6}.mdout -r equilib/${pre6}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre6}.rst7 -p ${parm} -i inputs/${pre7}.mdin -o equilib/${pre7}.mdout -r equilib/${pre7}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre7}.rst7 -p ${parm} -i inputs/${pre8}.mdin -o equilib/${pre8}.mdout -r equilib/${pre8}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7


    cp equilib/${pre8}.rst7 current/${rest}.rst7

    for tmpfile in current/${base}.mdout current/${base}.nc current/${base}.mdinfo logfile*; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

EOF
}
""")

    fh.write("""

# ---------------------------
# call to the main function
# ---------------------------

main

""")

    # make script executable by user
    import os
    import stat
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)
    

##########################################################
##########################################################
# NON-SOFTCORE INTERPOLATION AND ANALYSIS
##########################################################
##########################################################

def nonsoftcore_interpolation(param_0,param_1,lam):
    from math import sqrt

    Rc0 = param_0.LJ_radius
    Rc1 = param_1.LJ_radius
    E0  = param_0.LJ_depth
    E1  = param_1.LJ_depth

    if True:
        param = CopyParm( param_0 )
        
        ##########################################################
        for o,a,b in zip( param, param_0.atoms, param_1.atoms ):
            o.charge  = (1.-lam) * a.charge  + lam * b.charge
        ##########################################################
        pd = param.parm_data
        ntypes = param.pointers['NTYPES']
        for i in range(ntypes):
            for j in range(i, ntypes):
                index = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                if index < 0:
                    continue
                Rm0 = Rc0[i] + Rc0[j]
                Rm1 = Rc1[i] + Rc1[j]
                e0  = sqrt( E0[i]*E0[j] )
                e1  = sqrt( E1[i]*E1[j] )
                Rm  = (1.-lam) * Rm0 + lam * Rm1
                e   = (1.-lam) * e0  + lam * e1
                pd["LENNARD_JONES_ACOEF"][index] = e * Rm**12
                pd["LENNARD_JONES_BCOEF"][index] = 2 * e * Rm**6

        npar0 = len(param_0.bond_types)
        npar1 = len(param_1.bond_types)
        if npar0 != npar1:
            raise Exception("bond_types size mismatch %i %i"%(npar0,npar1))
        for i in range(npar0):
            t0=param_0.bond_types[i]
            t1=param_1.bond_types[i]
            param.bond_types[i].k   = (1.-lam)*t0.k   + lam*t1.k
            param.bond_types[i].req = (1.-lam)*t0.req + lam*t1.req
        npar0 = len(param_0.angle_types)
        npar1 = len(param_1.angle_types)
        if npar0 != npar1:
            raise Exception("angle_types size mismatch %i %i"%(npar0,npar1))
        for i in range(npar0):
            t0=param_0.angle_types[i]
            t1=param_1.angle_types[i]
            param.angle_types[i].k      = (1.-lam)*t0.k      + lam*t1.k
            param.angle_types[i].theteq = (1.-lam)*t0.theteq + lam*t1.theteq
        npar0 = len(param_0.dihedral_types)
        npar1 = len(param_1.dihedral_types)
        if npar0 != npar1:
            raise Exception("dihedral_types size mismatch %i %i"%(npar0,npar1))
        for i in range(npar0):
            t0=param_0.dihedral_types[i]
            t1=param_1.dihedral_types[i]
            param.dihedral_types[i].phi_k = (1.-lam)*t0.phi_k + lam*t1.phi_k
            if abs(param.dihedral_types[i].phi_k) < 1.e-10:
                param.dihedral_types[i].phi_k = 0
            if t0.per != t1.per:
                raise Exception("dihedral %i periodicity mismatch %i %i"%(i+1,t0.per,t1.per))
            if t0.phase != t1.phase:
                raise Exception("dihedral %i phase mismatch %i %i"%(i+1,t0.phase,t1.phase))
    return param


def nonsoftcore_tigradient(param_0,param_1,lam,dedparm_nc_file):
    from math import sqrt
    from scipy.io import netcdf
    import numpy as np

    nc = netcdf.NetCDFFile(dedparm_nc_file,'r')

    vdedq = nc.variables['q'] # kcal/mol per amber-charge
    nframe, npar = vdedq.shape
    if npar != len(param_0.atoms):
        raise Exception("nonsoftcore_tigradient num atom mismatch %i %i"%(npar,len(param_0.atoms)))
    dedlam = [0.]*nframe
    for ipar in range(npar):
        dpardlam = 18.2223 * ( param_1.atoms[ipar].charge - param_0.atoms[ipar].charge )
        for iframe in range(nframe):
            dedq = vdedq.data[iframe,ipar]
            dedlam[iframe] += dedq * dpardlam
#            if abs(dpardlam) > 0.0001:
#                print "%5i %8.4f %8.4f %8.4f q0 %12.6f q1 %12.6f"\
#                    %(ipar+1,dedlam[iframe],dedq,dpardlam,
#                      param_0.atoms[ipar].charge,
#                      param_1.atoms[ipar].charge)
    
    ntypes = param_0.ptr('NTYPES')
    nttyp  = ntypes*(ntypes+1)/2
    pd = param_0.parm_data

    Rc0 = param_0.LJ_radius
    Rc1 = param_1.LJ_radius
    E0 = param_0.LJ_depth
    E1 = param_1.LJ_depth

    vdeda = nc.variables['a']
    vdedb = nc.variables['b']
    nframe, npar = vdeda.shape

    for i in range(ntypes):
        for j in range(0, i+1):
            index = param_0.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            if index < 0:
                continue

            Rm0 = Rc0[i] + Rc0[j]
            Rm1 = Rc1[i] + Rc1[j]
            e0  = sqrt( E0[i]*E0[j] )
            e1  = sqrt( E1[i]*E1[j] )
            Rm  = (1.-lam) * Rm0 + lam * Rm1
            e   = (1.-lam) * e0  + lam * e1
            
            A     = e * Rm**12
            B     = 2 * e * Rm**6
            
            dRdL = Rm1-Rm0
            dedL = e1-e0
            
            dAde  = Rm**12
            dAdR  = e * 12 * Rm**11
            dBde  = 2 * Rm**6
            dBdR  = 2 * e * 6 * Rm**5
            dadl = dAde*dedL + dAdR*dRdL
            dbdl = dBde*dedL + dBdR*dRdL
            for iframe in range(nframe):
                deda = vdeda.data[iframe,index]
                dedb = vdedb.data[iframe,index]
                dedlam[iframe] += deda*dadl + dedb*dbdl

    vdedk = nc.variables['rk']
    vdedq = nc.variables['req']
    for var,par in [ ("rk","BOND_FORCE_CONSTANT"), 
                     ("req","BOND_EQUIL_VALUE"), 
                     ("tk","ANGLE_FORCE_CONSTANT"), 
                     ("teq","ANGLE_EQUIL_VALUE"), 
                     ("pk","DIHEDRAL_FORCE_CONSTANT") ]:
        vdedp = nc.variables[var]
        nframe, npar = vdedp.shape
        for ipar in range(npar):
            dpdl = param_1.parm_data[par][ipar] - param_0.parm_data[par][ipar]
            for iframe in range(nframe):
                dedp = vdedp.data[iframe,ipar]
                dedlam[iframe] += dedp*dpdl
    #print "final %8.4f"%(dedlam[0])
    return dedlam
    



##########################################################
##########################################################
# SOFTCORE INTERPOLATION AND ANALYSIS
##########################################################
##########################################################

# UTILITY FUNCTIONS

def Get_BondedCommonAtoms( parm, s1 ):
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    common = []
    for x in parm.bonds:
        has1 = x.atom1.idx in s1
        has2 = x.atom2.idx in s1
        if has1 and not has2:
            common.append( x.atom2.idx )
        elif has2 and not has1:
            common.append( x.atom1.idx )
    
    for x in parm.angles:
        has1 = x.atom1.idx in s1
        has2 = x.atom3.idx in s1
        if has1 and not has2:
            common.append( x.atom3.idx )
        elif has2 and not has1:
            common.append( x.atom1.idx )

    for x in parm.dihedrals:
        has1 = x.atom1.idx in s1
        has2 = x.atom4.idx in s1
        if has1 and not has2:
            common.append( x.atom4.idx )
        elif has2 and not has1:
            common.append( x.atom1.idx )
            
    return list(set(common))


def Get_TorsionedCommonAtoms( parm, s1 ):
    from collections import defaultdict as ddict
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    common = ddict( list )
    for x in parm.dihedrals:
        has1 = x.atom1.idx in s1
        has2 = x.atom4.idx in s1
        if has1 and not has2:
            common[x.atom1.idx].append( x.atom4.idx )
        elif has2 and not has1:
            common[x.atom4.idx].append( x.atom1.idx )
    for x in common:
        common[x] = list(set(common[x]))
    return common

def Get_TorsionedCommonAtomTypes( parm, s1 ):
    from collections import defaultdict as ddict
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    common = ddict( list )
    for x in parm.dihedrals:
        has1 = x.atom1.idx in s1
        has2 = x.atom4.idx in s1
        if has1 and not has2:
            common[x.atom1.nb_idx-1].append( x.atom4.nb_idx-1 )
        elif has2 and not has1:
            common[x.atom4.nb_idx-1].append( x.atom1.nb_idx-1 )
    for x in common:
        common[x] = list(set(common[x]))
    return common


def SetSCTI_BondedLJTypes( parm, common, dbgprint=False ):
    for i in common:
        sele = ListToSelection( [i] )
        if dbgprint:
            print "New LJ type for common atom (because it torsions):",sele
        parmed.tools.addLJType( parm, sele ).execute()

        
def SetSCTI_LJTypes( parm, s1, dbgprint=False ):
    from collections import defaultdict as ddict
    d1 = ddict( list )
    common = Get_TorsionedCommonAtoms( parm, s1 )
    t1 = []
    for iat in s1:
        skip=False
        if iat in common:
            if len(common[iat]) > 0:
                skip=True
                t1.append( iat ) # iat is a common atom that torsions with s1
        if not skip:
            d1[parm.atoms[iat].nb_idx].append( iat )
    for typ in d1:
        sele = ListToSelection( d1[typ] )
        if dbgprint:
            print "New LJ type for tireg atoms:",sele
        parmed.tools.addLJType( parm, sele ).execute()
    for i in t1:
        sele = ListToSelection( [i] )
        if dbgprint:
            print "New LJ type for tireg atoms (because it torsions):",sele
        parmed.tools.addLJType( parm, sele ).execute()

        
def PrintExclusions( parm, s1 ):
    # printing
    for iat in s1:
        a = parm.atoms[iat]
        aname = "%3s%04i:%-4s [%4i]"%(a.residue.name,a.residue.idx+1,a.name,a.idx+1)
        es = a.exclusion_partners
        ex = ListToSelection( [ b.idx for b in a.exclusion_partners ] )
        print "%s EXCLUDES: %s"%(aname,ex)

    
def SetSCTI_Exclusions( parm, s1, s2, dbgprint=False ):
    for iat in s1:
        a = parm.atoms[iat]
        for jat in s2:
            b = parm.atoms[jat]
            a.exclude(b)
    if dbgprint:
        PrintExclusions( parm, s1+s2 )
    
    

def SetSCTI_Regions( parm, s1, s2, dbgprint=False ):
    SetSCTI_LJTypes( parm, s1, dbgprint ) # define new LJ types for TI region 1
    SetSCTI_LJTypes( parm, s2, dbgprint ) # define new LJ types for TI region 2
    SetSCTI_Exclusions( parm, s1, s2, dbgprint ) # exclude pairs between the TI region
    common = Get_BondedCommonAtoms( parm, s1+s2 )
    SetSCTI_BondedLJTypes( parm, common )
    # You will need &ewald chngmask=0 &end in the mdin file
    nttypes = len(parm.parm_data['LENNARD_JONES_ACOEF'])
    # You will need to set &cntrl vdwmodel=2 &end in the mdin file
    parm.add_flag('LENNARD_JONES_LCOEF', '5E16.8', data=[0.]*nttypes)


    

def GetAtomsOfType( parm, ityp ):
    atms=[]
    for a in parm.atoms:
        if a.nb_idx-1 == ityp:
            atms.append( a.idx )
    return atms

def SetSCTI_LJ( parm, s1, s12, lam=0 ):
    #parm.recalculate_LJ()
    s1types  = [ parm.atoms[iat].nb_idx - 1 for iat in s1 ]
    s12types = [ parm.atoms[iat].nb_idx - 1 for iat in s12 ]
    torsionmap = Get_TorsionedCommonAtomTypes( parm, s12 )

    
    pd = parm.parm_data
    ntypes = parm.pointers['NTYPES']
    for i in range(ntypes):
        for j in range(i, ntypes):
            index = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            if index < 0:
                continue

            
            if (i in s1types and j not in s12types) or (j in s1types and i not in s12types):
                skip=False
                if i in s1types:
                    if i in torsionmap:
                        if j in torsionmap[i]:
                            skip=True
                if j in s1types:
                    if j in torsionmap:
                        if i in torsionmap[j]:
                            skip=True
                if skip:
                    iats = GetAtomsOfType(parm,i)
                    jats = GetAtomsOfType(parm,j)
                    #print "Skipping ",i,j,"due to torsions:",iats,jats
                if not skip:
                    A = pd["LENNARD_JONES_ACOEF"][index]
                    B = pd["LENNARD_JONES_BCOEF"][index]
                    if abs(B) > 1.e-6:
                        pd["LENNARD_JONES_LCOEF"][index] = (1.-lam) * ( A / (2.*B) )
                    else:
                        pd["LENNARD_JONES_LCOEF"][index] = 0.
                    pd["LENNARD_JONES_ACOEF"][index] *= lam
                    pd["LENNARD_JONES_BCOEF"][index] *= lam




def softcore_interpolation(inpparm,scmask1,scmask2,lam):
    parm = CopyParm( inpparm )
    if scmask1 == '@0':
        scmask1 = ''
    if scmask2 == '@0':
        scmask2 = ''
    s1 = GetSelectedAtomIndices( parm, scmask1 )
    s2 = GetSelectedAtomIndices( parm, scmask2 )
    s12 = s1+s2
    SetSCTI_Regions( parm, s1, s2, dbgprint=False )
    SetSCTI_LJ( parm, s1, s12, lam=(1-lam) )
    SetSCTI_LJ( parm, s2, s12, lam=lam )
    return parm


def softcore_tigradient( param_0, param_1, lam, dedparm_nc_file, scmask1, scmask2 ):
    from math import sqrt
    from scipy.io import netcdf
    import numpy as np
    from collections import defaultdict as ddict

    nc = netcdf.NetCDFFile(dedparm_nc_file,'r')

    vdedq = nc.variables['q'] # kcal/mol per amber-charge
    nframe, npar = vdedq.shape
    if npar != len(param_0.atoms):
        raise Exception("softcore_tigradient num atom mistmatch %i %i"%(npar,len(param_0.atoms)))
    
    dedlam = [0.]*nframe
    for ipar in range(npar):
        dpardlam = 18.2223 * ( param_1.atoms[ipar].charge - param_0.atoms[ipar].charge )
        for iframe in range(nframe):
            dedq = vdedq.data[iframe,ipar]
            dedlam[iframe] += dedq * dpardlam

    #print "dedq",dedlam[0]
    
    ntypes = param_0.ptr('NTYPES')
    nttyp  = ntypes*(ntypes+1)/2
    pd = param_0.parm_data

    Rc0 = param_0.LJ_radius
    Rc1 = param_1.LJ_radius
    E0 = param_0.LJ_depth
    E1 = param_1.LJ_depth

    vdeda = nc.variables['a']
    vdedb = nc.variables['b']
    vdedl = nc.variables['l']
    ljlrc = nc.variables['lrc']
    nframe, npar = vdeda.shape
#    print "nframe, npar = %i %i"%(nframe,npar)
#    print len(param_0.parm_data["LENNARD_JONES_ACOEF"])

    if False:
        disappearing_types = []
        appearing_types = []
        for i in range(ntypes):
            for j in range(0, i+1):
                index = param_0.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                if index < 0:
                    continue
                l0 = param_0.parm_data["LENNARD_JONES_LCOEF"][index]
                l1 = param_1.parm_data["LENNARD_JONES_LCOEF"][index]
                if abs(l1-l0) < 1.e-6:
                    continue
                if l1-l0 > 0:
                    #print "dis %3i %3i %12.3f %12.3f"%(i,j,l0,l1)
                    disappearing_types.append(i)
                    disappearing_types.append(j)
                else:
                    #print "app %3i %3i %12.3f %12.3f"%(i,j,l0,l1)
                    appearing_types.append(i)
                    appearing_types.append(j) 
        disappearing_types = set(disappearing_types)
        appearing_types = set(appearing_types)
        common_types = set()
        for i in range(ntypes):
            if i in disappearing_types and i in appearing_types:
                disappearing_types.remove(i)
                appearing_types.remove(i)
                common_types.add(i)

            
        disappearing_atoms = []
        appearing_atoms = []
        Ntype_disappearing = ddict( int )
        Ntype_appearing = ddict( int )
        Ntype_simulated = ddict( int )
        for a in param_0.atoms:
            i = a.nb_idx-1
            Ntype_simulated[ i ] += 1
            if i in disappearing_types:
                Ntype_disappearing[ i ] += 1
                disappearing_atoms.append( a.idx )
            elif i in appearing_types:
                Ntype_appearing[ i ] += 1
                appearing_atoms.append( a.idx )
            else:
                Ntype_disappearing[ i ] += 1
                Ntype_appearing[ i ] += 1
        #    print "disappearing atoms: %s"%(ListToSelection( disappearing_atoms ) )
        #    print "   appearing atoms: %s"%(ListToSelection( appearing_atoms ) )
        #    for a in disappearing_atoms:
        #        i = param_0.atoms[a].nb_idx - 1
        #        print "disappearing atom %i => %i"%(a,i)
    else:
        disappearing_atoms = GetSelectedAtomIndices( param_0, scmask1 )
        appearing_atoms = GetSelectedAtomIndices( param_0, scmask2 )
        Ntype_disappearing = ddict( int )
        Ntype_appearing = ddict( int )
        Ntype_simulated = ddict( int )
        for a in param_0.atoms:
            i = a.nb_idx-1
            Ntype_simulated[ i ] += 1
            if a.idx in disappearing_atoms:
                Ntype_disappearing[ i ] += 1
            elif a.idx in appearing_atoms:
                Ntype_appearing[ i ] += 1
            else:
                Ntype_disappearing[ i ] += 1
                Ntype_appearing[ i ] += 1
            
    # ==============================================================================
    # Sander computes a d/dB contribution from the LJ long-range-correction,
    # but the LRC incorrectly includes BOTH TI regions, so we need to recompute
    # and remove that particular contribution here.
    term_0 = 0.
    term_1 = 0.
    vdedb_correction = [0] * nframe * len(param_0.parm_data['NONBONDED_PARM_INDEX'])
    for i in range(ntypes):
        for j in range(0, i+1):
            index = param_0.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            if index < 0:
                continue
            B0 = param_0.parm_data["LENNARD_JONES_BCOEF"][index]
            B1 = param_1.parm_data["LENNARD_JONES_BCOEF"][index]
            B = B0 + lam*(B1-B0)
            NNC = -Ntype_simulated[i]*Ntype_simulated[j]
            NNC_0 = -Ntype_disappearing[i]*Ntype_disappearing[j]*B0
            NNC_1 = -Ntype_appearing[i]*Ntype_appearing[j]*B1
            if i != j:
                NNC *= 2.
                NNC_0 *= 2.
                NNC_1 *= 2.
            term_0 += NNC_0
            term_1 += NNC_1
            for iframe in range(nframe):
                vdedb_correction[iframe + index*nframe] = - ljlrc.data[iframe] * NNC

    # ==============================================================================
    # linear-TI transformation of the LJ long-range-correction term
    for iframe in range(nframe):
        dedlam[iframe] += ljlrc[iframe] * ( term_1-term_0 )

#    print "dedlrc",dedlam[0]

    for i in range(ntypes):
        for j in range(0, i+1):
            index = param_0.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            if index < 0:
                continue

            l0 = param_0.parm_data["LENNARD_JONES_LCOEF"][index]
            l1 = param_1.parm_data["LENNARD_JONES_LCOEF"][index]
            dldl = l1-l0
            a0 = param_0.parm_data["LENNARD_JONES_ACOEF"][index]
            a1 = param_1.parm_data["LENNARD_JONES_ACOEF"][index]
            dadl = a1-a0
            b0 = param_0.parm_data["LENNARD_JONES_BCOEF"][index]
            b1 = param_1.parm_data["LENNARD_JONES_BCOEF"][index]
            dbdl = b1-b0
            
            for iframe in range(nframe):
                #print len(param_0.parm_data["LENNARD_JONES_ACOEF"])
                deda = vdeda.data[iframe,index]
                dedb = vdedb.data[iframe,index] + vdedb_correction[iframe+index*nframe]
                dedl = vdedl.data[iframe,index]
                dedlam[iframe] += deda*dadl + dedb*dbdl + dedl*dldl

#    print "dedlj",dedlam[0]

    vdedk = nc.variables['rk']
    vdedq = nc.variables['req']
    for var,par in [ ("rk","BOND_FORCE_CONSTANT"), 
                     ("req","BOND_EQUIL_VALUE"), 
                     ("tk","ANGLE_FORCE_CONSTANT"), 
                     ("teq","ANGLE_EQUIL_VALUE"), 
                     ("pk","DIHEDRAL_FORCE_CONSTANT") ]:
        vdedp = nc.variables[var]
        nframe, npar = vdedp.shape
        for ipar in range(npar):
            dpdl = param_1.parm_data[par][ipar] - param_0.parm_data[par][ipar]
            for iframe in range(nframe):
                dedp = vdedp.data[iframe,ipar]
                dedlam[iframe] += dedp*dpdl
    return dedlam




##########################################################
##########################################################
# TESTING
##########################################################
##########################################################
    

def get_stdti_dvdl(parm,parmfilename,basename,rst,nc,lam,timask1,timask2,scmask1,scmask2):
    import subprocess
    mdin = Mdin( base = basename, inpcrd = rst )
    mdin.RST7 = "restart"
    mdin.NC = "mdcrd"
    mdin.INPNC = basename+".nc"
    mdin.PARM7 = parmfilename
    mdin.exe = "${AMBERHOME}/bin/pmemd"
    mdin.Set_NVT()
    mdin.Set_Restart(False)
    mdin.Set_PrintFreq(0)
    mdin.cntrl["nstlim"] = 0
    mdin.cntrl["ntpr"] = 1
    mdin.cntrl["icfe"] = 1
    mdin.cntrl["ifsc"] = 0
    mdin.cntrl["ntf"]  = 1
    mdin.cntrl["cut"]  = 10.
    
    mdin.cntrl["imin"]  = 6
    mdin.cntrl["ntmin"] = 1

    scsel1 = GetSelectedAtomIndices( parm, scmask1 )
    scsel2 = GetSelectedAtomIndices( parm, scmask2 )
    if len( scsel1+scsel2 ) > 0:
        mdin.cntrl["ifsc"] = 1
    mdin.cntrl["crgmask"] = "'%s'"%( ListToSelection(scsel1+scsel2) )
    
    mdin.cntrl["timask1"] = "'%s'"%(timask1)
    mdin.cntrl["timask2"] = "'%s'"%(timask2)
    mdin.cntrl["scmask1"] = "'%s'"%(scmask1)
    mdin.cntrl["scmask2"] = "'%s'"%(scmask2)
    mdin.cntrl["clambda"] = "%.8f"%(lam)

    mdin.WriteMdin()
    cmd = "%s %s"%(mdin.exe,mdin.CmdString())
    #print cmd
    subprocess.call(cmd,shell=True)
    
    fh = file(mdin.MDOUT,"r")
    dvdl=None
    for line in fh:
        if "DV/DL" in line:
            cs = line.strip().split()
            for i in range(len(cs)):
                if "DV/DL" in cs[i]:
                    dvdl = float(cs[i+2])
                    break
            break
    return dvdl


def calc_stdti_dvdls(tiobj):
    import os.path
    
    nlam = 21
    lams = [ float(i)/(nlam-1.) for i in range(nlam) ]



    deq=0.
    has_deq = False
    try:
        if tiobj.deqti is not None:
            has_deq = True
    except:
        pass
    
    sc=0.
    has_sc  = False
    try:
        if tiobj.scti is not None:
            has_sc = True
    except:
        pass
    
    req=0.
    has_req = False
    try:
        if tiobj.reqti is not None:
            has_req = True
    except:
        pass

    
    
    if has_deq:
        o        = tiobj.deqti
        parm     = o.tiparm
        dirname  = os.path.join(o.basename,"template")
        basename = "test.%s"%(o.basename)
        parmfilename = os.path.join(dirname,"%s.parm7"%(o.basename))
        timask1  = o.merged_timask1
        timask2  = o.merged_timask2
        scmask1  = o.merged_scmask1
        scmask2  = o.merged_scmask2
        
        rst =  basename+".rst7"
        nc  =  basename+".nc"
        parmed.tools.writeCoordinates(parm, rst).execute()
        parmed.tools.writeCoordinates(parm, nc).execute()

    

        if scmask1 == "" and scmask2 == "":
            dvdl0 = get_stdti_dvdl(parm,parmfilename,basename,rst,nc,0,timask1,timask2,scmask1,scmask2)
            dvdl = [dvdl0]*(nlam)
        else:
            dvdl = []
            for lam in lams:
                dvdl.append( get_stdti_dvdl(parm,parmfilename,basename,rst,nc,lam,timask1,timask2,scmask1,scmask2) )
        deq = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)


        print "# %4s %12.4f"%("deq",deq)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("# %4s %12.4f\n"%("deq",deq))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()
        
            
    if has_sc:
        o        = tiobj.scti
        parm     = o.tiparm
        dirname  = os.path.join(o.basename,"template")
        basename = "test.%s"%(o.basename)
        parmfilename = os.path.join(dirname,"%s.parm7"%(o.basename))
        timask1  = o.merged_timask1
        timask2  = o.merged_timask2
        scmask1  = o.merged_scmask1
        scmask2  = o.merged_scmask2
        rst =  basename+".rst7"
        nc  =  basename+".nc"
        parmed.tools.writeCoordinates(parm, rst).execute()
        parmed.tools.writeCoordinates(parm, nc).execute()

        dvdl = []
        for lam in lams:
            dvdl.append( get_stdti_dvdl(parm,parmfilename,basename,rst,nc,lam,timask1,timask2,scmask1,scmask2) )
        sc = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)
    
        print "# %4s %12.4f"%("sc",sc)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("\n# %4s %12.4f\n"%("sc",sc))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()

    if has_req:
        o        = tiobj.reqti
        parm     = o.tiparm
        dirname  = os.path.join(o.basename,"template")
        basename = "test.%s"%(o.basename)
        parmfilename = os.path.join(dirname,"%s.parm7"%(o.basename))
        timask1  = o.merged_timask1
        timask2  = o.merged_timask2
        scmask1  = o.merged_scmask1
        scmask2  = o.merged_scmask2
        rst =  basename+".rst7"
        nc  =  basename+".nc"
        parmed.tools.writeCoordinates(parm, rst).execute()
        parmed.tools.writeCoordinates(parm, nc).execute()

        dvdl = []
        for lam in lams:
            dvdl.append( get_stdti_dvdl(parm,parmfilename,basename,rst,nc,lam,timask1,timask2,scmask1,scmask2) )
        req = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)
    
        print "# %4s %12.4f"%("req",req)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("\n# %4s %12.4f\n"%("req",req))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()
    print "# %4s %12.4f"%("sum",deq+sc+req)

    fh.close()
        


def get_piti_dvdl(parm0,parm1,parmfile,basename,rst,nc,lam,scmask1,scmask2):
    import subprocess
    import copy
    import os.path

    
    mdin = Mdin( base = basename, inpcrd = rst )
    mdin.RST7 = "restart"
    mdin.PARM7 = parmfile
    mdin.INPNC = basename+".nc"
    mdin.exe = "${AMBERHOME}/bin/sander"
    mdin.Set_NVT()
    mdin.Set_Restart(False)
    mdin.Set_PrintFreq(0)
    mdin.cntrl["nstlim"] = 0
    mdin.cntrl["ntpr"] = 1
    mdin.cntrl["ntf"]  = 1
    mdin.cntrl["cut"]  = 10.

    mdin.cntrl["imin"] = 6 
    mdin.cntrl["ntmin"] = 1

    is_softcore = False
    sel1 = GetSelectedAtomIndices(parm0,scmask1)
    sel2 = GetSelectedAtomIndices(parm0,scmask2)
    if len(sel1+sel2) > 0:
        is_softcore = True
    
    if is_softcore:
        mdin.cntrl["vdwmodel"] = 2
        mdin.ewald["chngmask"] = 0
        
    always_use_sander = True
    
    if always_use_sander or not is_softcore:
        mdin.WriteMdin()
        cmd = "%s %s"%(mdin.exe,mdin.CmdString())
        #print cmd
        subprocess.call(cmd,shell=True)
        dedparm = mdin.INPNC + ".dedparm"
        if not is_softcore:
            dvdl = nonsoftcore_tigradient(parm0,parm1,lam,dedparm)
        else:
            dvdl = softcore_tigradient(parm0,parm1,lam,dedparm,scmask1,scmask2)

        
    return dvdl[0]



class bonded_term_collector(object):
    def __init__(self):
        self.tiparm_type_idx = None
        self.org_parm_objs = []
        


def calc_piti_dvdls(tiobj):
    import os.path


    deq=0.
    has_deq = False
    try:
        if tiobj.deqti is not None:
            has_deq = True
    except:
        pass
    
    sc=0.
    has_sc  = False
    try:
        if tiobj.scti is not None:
            has_sc = True
    except:
        pass
    
    req=0.
    has_req = False
    try:
        if tiobj.reqti is not None:
            has_req = True
    except:
        pass


    
    if has_deq:
        o        = tiobj.deqti
        lams     = list(sorted( lam for lam in tiobj.deq_parms ))
        nlam     = len(lams)
        parm0    = tiobj.deq_parms[lams[0]] 
        parm1    = tiobj.deq_parms[lams[-1]]
        basename = o.basename
        dirname  = os.path.join(o.basename,"template")
        rebase = "test.%s"%(o.basename)

        rst =  rebase+".rst7"
        nc  =  rebase+".nc"
        parmed.tools.writeCoordinates(parm0, rst).execute()
        parmed.tools.writeCoordinates(parm0, nc).execute()


        dvdl = []
        for lam in lams:
            parm = os.path.join(dirname,"%s_%.8f.parm7"%(o.basename,lam))
            dvdl.append( get_piti_dvdl(parm0,parm1,parm,rebase,rst,nc,lam,"","") )
        deq = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)

        print "# %4s %12.4f"%("deq",deq)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("# %4s %12.4f\n"%("deq",deq))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()

        
    if has_sc:
        o        = tiobj.scti
        lams     = list(sorted( lam for lam in tiobj.sc_parms ))
        nlam     = len(lams)
        parm0    = tiobj.sc_parms[lams[0]] 
        parm1    = tiobj.sc_parms[lams[-1]] 
        basename = o.basename
        dirname  = os.path.join(o.basename,"template")
        rebase = "test.%s"%(o.basename)
        rst =  rebase+".rst7"
        nc  =  rebase+".nc"
        parmed.tools.writeCoordinates(parm0, rst).execute()
        parmed.tools.writeCoordinates(parm0, nc).execute()
        
        dvdl = []
        scmask1 = tiobj.scti.merged_timask1
        scmask2 = tiobj.scti.merged_timask2
#        print "scmask1: %s"%(scmask1)
#        print "scmask2: %s"%(scmask2)
        for lam in lams:
            parm = os.path.join(dirname,"%s_%.8f.parm7"%(basename,lam))
            dvdl.append( get_piti_dvdl(parm0,parm1,parm,rebase,rst,nc,lam,scmask1,scmask2) )
        sc = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)

        print "# %4s %12.4f"%("sc",sc)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("\n# %4s %12.4f\n"%("sc",sc))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()
        
    if has_req:
        o        = tiobj.reqti
        lams     = list(sorted( lam for lam in tiobj.req_parms ))
        nlam     = len(lams)
        parm0    = tiobj.req_parms[lams[0]] 
        parm1    = tiobj.req_parms[lams[-1]] 
        basename = o.basename
        dirname  = os.path.join(o.basename,"template")
        rebase = "test.%s"%(o.basename)
        rst =  rebase+".rst7"
        nc  =  rebase+".nc"
        parmed.tools.writeCoordinates(parm0, rst).execute()
        parmed.tools.writeCoordinates(parm0, nc).execute()

        dvdl = []
        for lam in lams:
            parm = os.path.join(dirname,"%s_%.8f.parm7"%(basename,lam))
            dvdl.append( get_piti_dvdl(parm0,parm1,parm,rebase,rst,nc,lam,"","") )
        req = ( sum(dvdl)-0.5*(dvdl[0]+dvdl[-1]) ) / (nlam-1.)
    
        print "# %4s %12.4f"%("req",req)
        fh = file("dvdl.%s.dat"%(o.basename),"w")
        fh.write("\n# %4s %12.4f\n"%("req",req))
        for ilam in range(len(lams)):
            fh.write(" %.8f %12.4f\n"%(lams[ilam],dvdl[ilam]))
        fh.close()
    print "# %4s %12.4f"%("sum",deq+sc+req)

    fh.close()
        
##########################################################
##########################################################
# PARMED-LIKE TI-MERGE UTILITY
##########################################################
##########################################################
    
class stdtimerge(object):
    def __init__(self,parm,molmask1,molmask2,timask1,timask2, keep_bondparms_from_mol=1):
        import numpy as np
        
        self.tol = 0.01
        self.org_parm = CopyParm( parm )
        
        natom = len(self.org_parm.atoms)
        
        self.molmask1 = molmask1
        self.molsel1  = GetSelectedAtomIndices( self.org_parm, self.molmask1 )
        self.molbit1  = [ 1 if a in self.molsel1 else 0 for a in range(natom) ]

        self.molmask2 = molmask2
        self.molsel2  = GetSelectedAtomIndices( self.org_parm, self.molmask2 )
        self.molbit2  = [ 1 if a in self.molsel2 else 0 for a in range(natom) ]

        self.timask1  = timask1
        self.tisel1   = GetSelectedAtomIndices( self.org_parm, self.timask1 )
        self.tibit1   = [ 1 if a in self.tisel1  else 0 for a in range(natom) ]

        self.timask2  = timask2
        self.tisel2   = GetSelectedAtomIndices( self.org_parm, self.timask2 )
        self.tibit2   = [ 1 if a in self.tisel2  else 0 for a in range(natom) ]

        #print "timerge '%s' '%s' '%s' '%s'"%(self.molmask1,self.molmask2,self.timask1,self.timask2)
        

        if keep_bondparms_from_mol == 2:
            molbitRem  = self.molbit1
            tibitRem   = self.tibit1
            molbitKeep = self.molbit2
            tibitKeep  = self.tibit2
        else:
            molbitRem  = self.molbit2
            tibitRem   = self.tibit2
            molbitKeep = self.molbit1
            tibitKeep  = self.tibit1
            
        
        for i in range(natom):
            if i in self.tisel1 and i not in self.molsel1:
                raise Exception('timask1 must be a subset of molmask1')
            if i in self.tisel2 and i not in self.molsel2:
                raise Exception('timask2 must be a subset of molmask2')

        # Hmph!? How is this a restriction, exactly??
        if len(self.molsel2) > 0:
            if (( max(self.molsel1)+1 ) != min(self.molsel2)) \
               and (( max(self.molsel2)+1 ) != min(self.molsel1)):
                raise Exception('molmask1 and molmask2 must be adjacent in the topology')

        self.tiparm = CopyParm( self.org_parm )

        if len(self.molsel2) == 0:
            if set(self.tisel1) == set(self.molsel1):
                for i in range(natom):
                    self.tiparm.atoms[i].original_index = [i]
                    self.org_parm.atoms[i].original_index = [i]
                    self.tiparm.atoms[i].modified_index = i
                    self.org_parm.atoms[i].modified_index = i
                return
            
        coordinates = self.tiparm.coordinates
        if coordinates is None:
            raise Exception('Failure copying coordinates to self.tiparm')

        #
        # We will (eventually) remove all atoms from molmask1 that are NOT in scmask1
        #
        atoms_to_remove = []
        remove_map  = [0 for i in range(natom)] # key: org idx, value: sft idx
        new_atm_idx = 0
        for i in range(natom):
            self.tiparm.atoms[i].original_index = [i]
            self.org_parm.atoms[i].original_index = [i]
            if molbitRem[i] == 1 and tibitRem[i] == 0:
                atoms_to_remove.append(i)
                self.tiparm.atoms[i].modified_index = None
                self.org_parm.atoms[i].modified_index = None
            else:
                remove_map[i] = new_atm_idx
                self.tiparm.atoms[i].modified_index = new_atm_idx
                self.org_parm.atoms[i].modified_index = new_atm_idx
                new_atm_idx += 1

        remove_str = ListToSelection( atoms_to_remove )

        #
        # These will be the new atom idx of the softcore atoms
        #
        new_sc_atm1_int = []
        new_sc_atm2_int = []
        for i in range(natom):
            if self.tibit1[i]:
                new_sc_atm1_int.append(remove_map[i])
            elif self.tibit2[i]:
                new_sc_atm2_int.append(remove_map[i])
#        self.merged_timask1 = ListToSelection(new_sc_atm1_int)
#        self.merged_timask2 = ListToSelection(new_sc_atm2_int)
        

        
        mol1common = []
        for i in range(natom):
            if self.molbit1[i] == 1 and self.tibit1[i] == 0:
                mol1common.append(i)

        mol2common = []
        for i in range(natom):
            if self.molbit2[i] == 1 and self.tibit2[i] == 0:
                mol2common.append(i)

        if len(mol1common) != len(mol2common):
            raise Exception('The number of nonsoftcore atoms in '
                            'molmask1 and molmask2 must be the same.')

        mol2common_sort = []
        # reorder mol2common so that it matches mol1common
        for i in range(len(mol1common)):
            atm_i = mol1common[i]
            for j in range(len(mol2common)):
                atm_j = mol2common[j]
                diff = coordinates[atm_i]-coordinates[atm_j]
                if (np.abs(diff) < self.tol).sum() == 3:
                    mol2common_sort.append(atm_j)
                    break

        mol2common = mol2common_sort

        # check again if we didn't match all coords
        if len(mol1common) != len(mol2common):
            raise Exception('The number of nonsoftcore atoms in molmask1 '
                            'and molmask2 must be the same. Check the '
                            'masks. If these look correct try using a '
                            'larger tolerance.')

        #print len(mol1common)
        
        for i in range(len(mol1common)):
            atm_i = mol1common[i]
            atm_j = mol2common[i]
            #print "Common %6i %6i"%(atm_i,atm_j)
            if self.tiparm.atoms[atm_i].modified_index is None:
                self.tiparm.atoms[atm_i].modified_index = self.tiparm.atoms[atm_j].modified_index
                self.org_parm.atoms[atm_i].modified_index = self.org_parm.atoms[atm_j].modified_index
                oidx = self.tiparm.atoms[atm_i].original_index+self.tiparm.atoms[atm_j].original_index
                self.tiparm.atoms[atm_i].original_index = oidx
                self.tiparm.atoms[atm_j].original_index = oidx
            else:
                self.tiparm.atoms[atm_j].modified_index = self.tiparm.atoms[atm_i].modified_index
                oidx = self.tiparm.atoms[atm_i].original_index+self.tiparm.atoms[atm_j].original_index
                self.tiparm.atoms[atm_i].original_index = oidx
                self.tiparm.atoms[atm_j].original_index = oidx
                
            diff = coordinates[atm_i]-coordinates[atm_j]
            if (np.abs(diff) > self.tol).any():
                raise Exception('Common (nonsoftcore) atoms must have the ' # pragma: no cover
                                'same coordinates.')

        # We need to (temporarily) keep any common atom, i, if there is
        # a bond, angle, or torsion from the sc-region that touches it
        # ...so we can reset the parameters
        keep_mask = [0 for i in range(natom)]
        for i in range(natom):
            if molbitRem[i]:
                for j in range(natom):
                    if tibitRem[j]:
                        atm1 = self.tiparm.atoms[i]
                        atm2 = self.tiparm.atoms[j]

                        # if atom i is common and atom 2 is softcore, then
                        # don't throw away atom i yet
                        if (atm1 in atm2.bond_partners or
                            atm1 in atm2.angle_partners or
                            atm1 in atm2.dihedral_partners):
                            keep_mask[i] = 1
                            break

        for j in range(natom):
            if keep_mask[j] == 1 and tibitRem[j] == 0:
                atm = self.tiparm.atoms[j]
                if keep_bondparms_from_mol == 2:
                    idx = mol2common[mol1common.index(j)]
                else:
                    idx = mol1common[mol2common.index(j)]
                atm_new = self.tiparm.atoms[idx]
                    

                # What is happening here???
                for k in range(natom):
                    if tibitRem[k]:
                        atm2 = self.tiparm.atoms[k]
                        # update partners -- the exclusion list will be updated
                        # when the file is written out
                        if atm in atm2.bond_partners:
                            atm._bond_partners.remove(atm2)
                            atm2._bond_partners.remove(atm)
                            atm2.bond_to(atm_new)

                        if atm in atm2.angle_partners:
                            atm._angle_partners.remove(atm2)
                            atm2._angle_partners.remove(atm)
                            atm2.angle_to(atm_new)

                        if atm in atm2.dihedral_partners:
                            atm._dihedral_partners.remove(atm2)
                            atm2._dihedral_partners.remove(atm)
                            atm2.dihedral_to(atm_new)

                        # Now go through each array re-indexing the atoms
                        # Check to make sure that this is a bond/angle/dihed
                        # involving the common atom j and the softcore atom k

                        for bond in self.tiparm.bonds:
                            if (bond.atom1.idx == j and bond.atom2.idx == k):
                                bond.atom1 = atm_new
                            elif (bond.atom2.idx == j and bond.atom1.idx == k):
                                bond.atom2 = atm_new

                        for angle in self.tiparm.angles:
                            if angle.atom1.idx == j:
                                if angle.atom2.idx == k or angle.atom3.idx == k:
                                    angle.atom1 = atm_new
                            elif angle.atom2.idx == j:
                                if angle.atom1.idx == k or angle.atom3.idx == k:
                                    angle.atom2 = atm_new
                            elif angle.atom3.idx == j:
                                if angle.atom1.idx == k or angle.atom2.idx == k:
                                    angle.atom3 = atm_new

                        for dihed in self.tiparm.dihedrals:
                            if dihed.atom1.idx == j:
                                if (dihed.atom2.idx == k or dihed.atom3.idx == k
                                        or dihed.atom4.idx == k):
                                    dihed.atom1 = atm_new
                            elif dihed.atom2.idx == j:
                                if (dihed.atom1.idx == k or dihed.atom3.idx == k
                                        or dihed.atom4.idx == k):
                                    dihed.atom2 = atm_new
                            elif dihed.atom3.idx == j:
                                if (dihed.atom1.idx == k or dihed.atom2.idx == k
                                        or dihed.atom4.idx == k):
                                    dihed.atom3 = atm_new
                            elif dihed.atom4.idx == j:
                                if (dihed.atom1.idx == k or dihed.atom2.idx == k
                                        or dihed.atom3.idx == k):
                                    dihed.atom4 = atm_new

                        for imp in self.tiparm.impropers:
                            if imp.atom1.idx == j:
                                if (imp.atom2.idx == k or imp.atom3.idx == k
                                        or imp.atom4.idx == k):
                                    imp.atom1 = atm_new
                            elif imp.atom2.idx == j:
                                if (imp.atom1.idx == k or imp.atom3.idx == k
                                        or imp.atom4.idx == k):
                                    imp.atom2 = atm_new
                            elif imp.atom3.idx == j:
                                if (imp.atom1.idx == k or imp.atom2.idx == k
                                        or imp.atom4.idx == k):
                                    imp.atom3 = atm_new
                            elif imp.atom4.idx == j:
                                if (imp.atom1.idx == k or imp.atom2.idx == k
                                        or imp.atom3.idx == k):
                                    imp.atom4 = atm_new

                        for cmap in self.tiparm.cmaps:
                            if cmap.atom1.idx == j:
                                if (cmap.atom2.idx == k or cmap.atom3.idx == k
                                        or cmap.atom4.idx == k
                                        or cmap.atom5.idx == k):
                                    cmap.atom1 = atm_new
                            elif cmap.atom2.idx == j:
                                if (cmap.atom1.idx == k or cmap.atom3.idx == k
                                        or cmap.atom4.idx == k
                                        or cmap.atom5.idx == k):
                                    cmap.atom2 = atm_new
                            elif cmap.atom3.idx == j:
                                if (cmap.atom1.idx == k or cmap.atom2.idx == k
                                        or cmap.atom4.idx == k
                                        or cmap.atom5.idx == k):
                                    cmap.atom3 = atm_new
                            elif cmap.atom4.idx == j:
                                if (cmap.atom1.idx == k or cmap.atom2.idx == k
                                        or cmap.atom3.idx == k
                                        or cmap.atom5.idx == k):
                                    cmap.atom4 = atm_new
                            elif cmap.atom5.idx == j:
                                if (cmap.atom1.idx == k or cmap.atom2.idx == k
                                        or cmap.atom3.idx == k
                                        or cmap.atom4.idx == k):
                                    cmap.atom5 = atm_new
                            
        self.tiparm.atoms.changed = True
        if len(atoms_to_remove) > 0:
            self.tiparm.strip(remove_str)
            
        
    def save(self,basename,subdir=None):
        import os.path
        self.basename = basename
        if subdir is not None:
            basename = os.path.join(subdir,basename)
        SaveParm( self.tiparm, "%s.parm7"%(basename) )
        parmed.tools.writeCoordinates( self.tiparm, "%s.rst7"%(basename) ).execute()

        
    def remask(self,mask):
        oats = GetSelectedAtomIndices( self.org_parm, mask )
        mats = []
        for a in oats:
            mats.append( self.org_parm.atoms[a].modified_index )
            if a not in self.tiparm.atoms[ mats[-1] ].original_index:
                raise Exception("Could not find original index %i in merged parm atom %i while processing mask '%s'"%(a+1,mats[-1]+1,mask))
        return ListToSelection( mats )

    


    
def ti_init( self,               \
             parmfile, rstfile,  \
             molmask1, molmask2, \
             timask1,  timask2,  \
             scmask1,  scmask2   ):
    
    self.tol = 0.01
    self.parmfile = parmfile
    self.rstfile  = rstfile
    self.true_dihedral_deletion = True
    self.org_parm = OpenParm( parmfile, xyz=rstfile )
    self.sft_parm = None
    self.deq_parm = None
    self.req_parm = None

    
    if False:
        self.org_parm.parm_data["CHARGE"] = [ 0. ] * len(self.org_parm.parm_data["CHARGE"])
        for i in range(len(self.org_parm.atoms)):
            self.org_parm.atoms[i].charge = 0.
    if False:
        self.org_parm.parm_data["BOND_FORCE_CONSTANT"] = [ 0. ] * len(self.org_parm.parm_data["BOND_FORCE_CONSTANT"])
        for i in range(len(self.org_parm.bond_types)):
            self.org_parm.bond_types[i].k = 0
    if False:
        self.org_parm.parm_data["ANGLE_FORCE_CONSTANT"] = [ 0. ] * len(self.org_parm.parm_data["ANGLE_FORCE_CONSTANT"])
        for i in range(len(self.org_parm.angle_types)):
            self.org_parm.angle_types[i].k = 0
    if False:
        self.org_parm.parm_data["DIHEDRAL_FORCE_CONSTANT"] = [ 0. ] * len(self.org_parm.parm_data["DIHEDRAL_FORCE_CONSTANT"])
        for i in range(len(self.org_parm.dihedral_types)):
            self.org_parm.dihedral_types[i].phi_k = 0
    if False:
        self.org_parm.parm_data["LENNARD_JONES_ACOEF"] = [ 0. ] * len(self.org_parm.parm_data["LENNARD_JONES_ACOEF"])
        self.org_parm.parm_data["LENNARD_JONES_BCOEF"] = [ 0. ] * len(self.org_parm.parm_data["LENNARD_JONES_BCOEF"])
        self.org_parm.LJ_radius = [0.]*len(self.org_parm.LJ_radius)
        self.org_parm.LJ_depth  = [0.]*len(self.org_parm.LJ_depth)
#        for i in range(len(self.org_parm.atoms)):
#            self.org_parm.atoms[i].rmin = 0.
#            self.org_parm.atoms[i]. = 0.
        self.org_parm.recalculate_LJ()

    
            


            
    if self.org_parm.coordinates is None:
        raise Exception('Load coordinates before merging topology')
    
    natom = len(self.org_parm.atoms)

    self.molmask1 = molmask1
    self.molsel1  = GetSelectedAtomIndices( self.org_parm, self.molmask1 )
    self.molbit1  = [ 1 if a in self.molsel1 else 0 for a in range(natom) ]

    self.molmask2 = molmask2
    self.molsel2  = GetSelectedAtomIndices( self.org_parm, self.molmask2 )
    self.molbit2  = [ 1 if a in self.molsel2 else 0 for a in range(natom) ]

    self.timask1  = timask1
    self.tisel1   = GetSelectedAtomIndices( self.org_parm, self.timask1 )
    self.tibit1   = [ 1 if a in self.tisel1  else 0 for a in range(natom) ]

    self.timask2  = timask2
    self.tisel2   = GetSelectedAtomIndices( self.org_parm, self.timask2 )
    self.tibit2   = [ 1 if a in self.tisel2  else 0 for a in range(natom) ]

    self.scmask1  = scmask1
    self.scsel1   = GetSelectedAtomIndices( self.org_parm, self.scmask1 )
    self.scbit1   = [ 1 if a in self.scsel1  else 0 for a in range(natom) ]

    self.scmask2  = scmask2
    self.scsel2   = GetSelectedAtomIndices( self.org_parm, self.scmask2 )
    self.scbit2   = [ 1 if a in self.scsel2  else 0 for a in range(natom) ]

        
    for i in range(natom):
        if i in self.tisel1 and i not in self.molsel1:
            raise Exception('timask1 must be a subset of molmask1')
        if i in self.scsel1 and i not in self.tisel1:
            raise Exception('scmask1 must be a subset of timask1')
        if i in self.tisel2 and i not in self.molsel2:
            raise Exception('timask2 must be a subset of molmask2')
        if i in self.scsel2 and i not in self.tisel2:
            raise Exception('scmask2 must be a subset of timask2')

    # Hmph!? How is this a restriction, exactly??
    if len(self.molsel2) > 0:
        if (( max(self.molsel1)+1 ) != min(self.molsel2)) \
           and (( max(self.molsel2)+1 ) != min(self.molsel1)):
            raise Exception('molmask1 and molmask2 must be adjacent in the topology')

    if False:
        ti1 = []
        for i in self.tisel1:
            if i not in self.scsel1:
                ti1.append( i )
        ti2 = []
        for i in self.tisel2:
            if i not in self.scsel2:
                ti2.append( i )
        ti12 = ti2
        for i in range(len(self.org_parm.dihedrals)):
            d = self.org_parm.dihedrals[i]
            if d.atom1.idx in ti12 and \
               d.atom2.idx in ti12 and \
               d.atom3.idx in ti12 and \
               d.atom4.idx in ti12:
                tidx = d.type.idx
                d.type.phi_k = 0.
                self.org_parm.parm_data["DIHEDRAL_FORCE_CONSTANT"][tidx] = 0.




##########################################################
##########################################################
# STANDARD TI PATHWAY
##########################################################
##########################################################
    
class sctisetup(object):
    
    ###################################################################################
    ###################################################################################
    def __init__( self,               \
                  parmfile, rstfile,  \
                  molmask1, molmask2, \
                  timask1,  timask2,  \
                  scmask1,  scmask2, \
                  refit="" ):
        #print "scmask1=",scmask1
        self.refit=refit
        ti_init(self,parmfile,rstfile,molmask1,molmask2,timask1,timask2,scmask1,scmask2)
        self.true_dihderal_deletion = False # not possible otherwise
        self.prepare_common_atom_maps()
        self.pathway = "std"

    def test(self):
        calc_stdti_dvdls(self)

    def execute(self,lams,compinfo,sclams=None):
        #from parmutils import ListToSelection as l2s
        #from parmutils import GetSelectedAtomIndices as s2l
        l2s=ListToSelection
        s2l=GetSelectedAtomIndices

        
        if len(self.scsel1 + self.scsel2) > 0 and len(self.scsel1) == 0:
            pass
        else:
            basename=self.pathway + "deq"
            print "prepare_decharge"
            self.prepare_decharge(lams,save=True,basename=basename)
            o = self.deqti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(lams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))

        if len(self.scsel1 + self.scsel2) > 0:
            if sclams is None:
                sclams = lams
            basename=self.pathway + "sc"
            print "prepare_softcore"
            self.prepare_softcore(sclams,save=True,basename=basename)
            o = self.scti
            c1 = s2l( o.tiparm, o.merged_scmask1 )
            c2 = s2l( o.tiparm, o.merged_scmask2 )
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(sclams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))

        if len(self.scsel2) > 0:
            basename=self.pathway + "req"
            print "prepare_recharge"
            self.prepare_recharge(lams,save=True,basename=basename)
            o = self.reqti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(lams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))
        
    def prepare_common_atom_maps(self):
        from collections import defaultdict as ddict
        
        self.o1_to_o2 = ddict( int )
        self.o2_to_o1 = ddict( int )
        self.ti_o1_to_o2 = ddict( int )
        self.ti_o2_to_o1 = ddict( int )


        if len( self.molsel2 ) > 0:
        
            common = stdtimerge( self.org_parm,
                                 self.molmask1,self.molmask2,
                                 self.scmask1,self.scmask2,
                                 keep_bondparms_from_mol = 2 )
        
            for a in common.tiparm.atoms:
                #print "%6i %s"%(a.idx,str(a.original_index))
                if len(a.original_index) > 1:
                    first = min(a.original_index)
                    second = max(a.original_index)
                    self.o1_to_o2[ first ] = second
                    self.o2_to_o1[ second ] = first
                    if first in self.tisel1 and second in self.tisel2:
                        self.ti_o1_to_o2[ first ] = second
                        self.ti_o2_to_o1[ second ] = first
                    else:
                        continue

        if len(self.refit) > 0:
            refitq=False
            refith=False
            if "H" in self.refit:
                refith=True
            elif "Q" in self.refit:
                refitq=True

        self.ti1_charges = []
        for k in sorted(self.tisel1):
            q = self.org_parm.atoms[k].charge
            if k in self.scsel1:
                q=0
            self.ti1_charges.append(q)

        if len(self.refit) > 0:
            if len(self.tisel1) - len(self.scsel1) > 0 and len(self.scsel1) > 0:
                print "refitting ti1_charges"
                self.ti1_charges = RefitCharges( self.org_parm, self.timask1, self.scmask1, equivH=refith, equivQ=refitq )


                
        self.ti2_charges = []
        for k in sorted(self.tisel2):
            q = self.org_parm.atoms[k].charge
            if k in self.scsel2:
                q=0
            self.ti2_charges.append(q)
        if len(self.refit) > 0:
            if len(self.tisel2) - len(self.scsel2) > 0 and  len(self.scsel2) > 0:
                print "refitting ti2_charges"
                self.ti2_charges = RefitCharges( self.org_parm, self.timask2, self.scmask2, equivH=refith, equivQ=refitq )        
            
                #print "%6i <-> %6i"%(first,second)

                
    ###################################################################################
    ###################################################################################
    def prepare_softcore(self,lams,save=True,basename="stdsc"):

        if len(self.molsel2) > 0:
            self.scti = stdtimerge( self.org_parm,
                                    self.molmask1,self.molmask2,
                                    self.timask1,self.timask2,
                                    keep_bondparms_from_mol = 2 )

            self.scti.merged_timask1 = self.scti.remask( self.timask1 )
            self.scti.merged_timask2 = self.scti.remask( self.timask2 )
            self.scti.merged_scmask1 = self.scti.remask( self.scmask1 )
            self.scti.merged_scmask2 = self.scti.remask( self.scmask2 )
            
            sc_atoms = GetSelectedAtomIndices(
                self.scti.tiparm,
                "(%s)|(%s)"%(self.scti.merged_scmask1,self.scti.merged_scmask2) )

            ###
            print "ti1_charges",self.ti1_charges
            print "ti2_charges",self.ti2_charges
            ti_atoms = GetSelectedAtomIndices( self.scti.tiparm, self.scti.merged_timask1 )
            for iat in range(len(ti_atoms)):
                a = ti_atoms[iat]
                self.scti.tiparm.atoms[a].charge = self.ti1_charges[iat]
                self.scti.tiparm.parm_data["CHARGE"][a] = self.ti1_charges[iat]
                
            ti_atoms = GetSelectedAtomIndices( self.scti.tiparm, self.scti.merged_timask2 )
            for iat in range(len(ti_atoms)):
                a = ti_atoms[iat]
                self.scti.tiparm.atoms[a].charge = self.ti2_charges[iat]
                self.scti.tiparm.parm_data["CHARGE"][a] = self.ti2_charges[iat]
            ###
        

            
        else:
            self.scti = stdtimerge( self.org_parm,
                                    self.molmask1,self.molmask2,
                                    self.timask1,self.timask2,
                                    keep_bondparms_from_mol = 1 )
            
            self.scti.merged_timask1 = self.scti.remask( self.timask1 )
            self.scti.merged_timask2 = ""
            self.scti.merged_scmask1 = self.scti.remask( self.scmask1 )
            self.scti.merged_scmask2 = ""
            
            sc_atoms = GetSelectedAtomIndices(
                self.scti.tiparm,
                "%s"%( self.scti.merged_scmask1 ) )

            
            
        for a in sc_atoms:
            self.scti.tiparm.atoms[a].charge = 0.
            self.scti.tiparm.parm_data["CHARGE"][a] = 0.
            
        self.scti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.scti.save( basename, subdir=dirname )
        
        #print "scstep: icfe=1, ifsc=1, timask1='%s', timask2='%s', scmask1='%s', scmask2='%s', crgmask=''"%(self.scti.merged_timask1,self.scti.merged_timask2,self.scti.merged_scmask1,self.scti.merged_scmask2)


            
    
    ###################################################################################
    ###################################################################################
    def prepare_decharge(self,lams,save=True,basename="stddeq"):
        import numpy as np
        
        natom = len(self.org_parm.atoms)
        self.deq_parm = CopyParm( self.org_parm )
        if self.deq_parm.coordinates is None:
            raise Exception('Failure copying coordinates to self.deq_parm')

        premol_init_common_idx  = 0
        premol_last_common_idx  = min( self.molsel1 + self.molsel2 ) - 1
        postmol_init_common_idx = max( self.molsel1 + self.molsel2 ) + 1
        postmol_last_common_idx = natom - 1

        mol1_natom = len(self.molsel1)

        for i,a in enumerate(self.deq_parm.atoms):
            a.master_index = i
        mol = Extract( self.deq_parm, self.molmask1 )
        postmask = "@%i-%i"%( postmol_init_common_idx+1, postmol_last_common_idx+1 )
        postmol  = Extract( self.deq_parm, postmask )
        if premol_last_common_idx >= 0:
            premask = "@%i-%i"%( premol_init_common_idx+1, premol_last_common_idx+1 )
            premol  = Extract( self.deq_parm, premask )
            
            master_array = GetSelectedAtomIndices( self.deq_parm,premask ) \
                           + GetSelectedAtomIndices( self.deq_parm,self.molmask1 ) \
                           + GetSelectedAtomIndices( self.deq_parm,self.molmask1 ) \
                           + GetSelectedAtomIndices( self.deq_parm,postmask )
            master_array = [ self.deq_parm.atoms[idx].master_index for idx in master_array ]
            self.deq_parm = Join(
                premol, Join( mol, Join( mol, postmol ) ) )
        else:
            master_array = GetSelectedAtomIndices( self.deq_parm,self.molmask1 ) \
                           + GetSelectedAtomIndices( self.deq_parm,self.molmask1 ) \
                           + GetSelectedAtomIndices( self.deq_parm,postmask )
            master_array = [ self.deq_parm.atoms[idx].master_index for idx in master_array ]
            self.deq_parm = Join( mol, Join( mol, postmol ) )

        for idx,a in zip( master_array, self.deq_parm.atoms):
            a.master_index = idx

        old_mol1_atom1 = min(self.molsel1)
        new_mol1_atom1 = premol_last_common_idx + 1
        new_mol2_atom1 = new_mol1_atom1 + mol1_natom
        new_molmask1   = "@%i-%i"%( new_mol1_atom1+1, new_mol1_atom1+mol1_natom )
        new_molmask2   = "@%i-%i"%( new_mol2_atom1+1, new_mol2_atom1+mol1_natom )
        new_tisel1     = [ i-old_mol1_atom1+new_mol1_atom1 for i in self.tisel1 ]
        new_tisel2     = [ i-old_mol1_atom1+new_mol2_atom1 for i in self.tisel1 ]
        new_timask1    = ListToSelection( new_tisel1 )
        new_timask2    = ListToSelection( new_tisel2 )
        new_scsel1     = [ i-old_mol1_atom1+new_mol1_atom1 for i in self.scsel1 ]
        new_scsel2     = [ i-old_mol1_atom1+new_mol2_atom1 for i in self.scsel1 ]
        new_scmask1    = ListToSelection( new_scsel1 )
        new_scmask2    = ListToSelection( new_scsel2 )

        ti = stdtimerge( self.deq_parm,
                         new_molmask1, new_molmask2,
                         new_timask1,  new_timask2,
                         keep_bondparms_from_mol = 1 )
        
        ti.merged_timask1 = ti.remask( new_timask1 )
        ti.merged_timask2 = ti.remask( new_timask2 )
        ti.merged_scmask1 = ""
        ti.merged_scmask2 = ""

        sc_atoms = GetSelectedAtomIndices( ti.tiparm, ti.remask( new_scmask2 ) )
        for a in sc_atoms:
            ti.tiparm.atoms[a].charge = 0.
            ti.tiparm.parm_data["CHARGE"][a] = 0.

        ###
        ti_atoms = GetSelectedAtomIndices( ti.tiparm, ti.remask( new_timask2 ) )
        for iat in range(len(ti_atoms)):
            a=ti_atoms[iat]
            ti.tiparm.atoms[a].charge = self.ti1_charges[iat]
            ti.tiparm.parm_data["CHARGE"][a] = self.ti1_charges[iat]
        ###
        

        tiorg = CopyParm( self.org_parm )
        for iat,a in enumerate(tiorg.atoms):
            a.master_index = iat
            a.original_index = [iat]
            a.modified_index = []
        for iat,a in enumerate(ti.tiparm.atoms):
            #print iat,a.original_index
            for i in range(len(a.original_index)):
                o = a.original_index[i]
                #print iat,o,len(self.deq_parm.atoms)
                m = self.deq_parm.atoms[o].master_index
                a.original_index[i] = m
                tiorg.atoms[m].modified_index.append(iat)
        for a in tiorg.atoms:
            if len(a.modified_index) == 0:
                a.modified_index = None
            elif len(a.modified_index) == 1:
                a.modified_index = a.modified_index[0]
        ti.org_parm = tiorg

        #for i,a in enumerate(ti.org_parm.atoms):
        #    print "%6i %4s:%4s %s"%( i+1, a.residue.name, a.name, str(a.modified_index) )



        
        self.deqti = ti
        self.deqti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.deqti.save( basename, subdir=dirname )

        #print "deqstep: icfe=1, ifsc=0, timask1='%s', timask2='%s', scmask1='', scmask2='', crgmask=''"%(ti.merged_timask1,ti.merged_timask2)

        
    ###################################################################################
    ###################################################################################
    def prepare_recharge(self,lams,save=True,basename="stdreq"):
        import numpy as np
        
        natom = len(self.org_parm.atoms)
        self.req_parm = CopyParm( self.org_parm )
        if self.req_parm.coordinates is None:
            raise Exception('Failure copying coordinates to self.req_parm')

        premol_init_common_idx  = 0
        premol_last_common_idx  = min( self.molsel1 + self.molsel2 ) - 1
        postmol_init_common_idx = max( self.molsel1 + self.molsel2 ) + 1
        postmol_last_common_idx = natom - 1

        mol2_natom = len(self.molsel2)

        for i,a in enumerate(self.req_parm.atoms):
            a.master_index = i
        mol = Extract( self.req_parm, self.molmask2 )
        postmask = "@%i-%i"%( postmol_init_common_idx+1, postmol_last_common_idx+1 )
        postmol  = Extract( self.req_parm, postmask )
        if premol_last_common_idx >= 0:
            premask = "@%i-%i"%( premol_init_common_idx+1, premol_last_common_idx+1 )
            premol  = Extract( self.req_parm, premask )
            master_array = GetSelectedAtomIndices( self.req_parm,premask ) \
                           + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                           + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                           + GetSelectedAtomIndices( self.req_parm,postmask )
            master_array = [ self.req_parm.atoms[idx].master_index for idx in master_array ]
            
            #print premask
            self.req_parm = Join(
                premol, Join( mol, Join( mol, postmol ) ) )
        else:
            master_array = GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                           + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                           + GetSelectedAtomIndices( self.req_parm,postmask )
            master_array = [ self.req_parm.atoms[idx].master_index for idx in master_array ]
            self.req_parm = Join( mol, Join( mol, postmol ) )

        for idx,a in zip( master_array, self.req_parm.atoms):
            a.master_index = idx

        old_mol2_atom1 = min(self.molsel2)
        new_mol1_atom1 = premol_last_common_idx + 1
        new_mol2_atom1 = new_mol1_atom1 + mol2_natom
        new_molmask1   = "@%i-%i"%( new_mol1_atom1+1, new_mol1_atom1+mol2_natom )
        new_molmask2   = "@%i-%i"%( new_mol2_atom1+1, new_mol2_atom1+mol2_natom )
        new_tisel1     = [ i-old_mol2_atom1+new_mol1_atom1 for i in self.tisel2 ]
        new_tisel2     = [ i-old_mol2_atom1+new_mol2_atom1 for i in self.tisel2 ]


        
        # for idx in range(new_mol1_atom1):
        #     self.req_parm.atoms
        # for newidx,orgidx in zip(new_tisel1,self.tisel2):
        #     self.req_parm.atoms[newidx].master_index = orgidx
        # for newidx,orgidx in zip(new_tisel2,self.tisel2):
        #     self.req_parm.atoms[newidx].master_index = orgidx

        new_timask1    = ListToSelection( new_tisel1 )
        new_timask2    = ListToSelection( new_tisel2 )
        new_scsel1     = [ i-old_mol2_atom1+new_mol1_atom1 for i in self.scsel2 ]
        new_scsel2     = [ i-old_mol2_atom1+new_mol2_atom1 for i in self.scsel2 ]
        new_scmask1    = ListToSelection( new_scsel1 )
        new_scmask2    = ListToSelection( new_scsel2 )

        ti = stdtimerge( self.req_parm,
                         new_molmask1, new_molmask2,
                         new_timask1,  new_timask2,
                         keep_bondparms_from_mol = 2 )

        ti.merged_timask1 = ti.remask( new_timask1 )
        ti.merged_timask2 = ti.remask( new_timask2 )
        ti.merged_scmask1 = ""
        ti.merged_scmask2 = ""
        
        sc_atoms = GetSelectedAtomIndices( ti.tiparm, ti.remask( new_scmask1 ) )
        for a in sc_atoms:
            ti.tiparm.atoms[a].charge = 0.
            ti.tiparm.parm_data["CHARGE"][a] = 0.


        ###
        ti_atoms = GetSelectedAtomIndices( ti.tiparm, ti.remask( new_timask1 ) )
        for iat in range(len(ti_atoms)):
            a=ti_atoms[iat]
            ti.tiparm.atoms[a].charge = self.ti2_charges[iat]
            ti.tiparm.parm_data["CHARGE"][a] = self.ti2_charges[iat]
        ###
        

            

        tiorg = CopyParm( self.org_parm )
        for iat,a in enumerate(tiorg.atoms):
            a.master_index = iat
            a.original_index = [iat]
            a.modified_index = []
        for iat,a in enumerate(ti.tiparm.atoms):
            for i in range(len(a.original_index)):
                o = a.original_index[i]
                m = self.req_parm.atoms[o].master_index
                a.original_index[i] = m
                tiorg.atoms[m].modified_index.append(iat)
        for a in tiorg.atoms:
            if len(a.modified_index) == 0:
                a.modified_index = None
            elif len(a.modified_index) == 1:
                a.modified_index = a.modified_index[0]
        ti.org_parm = tiorg
        
        self.reqti = ti
        self.reqti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.reqti.save( basename, subdir=dirname )

        
        #print "reqstep: icfe=1, ifsc=0, timask1='%s', timask2='%s', scmask1='', scmask2='', crgmask=''"%(ti.merged_timask1,ti.merged_timask2)



   

def AtmStr(a):
    if isinstance(a,list):
        return "(%s)"%( ", ".join( [ AtmStr(x) for x in a ] ) )
    else:
        return "%4i:%4s@%-4s"%(a.residue.idx+1,a.residue.name,a.name)

def GetModifiedIndexs(tiobj,orgmask,imodregion=1):
    osele = GetSelectedAtomIndices( tiobj.org_parm, orgmask )
    msele = []
    for o1 in osele:
        oatm = tiobj.org_parm.atoms[ o1 ]
        midx = None
        if isinstance(oatm.modified_index,list):
            midx = oatm.modified_index[ imodregion ]
        else:
            midx = oatm.modified_index
        msele.append( midx )
    return msele

def GetModifiedMask( tiobj, orgmask, imodregion=1 ):
    msele = GetModifiedIndexs(tiobj,orgmask,imodregion)
    return "@%s"%( ",".join( [ str(x+1) for x in msele ] ) )


def Org2ModifiedIndexs(tiobj,osele,imodregion=1):
    msele = []
    for o1 in osele:
        oatm = tiobj.org_parm.atoms[ o1 ]
        midx = None
        if isinstance(oatm.modified_index,list):
            midx = oatm.modified_index[ imodregion ]
        else:
            midx = oatm.modified_index
        msele.append( midx )
    return msele


##########################################################
##########################################################
# SERIES TI PATHWAY (LIKE STD, BUT REMOVES SOFTCORE REGION
# BY HEAVY ATOM INSTEAD OF PARAMETER INTERPOLATION)
##########################################################
##########################################################
    
class sersctisetup(sctisetup):
    
    ###################################################################################
    ###################################################################################
    def __init__( self,               \
                  parmfile, rstfile,  \
                  molmask1, molmask2, \
                  timask1,  timask2,  \
                  scmask1,  scmask2,  \
                  deletedihedrals=True   ):
        super(sersctisetup,self).__init__(parmfile,rstfile,molmask1,molmask2,timask1,timask2,scmask1,scmask2)
        self.pathway = "ser"


    def test(self):
        calc_stdti_dvdls(self)


    def execute(self,lams,compinfo,sclams=None):
        import copy
        import os
        #from parmutils import ListToSelection as l2s
        #from parmutils import GetSelectedAtomIndices as s2l
        l2s=ListToSelection
        s2l=GetSelectedAtomIndices

        
        if len(self.scsel1 + self.scsel2) > 0 and len(self.scsel1) == 0:
            pass
        else:
            basename=self.pathway + "deq"
            print "prepare_decharge"
            self.prepare_decharge(lams,save=True,basename=basename)
            o = self.deqti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(lams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))

        if len(self.scsel1 + self.scsel2) > 0:
            if sclams is None:
                sclams = lams
            istage = 0
            if len(self.scsel1) > 0:
                
                p = CopyParm( self.org_parm )
                for ik,k in enumerate( sorted(self.tisel1) ):
                    p.atoms[k].charge = self.ti1_charges[ik]
                    p.parm_data["CHARGE"][k] = self.ti1_charges[ik]
                tisel1 = copy.deepcopy( self.tisel1 )
                scsel1 = copy.deepcopy( self.scsel1 )
                molsel2 = copy.deepcopy( self.molsel2 )
                tisel1 = ReindexFromDeletions( tisel1, molsel2 )
                scsel1 = ReindexFromDeletions( scsel1, molsel2 )
                p = Strip( p, self.molmask2 )
                fseries = find_atoms_to_delete( p, scsel1 )
                for iser in range(len(fseries)):
                    for pser in range(iser):
                        fseries[iser] = ReindexFromDeletions( fseries[iser], fseries[pser] )
                for iser in range(len(fseries)):
                    istage += 1
                    basename=self.pathway + "sc" + "%02i"%(istage)
                    print("scsel",[p.atoms[idx].name for idx in scsel1])
                    newdels = fseries[iser]
                    print([p.atoms[idx].name for idx in newdels])

                    newdelsel = ListToSelection(newdels)
                    newtisel = ListToSelection(tisel1)
                    newscsel = ListToSelection(scsel1)
                    dirname = os.path.join(basename,"template")
                    if not os.path.exists(dirname):
                        os.makedirs(dirname)
                    crgmask = l2s(scsel1)
                    noshakemask = newtisel
                    write_stdti_mdin_script(sclams,basename,compinfo,
                                            crgmask,
                                            noshakemask,
                                            newdelsel,"",
                                            newdelsel,"",
                                            is_gas=(self.org_parm.box is None))
                    SaveParm( p, os.path.join(dirname,basename+".parm7") )
                    parmed.tools.writeCoordinates( p, os.path.join(dirname,"%s.rst7"%(basename)) ).execute()
                    p = Strip( p, newdelsel )
                    tisel1  = ReindexFromDeletions( tisel1, newdels )
                    scsel1  = ReindexFromDeletions( scsel1, newdels )


            if True:
                istage += 1
                basename=self.pathway + "sc" + "%02i"%(istage)
                self.prepare_softcore(sclams,save=False,basename=basename)
                o = self.scti
                c1 = s2l( o.tiparm, o.merged_scmask1 )
                c2 = s2l( o.tiparm, o.merged_scmask2 )
                todelete = c1 + c2
                merged_timask1 = ReindexFromDeletions( s2l( o.tiparm, o.merged_timask1 ), todelete )
                merged_timask2 = ReindexFromDeletions( s2l( o.tiparm, o.merged_timask2 ), todelete )
                p = Strip( o.tiparm, l2s(todelete) )
                noshakemask = l2s( merged_timask1 + merged_timask2 )
                dirname = os.path.join(basename,"template")
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
                SaveParm( p, os.path.join(dirname,basename+".parm7") )
                parmed.tools.writeCoordinates( p, os.path.join(dirname,"%s.rst7"%(basename)) ).execute()

                write_stdti_mdin_script(sclams,basename,compinfo,
                                        "",
                                        noshakemask,
                                        l2s(merged_timask1),l2s(merged_timask2),
                                        "","",
                                        is_gas=(o.org_parm.box is None))
                    
            if len(self.scsel2) > 0:
                
                p = CopyParm( self.org_parm )
                for ik,k in enumerate( sorted(self.tisel2) ):
                    p.atoms[k].charge = self.ti2_charges[ik]
                    p.parm_data["CHARGE"][k] = self.ti2_charges[ik]
                tisel1 = copy.deepcopy( self.tisel2 )
                scsel1 = copy.deepcopy( self.scsel2 )
                molsel2 = copy.deepcopy( self.molsel1 )
                tisel1 = ReindexFromDeletions( tisel1, molsel2 )
                scsel1 = ReindexFromDeletions( scsel1, molsel2 )
                p = Strip( p, self.molmask1 )
                fseries = find_atoms_to_delete( p, scsel1 )
                for iser in range(len(fseries)):
                    for pser in range(iser):
                        fseries[iser] = ReindexFromDeletions( fseries[iser], fseries[pser] )
                istage += len(fseries)+1
                for iser in range(len(fseries)):
                    istage -= 1
                    basename=self.pathway + "sc" + "%02i"%(istage)
                    print("scsel",[p.atoms[idx].name for idx in scsel1])
                    newdels = fseries[iser]
                    print([p.atoms[idx].name for idx in newdels])
                    newdelsel = ListToSelection(newdels)
                    newtisel = ListToSelection(tisel1)
                    newscsel = ListToSelection(scsel1)
                    dirname = os.path.join(basename,"template")
                    if not os.path.exists(dirname):
                        os.makedirs(dirname)
                    crgmask = l2s(scsel1)
                    noshakemask = newtisel
                    write_stdti_mdin_script(sclams,basename,compinfo,
                                            crgmask,
                                            noshakemask,
                                            "",newdelsel,
                                            "",newdelsel,
                                            is_gas=(o.org_parm.box is None))
                    SaveParm( p, os.path.join(dirname,basename+".parm7") )
                    parmed.tools.writeCoordinates( p, os.path.join(dirname,"%s.rst7"%(basename)) ).execute()
                    p = Strip( p, newdelsel )
                    tisel1  = ReindexFromDeletions( tisel1, newdels )
                    scsel1  = ReindexFromDeletions( scsel1, newdels )
                    

        if len(self.scsel2) > 0:
            basename=self.pathway + "req"
            print "prepare_recharge"
            self.prepare_recharge(lams,save=True,basename=basename)
            o = self.reqti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(lams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))



        
##########################################################
##########################################################
# ALTERNATE TI PATHWAY (LIKE PISCTI, BUT USING LINEAR
# TRANSFORMATIONS INSTEAD OF PARAMETER INTERPOLATION)
##########################################################
##########################################################
    
class altsctisetup(sctisetup):
    
    ###################################################################################
    ###################################################################################
    def __init__( self,               \
                  parmfile, rstfile,  \
                  molmask1, molmask2, \
                  timask1,  timask2,  \
                  scmask1,  scmask2,  \
                  deletedihedrals=True   ):
        super(altsctisetup,self).__init__(parmfile,rstfile,molmask1,molmask2,timask1,timask2,scmask1,scmask2)
        self.pathway = "alt"
        self.true_dihedral_deletion = deletedihedrals
        self.assign_unique_lj_types_to_timask1()
        self.assign_unique_bond_types_to_timask1()
        self.assign_unique_angle_types_to_timask1()
        self.assign_unique_dihed_types_to_timask1()


    def test(self):
        calc_stdti_dvdls(self)

    def assign_unique_lj_types_to_timask1( self ):
        from collections import defaultdict as ddict
        
        self.new_lj_type_groups = []
        
        nbgrps = ddict(list)
        for o1 in self.ti_o1_to_o2:
            o2  = self.ti_o1_to_o2[o1]
            lj1 = self.org_parm.atoms[o1].nb_idx
            lj2 = self.org_parm.atoms[o2].nb_idx
            if lj1 != lj2:
                nbgrps[(lj1,lj2)].append( o1 )
                
        for lj in nbgrps:
            grp_members = []
            for o in nbgrps[lj]:
                oatm = self.org_parm.atoms[ o ]
                grp_members.append( oatm.idx )
            sele = ListToSelection( grp_members )
            #print "New LJ type for timask1 common atoms:",sele
            parmed.tools.addLJType( self.org_parm, sele ).execute()
            self.new_lj_type_groups.append( sele )



    def assign_unique_bond_types_to_timask1(self):
        from collections import defaultdict as ddict
        o1terms = []
        o2terms = []
        for bond in self.org_parm.bonds:
            if bond.atom1.idx in self.ti_o1_to_o2 and \
               bond.atom2.idx in self.ti_o1_to_o2:
                o1terms.append( bond )
            elif bond.atom1.idx in self.ti_o2_to_o1 and \
                 bond.atom2.idx in self.ti_o2_to_o1:
                o2terms.append( bond )

        tgrps = ddict(list)
        for o1t in o1terms:
            str1 = AtmStr([o1t.atom1,
                           o1t.atom2])
            fwdkey_target = ( self.ti_o1_to_o2[ o1t.atom1.idx ],
                              self.ti_o1_to_o2[ o1t.atom2.idx ] )
            revkey_target = tuple(reversed( fwdkey_target ))
            found = False
            for o2t in o2terms:
                str2 = AtmStr([o2t.atom1,
                               o2t.atom2])
                key = ( o2t.atom1.idx,
                        o2t.atom2.idx )
                if key == fwdkey_target or key == revkey_target:
                    #print "bond equiv %s <-> %s"%(str1,str2)
                    found=True
                    break
            if not found:
                raise Exception("Could not find an equivalent bond %s in timask2"%(str1))
            elif o1t.type.idx != o2t.type.idx:
                tgrps[ (o1t.type.idx,o2t.type.idx) ].append( o1t )

        self.new_bond_type_groups = []
        
        for tgrp in tgrps:
            term = tgrps[tgrp]
            sele1 = [ x.atom1.idx for x in term ]
            sele2 = [ x.atom2.idx for x in term ]
            mask1 = "@%s"%( ",".join( [ str(x+1) for x in sele1 ] ) )
            mask2 = "@%s"%( ",".join( [ str(x+1) for x in sele2 ] ) )
            #print "New bond type for timask1. %s <> %s"%(mask1,mask2)
            AddNewBondType( self.org_parm, sele1, sele2 )
            self.new_bond_type_groups.append( (sele1,sele2) )



    def assign_unique_angle_types_to_timask1(self):
        from collections import defaultdict as ddict
        o1terms = []
        o2terms = []
        for angle in self.org_parm.angles:
            if angle.atom1.idx in self.ti_o1_to_o2 and \
               angle.atom2.idx in self.ti_o1_to_o2 and \
               angle.atom3.idx in self.ti_o1_to_o2:
                o1terms.append( angle )
            elif angle.atom1.idx in self.ti_o2_to_o1 and \
                 angle.atom2.idx in self.ti_o2_to_o1 and \
                 angle.atom3.idx in self.ti_o2_to_o1:
                o2terms.append( angle )

        tgrps = ddict(list)
        for o1t in o1terms:
            str1 = AtmStr([o1t.atom1,
                           o1t.atom2,
                           o1t.atom3])
            fwdkey_target = ( self.ti_o1_to_o2[ o1t.atom1.idx ],
                              self.ti_o1_to_o2[ o1t.atom2.idx ],
                              self.ti_o1_to_o2[ o1t.atom3.idx ] )
            revkey_target = tuple(reversed( fwdkey_target ))
            found = False
            for o2t in o2terms:
                str2 = AtmStr([o2t.atom1,
                               o2t.atom2,
                               o2t.atom3])
                key = ( o2t.atom1.idx,
                        o2t.atom2.idx,
                        o2t.atom3.idx )
                if key == fwdkey_target or key == revkey_target:
                    #print "angle equiv %s <-> %s"%(str1,str2)
                    found=True
                    break
            if not found:
                raise Exception("Could not find an equivalent angle %s in timask2"%(str1))
            elif o1t.type.idx != o2t.type.idx:
                tgrps[ (o1t.type.idx,o2t.type.idx) ].append( o1t )

        self.new_angle_type_groups = []
        
        for tgrp in tgrps:
            term = tgrps[tgrp]
            sele1 = [ x.atom1.idx for x in term ]
            sele2 = [ x.atom2.idx for x in term ]
            sele3 = [ x.atom3.idx for x in term ]
            mask1 = "@%s"%( ",".join( [ str(x+1) for x in sele1 ] ) )
            mask2 = "@%s"%( ",".join( [ str(x+1) for x in sele2 ] ) )
            mask3 = "@%s"%( ",".join( [ str(x+1) for x in sele3 ] ) )
            #print "New angl type for timask1. %s <> %s <> %s"%(mask1,mask2,mask3)
            AddNewAngleType( self.org_parm, sele1, sele2, sele3 )
            self.new_angle_type_groups.append( (sele1,sele2,sele3) )



    def assign_unique_dihed_types_to_timask1(self):
        from parmed.topologyobjects import BondType,AngleType,DihedralType
        from collections import defaultdict as ddict

        self.o1_diheds_to_remove = ddict(list)
        self.o2_diheds_to_insert = ddict(list)

        
        o1terms = []
        o2terms = []
        for angle in self.org_parm.dihedrals:
            if angle.atom1.idx in self.ti_o1_to_o2 and \
               angle.atom2.idx in self.ti_o1_to_o2 and \
               angle.atom3.idx in self.ti_o1_to_o2 and \
               angle.atom4.idx in self.ti_o1_to_o2:
                o1terms.append( angle )
            elif angle.atom1.idx in self.ti_o2_to_o1 and \
                 angle.atom2.idx in self.ti_o2_to_o1 and \
                 angle.atom3.idx in self.ti_o2_to_o1 and \
                 angle.atom4.idx in self.ti_o2_to_o1:
                o2terms.append( angle )

        for o1t in o1terms:
            o1key = [o1t.atom1.idx,o1t.atom2.idx,o1t.atom3.idx,o1t.atom4.idx]
            str1 = AtmStr([o1t.atom1,o1t.atom2,o1t.atom3,o1t.atom4])
            found = False
            for o2t in o2terms:
                o2key = [o2t.atom1.idx,o2t.atom2.idx,o2t.atom3.idx,o2t.atom4.idx]
                o2key = [ self.ti_o2_to_o1[x] for x in o2key ]
                o2rev = list(reversed(o2key))
                str2 = AtmStr([o2t.atom1,o2t.atom2,o2t.atom3,o2t.atom4])
                if (o1key == o2key or o1key == o2rev ) and o1t.type == o2t.type:
                    found=True
                    break
            if not found:
                # this dihedral should be removed during the transformation
                # print "dihe to be removed %s %s"%(str1,str(o1t.type))
                self.o1_diheds_to_remove[o1t.type].append( o1t )

        for o2t in o2terms:
            o2s = [o2t.atom1.idx,o2t.atom2.idx,o2t.atom3.idx,o2t.atom4.idx]
            o1s = [ self.ti_o2_to_o1[x] for x in o2s ]
            str2 = AtmStr([self.org_parm.atoms[x] for x in o2s])
            str1 = AtmStr([self.org_parm.atoms[x] for x in o1s])
            found = False
            for o1t in o1terms:
                o1key = [o1t.atom1.idx,o1t.atom2.idx,o1t.atom3.idx,o1t.atom4.idx]
                o1rev = list(reversed(o1key))
                if (o1key == o1s or o1rev == o1s) and o1t.type == o2t.type:
                    found=True
                    break
            if not found:
                # this dihedral should be added during the transformation
                # print "dihe to be added   %s <-- %s"%(str1,str2)
                self.o2_diheds_to_insert[o2t.type].append( o2t )
    
        retyped_removals = ddict(list)
        #retyped_removals = ddict( bonded_term_collector )

        for t in self.o1_diheds_to_remove:
            newtype = DihedralType(t.phi_k,t.per,t.phase,t.scee,t.scnb,self.org_parm.dihedral_types)
            self.org_parm.dihedral_types.append( newtype )
            # print "Created new dihe type %s"%(str(self.org_parm.dihedral_types[-1]))
            for o1 in self.o1_diheds_to_remove[t]:
                o1sel = [o1.atom1.idx,o1.atom2.idx,o1.atom3.idx,o1.atom4.idx]
                mstr = AtmStr([self.org_parm.atoms[x] for x in o1sel])
                # print "dihe type reassigned %s"%(mstr)
                o1.type = newtype
                #retyped_removals[self.org_parm.dihedral_types[-1]].org_parm_objs.append( o1 )
                retyped_removals[newtype].append( o1 )
        self.o1_diheds_to_remove = retyped_removals

        insertions = ddict( bonded_term_collector )
        for t in self.o2_diheds_to_insert:
            for o in self.o2_diheds_to_insert[t]:
                insertions[t].org_parm_objs.append( o )
        self.o2_diheds_to_insert = insertions

        removals = ddict( bonded_term_collector )
        for t in self.o1_diheds_to_remove:
            for o in self.o1_diheds_to_remove[t]:
                removals[t].org_parm_objs.append( o )
        self.o1_diheds_to_remove = removals



            
            
    def reassign_unique_types_to_timask1_copy( self, tiobj, imodregion=1 ):
        for grpmask in self.new_lj_type_groups:
            sele = GetModifiedMask( tiobj, grpmask, imodregion=imodregion )
            #print "New LJ   type for copied timask1. idx1: %s"%(sele)
            parmed.tools.addLJType( tiobj.tiparm, sele ).execute()
        for grpmasks in self.new_bond_type_groups:
            sele1 = Org2ModifiedIndexs( tiobj, grpmasks[0], imodregion=imodregion )
            sele2 = Org2ModifiedIndexs( tiobj, grpmasks[1], imodregion=imodregion )
            mask1 = "@%s"%( ",".join( [ str(x+1) for x in sele1 ] ) )
            mask2 = "@%s"%( ",".join( [ str(x+1) for x in sele2 ] ) )
            #print "New bond type for copied timask1. %s <> %s"%(mask1,mask2)
            AddNewBondType(tiobj.tiparm, sele1, sele2)
        for grpmasks in self.new_angle_type_groups:
            sele1 = Org2ModifiedIndexs( tiobj, grpmasks[0], imodregion=imodregion )
            sele2 = Org2ModifiedIndexs( tiobj, grpmasks[1], imodregion=imodregion )
            sele3 = Org2ModifiedIndexs( tiobj, grpmasks[2], imodregion=imodregion )
            mask1 = "@%s"%( ",".join( [ str(x+1) for x in sele1 ] ) )
            mask2 = "@%s"%( ",".join( [ str(x+1) for x in sele2 ] ) )
            mask3 = "@%s"%( ",".join( [ str(x+1) for x in sele3 ] ) )
            #print "New angl type for copied timask1. %s <> %s <> %s"%(mask1,mask2,mask3)
            AddNewAngleType(tiobj.tiparm, sele1, sele2, sele3)
         

            
            
            

            
    def copy_params( self, tiobj, imodregion=1, lam=None ):
        self.copy_lj_params( tiobj, imodregion, lam=lam )
        self.copy_bond_params( tiobj, imodregion, lam=lam )
        self.copy_angle_params( tiobj, imodregion, lam=lam )

                
    def copy_lj_params( self, tiobj, imodregion=1, lam=None ):
        for o1 in self.ti_o1_to_o2:
            o2  = self.ti_o1_to_o2[o1]
            lj1 = tiobj.org_parm.atoms[o1].nb_idx
            lj2 = tiobj.org_parm.atoms[o2].nb_idx
            if lj1 == lj2:
                continue
            oatm = tiobj.org_parm.atoms[ o1 ]
            
            rmin0 = tiobj.org_parm.LJ_radius[oatm.nb_idx-1]
            eps0  = tiobj.org_parm.LJ_depth[oatm.nb_idx-1]
            
            midx = None
            if isinstance(oatm.modified_index,list):
                midx = oatm.modified_index[ imodregion ]
            else:
                midx = oatm.modified_index
            matm = tiobj.tiparm.atoms[ midx ] # the copied atom
            oatm = tiobj.org_parm.atoms[ o2 ]
            rmin = tiobj.org_parm.LJ_radius[oatm.nb_idx-1]
            eps  = tiobj.org_parm.LJ_depth[oatm.nb_idx-1]
            if lam is not None:
                rmin = rmin0 + lam * (rmin - rmin0)
                eps  = eps0  + lam * (eps  - eps0)
            try:
                depth14  = self.LJ_14_depth[atom.nb_idx-1]
                radius14 = self.LJ_14_radius[atom.nb_idx-1]
            except:
                depth14 = radius14 = None
            #print "%4s(%5i):%-4s LJ change from (%7.4f,%7.4f)"\
            #    %(matm.residue.name,matm.residue.idx+1,
            #      matm.name,matm.rmin,matm.epsilon),
            tiobj.tiparm.LJ_radius[matm.nb_idx-1] = rmin
            tiobj.tiparm.LJ_depth[matm.nb_idx-1]  = eps
            matm.atom_type.set_lj_params(eps,rmin, depth14, radius14)
            #print " to (%7.4f,%7.4f)"%(matm.rmin,matm.epsilon)
        tiobj.tiparm.recalculate_LJ()

      
    def copy_bond_params( self, tiobj, imodregion=1, lam=None ):
        import copy
        from parmed.topologyobjects import BondType,AngleType,DihedralType

        for grpmasks in self.new_bond_type_groups:
            o1sel1  = grpmasks[0]
            msel1   = Org2ModifiedIndexs( tiobj, o1sel1, imodregion=imodregion )
            o2sel1  = [ self.ti_o1_to_o2[x] for x in o1sel1 ]
            
            o1sel2  = grpmasks[1]
            msel2   = Org2ModifiedIndexs( tiobj, o1sel2, imodregion=imodregion )
            o2sel2  = [ self.ti_o1_to_o2[x] for x in o1sel2 ]

            # we want to copy the parameters from the bonds in
            # (o2sel1-o2sel2) to the bonds in (msel1-msel2)
            
            mbond = list( set( tiobj.tiparm.atoms[ msel1[0] ].bonds ) & set( tiobj.tiparm.atoms[ msel2[0] ].bonds ) )
            if len(mbond) != 1:
                raise Exception("Could not find bond in copy_bond_parms %i and %i because len is %i"%(msel1[0],msel2[0],len(mbond)))
            mbond = mbond[0]

            obond = list( set( tiobj.org_parm.atoms[ o2sel1[0] ].bonds ) & set( tiobj.org_parm.atoms[ o2sel2[0] ].bonds ) )
            if len(obond) != 1:
                raise Exception("Could not find bond in copy_bond_parms %i and %i because len is %i"%(o2sel1,o2sel2,len(obond)))
            obond = obond[0]

            zbond = list( set( tiobj.org_parm.atoms[ o1sel1[0] ].bonds ) & set( tiobj.org_parm.atoms[ o1sel2[0] ].bonds ) )
            if len(zbond) != 1:
                raise Exception("Could not find bond in copy_bond_parms %i and %i because len is %i"%(o1sel1,o1sel2,len(obond)))
            zbond = zbond[0]

            
            
            req = obond.type.req
            k   = obond.type.k
            req0= zbond.type.req
            k0  = zbond.type.k
            
            if lam is not None:
                req = req0 + lam*(req-req0)
                k   = k0   + lam*(k-k0)
#                print "changing bond parm from %i %.4f %.4f to"%(obond.type.idx,req0,k0),
#                print "%i %.4f %.4f using lam: %s"%(mbond.type.idx,req,k,str(lam))
            else:
                pass
#                print "changing bond parm from %i %.4f %.4f to"%(obond.type.idx,req0,k0),
#                print "%i %.4f %.4f using lam: %s"%(mbond.type.idx,req,k,str(lam))

            mbond.type.req = req
            mbond.type.k   = k



    def copy_angle_params( self, tiobj, imodregion=1, lam=None ):
        import copy
        from parmed.topologyobjects import BondType,AngleType,DihedralType

        for grpmasks in self.new_angle_type_groups:
            o1sel1  = grpmasks[0]
            msel1   = Org2ModifiedIndexs( tiobj, o1sel1, imodregion=imodregion )
            o2sel1  = [ self.ti_o1_to_o2[x] for x in o1sel1 ]
            
            o1sel2  = grpmasks[1]
            msel2   = Org2ModifiedIndexs( tiobj, o1sel2, imodregion=imodregion )
            o2sel2  = [ self.ti_o1_to_o2[x] for x in o1sel2 ]

            o1sel3  = grpmasks[2]
            msel3   = Org2ModifiedIndexs( tiobj, o1sel3, imodregion=imodregion )
            o2sel3  = [ self.ti_o1_to_o2[x] for x in o1sel3 ]

            mangle = None
            for t in tiobj.tiparm.atoms[ msel1[0] ].angles:
                if (t.atom1.idx == msel1[0] and \
                    t.atom2.idx == msel2[0] and \
                    t.atom3.idx == msel3[0]) or \
                    (t.atom1.idx == msel3[0] and \
                    t.atom2.idx == msel2[0] and \
                    t.atom3.idx == msel1[0]):
                    mangle = t
                    break

            oangle = None
            for t in tiobj.org_parm.atoms[ o2sel1[0] ].angles:
                if (t.atom1.idx == o2sel1[0] and \
                    t.atom2.idx == o2sel2[0] and \
                    t.atom3.idx == o2sel3[0]) or \
                    (t.atom1.idx == o2sel3[0] and \
                    t.atom2.idx == o2sel2[0] and \
                    t.atom3.idx == o2sel1[0]):
                    oangle = t
                    break


            zangle = None
            for t in tiobj.org_parm.atoms[ o1sel1[0] ].angles:
                if (t.atom1.idx == o1sel1[0] and \
                    t.atom2.idx == o1sel2[0] and \
                    t.atom3.idx == o1sel3[0]) or \
                    (t.atom1.idx == o1sel3[0] and \
                    t.atom2.idx == o1sel2[0] and \
                    t.atom3.idx == o1sel1[0]):
                    zangle = t
                    break
                
            
            theteq = oangle.type.theteq
            k      = oangle.type.k
            theteq0= zangle.type.theteq
            k0     = zangle.type.k

            if lam is not None:
                theteq = theteq0 + lam*(theteq-theteq0)
                k      = k0      + lam*(k-k0)
#                print "changing angle parm from %i %.4f %.4f to"%(oangle.type.idx,theteq0,k0),
#                print "%i %.4f %.4f using lam: %s"%(mangle.type.idx,theteq,k,str(lam))

            else:
                pass
#                print "changing angle parm from %i %.4f %.4f to"%(oangle.type.idx,theteq0,k0),
#                print "%i %.4f %.4f using lam: %s"%(mangle.type.idx,theteq,k,str(lam))
            mangle.type.theteq = theteq
            mangle.type.k      = k


    def delete_dihe_params( self, tiobj, imodregion=1, lam=None, onlymarkidx=False ):
        from parmed.topologyobjects import BondType,AngleType,DihedralType
        
        for t in self.o1_diheds_to_remove:
            if lam is None:
                zero_phi_k = False
                for o1 in self.o1_diheds_to_remove[t].org_parm_objs:
                    o1sel = [o1.atom1.idx,o1.atom2.idx,o1.atom3.idx,o1.atom4.idx]
                    msel  = Org2ModifiedIndexs( tiobj, o1sel, imodregion=imodregion )
                    mstr = AtmStr([tiobj.tiparm.atoms[x] for x in msel])
                    idx = None
                    for d in tiobj.tiparm.dihedrals:
                        if d.same_atoms( msel ):
                            if d.type == t:
                                idx = d
                                break
                    if idx is not None:
                        #print "deleting %s %s"%(mstr,str(idx.type))
                        self.o1_diheds_to_remove[t].tiparm_type_idx = idx.type.idx
                        if not onlymarkidx:
                            if not self.true_dihedral_deletion:
                                # delay the zeroing until we finish looking at diheds
                                zero_phi_k = True
                            else:
                                idx.delete()
                                tiobj.tiparm.dihedrals.remove(idx)
                    else:
                        print "FAILED TO DELETE %s %s"%(mstr,str(d.type))
                        print "Looking for quartet %s"%(mstr)
                        for d in tiobj.tiparm.dihedrals:
                            #astr = AtmStr( [d.atom1,d.atom2,d.atom3,d.atom4] )
                            if d.same_atoms( msel ):
                                print "ONLY FOUND %s"%(str(d.type))
                if zero_phi_k:
                    idx.type.phi_k = 0.

            else:
                idx = self.o1_diheds_to_remove[t].tiparm_type_idx
                dtype = tiobj.tiparm.dihedral_types[idx]
                dtype.phi_k *= (1.-lam)



                    
                    
                    
    def insert_dihe_params( self, tiobj, imodregion=1, lam=None ):
        from parmed.topologyobjects import Dihedral,DihedralType
        
        for t in self.o2_diheds_to_insert:
            if lam is None:
                tiobj.tiparm.dihedral_types.append \
                    (DihedralType(t.phi_k,t.per,t.phase,t.scee,t.scnb,tiobj.tiparm.dihedral_types))
                self.o2_diheds_to_insert[t].tiparm_type_idx = tiobj.tiparm.dihedral_types[-1].idx
                for o2 in self.o2_diheds_to_insert[t].org_parm_objs:
                    o2sel = [o2.atom1.idx,o2.atom2.idx,o2.atom3.idx,o2.atom4.idx]
                    o1sel = [ self.ti_o2_to_o1[x] for x in o2sel ]
                    msel  = Org2ModifiedIndexs( tiobj, o1sel, imodregion=imodregion )
                    atm1 = tiobj.tiparm.atoms[ msel[0] ]
                    atm2 = tiobj.tiparm.atoms[ msel[1] ]
                    atm3 = tiobj.tiparm.atoms[ msel[2] ]
                    atm4 = tiobj.tiparm.atoms[ msel[3] ]
                    mstr = AtmStr([tiobj.tiparm.atoms[x] for x in msel])
                    #print "adding   %s %s"%(mstr,str(tiobj.tiparm.dihedral_types[-1]))
                    tiobj.tiparm.dihedrals.append(
                        Dihedral(atm1, atm2, atm3, atm4, improper=o2.improper,
                                 ignore_end=o2.ignore_end, type=tiobj.tiparm.dihedral_types[-1])
                    )
            else:
                idx = self.o2_diheds_to_insert[t].tiparm_type_idx
                dtype = tiobj.tiparm.dihedral_types[idx]
                dtype.phi_k *= lam


    def replace_dihe_params( self, tiobj, imodregion=1, lam=None ):
        self.delete_dihe_params( tiobj, imodregion, lam=lam )
        self.insert_dihe_params( tiobj, imodregion, lam=lam )

           
    def pure_softcore_backend(self):
        
        self.scti = stdtimerge( self.org_parm,
                                self.molmask1,self.molmask2,
                                self.scmask1,self.scmask2,
                                keep_bondparms_from_mol = 2 )

        self.scti.merged_timask1 = self.scti.remask( self.scmask1 )
        self.scti.merged_timask2 = self.scti.remask( self.scmask2 )
        self.scti.merged_scmask1 = self.scti.remask( self.scmask1 )
        self.scti.merged_scmask2 = self.scti.remask( self.scmask2 )

        sc_atoms = GetSelectedAtomIndices(
            self.scti.tiparm,
            "(%s)|(%s)"%( self.scti.merged_scmask1, self.scti.merged_scmask2 ) )
        
        for a in sc_atoms:
            self.scti.tiparm.atoms[a].charge = 0.
            self.scti.tiparm.parm_data["CHARGE"][a] = 0.
            

        
    def prepare_softcore(self,lams,save=True,basename="altsc"):
        import copy
        from collections import defaultdict as ddict

        if False:
            super(altsctisetup,self).prepare_softcore(lams,save=save,basename=basename)
            self.copy_params( self.scti )
            self.delete_dihe_params( self.scti, imodregion=0 )
            self.insert_dihe_params( self.scti, imodregion=0 )

            ####
            # o1ds = ddict( list )
            # for d in self.scti.org_parm.dihedrals:
            #     if d.atom1.idx in self.ti_o1_to_o2 and \
            #        d.atom2.idx in self.ti_o1_to_o2 and \
            #        d.atom3.idx in self.ti_o1_to_o2 and \
            #        d.atom4.idx in self.ti_o1_to_o2:
            #         key = (d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx)
            #         o1ds[key].append( d )

            # o2ds = ddict( list )
            # for d in self.scti.org_parm.dihedrals:
            #     if d.atom1.idx in self.ti_o2_to_o1 and \
            #        d.atom2.idx in self.ti_o2_to_o1 and \
            #        d.atom3.idx in self.ti_o2_to_o1 and \
            #        d.atom4.idx in self.ti_o2_to_o1:
            #         key = tuple( [ self.ti_o2_to_o1[x] for x in [d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx] ] )
            #         o2ds[key].append( d )
                    
                    
            # m1ds = ddict( list )
            # m1sel = Org2ModifiedIndexs( self.scti, [x for x in self.ti_o1_to_o2], imodregion=0 )
            # for d in self.scti.tiparm.dihedrals:
            #     if d.atom1.idx in m1sel and \
            #        d.atom2.idx in m1sel and \
            #        d.atom3.idx in m1sel and \
            #        d.atom4.idx in m1sel:
            #         key = (d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx)
            #         m1ds[key].append( d )
                    
            # for o1key in o1ds:
            #     key = tuple( Org2ModifiedIndexs( self.scti, list(o1key), imodregion=0 ) )
            #     if key not in m1ds:
            #         key = tuple(list(reversed(list(key))))
            #         if key not in m1ds:
            #             print "Could not find o1ds key %s in m1ds"%(str(o1key))
            #             continue
            #     for t in o1ds[o1key]:
            #         found=False
            #         same_type=False
            #         same_parm=False
            #         for s in m1ds[key]:
            #             if t.type.idx == s.type.idx:
            #                 found=True
            #                 same_type=True
            #                 same_parm=True
            #                 break
            #             elif t.type.phi_k == s.type.phi_k and \
            #                  t.type.per == s.type.per and \
            #                  t.type.phase == s.type.phase and \
            #                  t.type.scee == s.type.scee and \
            #                  t.type.scnb == s.type.scnb:
            #                 found=True
            #                 same_parm=True
            #                 break
            #         print "o1d in m1d %8s   %15s   same_type=%8s same_parm=%8s"%(found,str(key),same_type,same_parm)
            # ####
            # for o1key in o2ds:
            #     key = tuple( Org2ModifiedIndexs( self.scti, list(o1key), imodregion=0 ) )
            #     if key not in m1ds:
            #         key = tuple(list(reversed(list(key))))
            #         if key not in m1ds:
            #             print "Could not find o2ds key %s in m1ds"%(str(o1key))
            #             continue
            #     for t in o2ds[o1key]:
            #         found=False
            #         same_type=False
            #         same_parm=False
            #         for s in m1ds[key]:
            #             if t.type.idx == s.type.idx:
            #                 found=True
            #                 same_type=True
            #                 same_parm=True
            #                 break
            #             elif t.type.phi_k == s.type.phi_k and \
            #                  t.type.per == s.type.per and \
            #                  t.type.phase == s.type.phase and \
            #                  t.type.scee == s.type.scee and \
            #                  t.type.scnb == s.type.scnb:
            #                 found=True
            #                 same_parm=True
            #                 break
            #         print "o2d in m1d %8s   %15s   same_type=%8s same_parm=%8s"%(found,str(key),same_type,same_parm)
            # ####

            
            
            #self.scti.tiparm.dihedrals.changed = True
            #self.scti.tiparm.update_dihedral_exclusions()

            for oatm_index in self.ti_o1_to_o2:
                # the original atom
                oatm = self.scti.org_parm.atoms[ oatm_index ]
                if oatm_index in self.ti_o1_to_o2:
                    catm_index = self.ti_o1_to_o2[oatm_index]
                elif oatm_index in self.ti_o2_to_o1:
                    catm_index = self.ti_o2_to_o1[oatm_index]
                else:
                    continue
                # the atom common to this original atom
                catm = self.scti.org_parm.atoms[ catm_index ]

                # if there is no common atom in the SC step, then we don't get this far
            
                if oatm.modified_index is not None:
                    # this atom was not deleted from this ti transformation
                    matm = self.scti.tiparm.atoms[ oatm.modified_index ] 
                    print "CopyQ from %5i:%4s:%-4s (%7.4f) to %5i:%4s:%-4s (%7.4f)"% \
                        (catm.idx,catm.residue.name,catm.name,catm.charge,\
                         matm.idx,matm.residue.name,matm.name,matm.charge)
                    matm.charge = catm.charge
                    self.scti.tiparm.parm_data["CHARGE"][matm.idx] = catm.charge
        else:
            self.pure_softcore_backend()

        ti = self.scti
        if True:
            for ires in range(len(ti.tiparm.residues)):
                 res = ti.tiparm.residues[ires]
                 if len(res.atoms) == 1:
                    atm = res.atoms[0]
                    if atm.atomic_number == 1:
                       if len(atm.bond_partners) > 0:
                          jres = atm.bond_partners[0].residue.idx
                       else:
                          raise Exception("PMEMD does not support residues consisting solely of an H atom")
                       pres = ti.tiparm.residues[jres]
                       if jres < ires:
                          print "merging :%i into :%i"%(res.idx+1,pres.idx+1)
                          pres.atoms.append( atm )
                          for a in pres.atoms:
                                a.residue = pres
                          res.atoms = []
                       else:
                          print "merging :%i into :%i"%(pres.idx+1,res.idx+1)
                          res.atoms.extend( pres.atoms )
                          for a in res.atoms:
                                a.residue = res
                          pres.atoms = []


        self.scti = ti
        self.scti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.scti.save(basename, subdir=dirname)
                

            
    def prepare_decharge(self,lams,save=True,basename="altdeq"):
        from collections import defaultdict as ddict

        super(altsctisetup,self).prepare_decharge(lams,save=False,basename=basename)
        self.reassign_unique_types_to_timask1_copy( self.deqti )
        self.copy_params( self.deqti )
        self.delete_dihe_params( self.deqti, imodregion=1 )
        self.insert_dihe_params( self.deqti, imodregion=1 )
        
                
        for oatm_index in range(len(self.deqti.org_parm.atoms)):
            # the original atom
            oatm = self.deqti.org_parm.atoms[ oatm_index ]
            if oatm_index in self.ti_o1_to_o2:
                catm_index = self.ti_o1_to_o2[oatm_index]
            elif oatm_index in self.ti_o2_to_o1:
                catm_index = self.ti_o2_to_o1[oatm_index]
            else:
                continue
            # the atom common to this original atom
            catm = self.deqti.org_parm.atoms[ catm_index ]

            # if there is no common atom in the SC step, then we don't get this far
            
            if oatm.modified_index is not None:
                # this atom was not deleted from this ti transformation
                if len(oatm.modified_index) > 1:
                    # this atom was copied in this ti transformation
                    matm = self.deqti.tiparm.atoms[ oatm.modified_index[1] ] # the copied atom
                    # print "CopyQ from %5i:%4s:%-4s (%7.4f) to %5i:%4s:%-4s (%7.4f)"%(\
                    #     catm.idx,catm.residue.name,catm.name,catm.charge,\
                    #     matm.idx,matm.residue.name,matm.name,matm.charge)
                    matm.charge = catm.charge
                    self.deqti.tiparm.parm_data["CHARGE"][matm.idx] = catm.charge

        ti = self.deqti
        if True:
            for ires in range(len(ti.tiparm.residues)):
                 res = ti.tiparm.residues[ires]
                 if len(res.atoms) == 1:
                    atm = res.atoms[0]
                    if atm.atomic_number == 1:
                       if len(atm.bond_partners) > 0:
                          jres = atm.bond_partners[0].residue.idx
                       else:
                          raise Exception("PMEMD does not support residues consisting solely of an H atom")
                       pres = ti.tiparm.residues[jres]
                       if jres < ires:
                          print "merging :%i into :%i"%(res.idx+1,pres.idx+1)
                          pres.atoms.append( atm )
                          for a in pres.atoms:
                                a.residue = pres
                          res.atoms = []
                       else:
                          print "merging :%i into :%i"%(pres.idx+1,res.idx+1)
                          res.atoms.extend( pres.atoms )
                          for a in res.atoms:
                                a.residue = res
                          pres.atoms = []


        self.deqti = ti                    
        self.deqti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.deqti.save(basename, subdir=dirname)


            
    def prepare_recharge(self,lams,save=True,basename="altreq"):
        import numpy as np
        if False:
            super(altsctisetup,self).prepare_recharge(lams,save=False,basename=basename)
        else:
            natom = len(self.org_parm.atoms)
            self.req_parm = CopyParm( self.org_parm )
            if self.req_parm.coordinates is None:
                raise Exception('Failure copying coordinates to self.req_parm')

            premol_init_common_idx  = 0
            premol_last_common_idx  = min( self.molsel1 + self.molsel2 ) - 1
            postmol_init_common_idx = max( self.molsel1 + self.molsel2 ) + 1
            postmol_last_common_idx = natom - 1

            mol2_natom = len(self.molsel2)

            for i,a in enumerate(self.req_parm.atoms):
                a.master_index = i
                
            mol = Extract( self.req_parm, self.molmask2 )
            postmask = "@%i-%i"%( postmol_init_common_idx+1, postmol_last_common_idx+1 )
            postmol  = Extract( self.req_parm, postmask )
            if premol_last_common_idx >= 0:
                premask = "@%i-%i"%( premol_init_common_idx+1, premol_last_common_idx+1 )
                premol  = Extract( self.req_parm, premask )
                master_array = GetSelectedAtomIndices( self.req_parm,premask ) \
                               + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                               + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                               + GetSelectedAtomIndices( self.req_parm,postmask )
                master_array = [ self.req_parm.atoms[idx].master_index for idx in master_array ]
            
                self.req_parm = Join(
                    premol, Join( mol, Join( mol, postmol ) ) )
            else:
                master_array = GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                               + GetSelectedAtomIndices( self.req_parm,self.molmask2 ) \
                               + GetSelectedAtomIndices( self.req_parm,postmask )
                master_array = [ self.req_parm.atoms[idx].master_index for idx in master_array ]
                self.req_parm = Join( mol, Join( mol, postmol ) )

            for idx,a in zip( master_array, self.req_parm.atoms):
                a.master_index = idx

            old_mol2_atom1 = min(self.molsel2)
            new_mol1_atom1 = premol_last_common_idx + 1
            new_mol2_atom1 = new_mol1_atom1 + mol2_natom
            new_molmask1   = "@%i-%i"%( new_mol1_atom1+1, new_mol1_atom1+mol2_natom )
            new_molmask2   = "@%i-%i"%( new_mol2_atom1+1, new_mol2_atom1+mol2_natom )
            new_tisel1     = [ i-old_mol2_atom1+new_mol1_atom1 for i in self.tisel2 ]
            new_tisel2     = [ i-old_mol2_atom1+new_mol2_atom1 for i in self.tisel2 ]

            new_timask1    = ListToSelection( new_tisel1 )
            new_timask2    = ListToSelection( new_tisel2 )
            new_scsel1     = [ i-old_mol2_atom1+new_mol1_atom1 for i in self.scsel2 ]
            new_scsel2     = [ i-old_mol2_atom1+new_mol2_atom1 for i in self.scsel2 ]
            new_scmask1    = ListToSelection( new_scsel1 )
            new_scmask2    = ListToSelection( new_scsel2 )

            ti = stdtimerge( self.req_parm,
                             new_molmask1, new_molmask2,
                             new_scmask1,  new_scmask2,
                             keep_bondparms_from_mol = 2 )

            ti.merged_timask1 = ti.remask( new_scmask1 )
            ti.merged_timask2 = ti.remask( new_scmask2 )
            ti.merged_scmask1 = ""
            ti.merged_scmask2 = ""
        
            sc_atoms = GetSelectedAtomIndices( ti.tiparm, ti.remask( new_scmask1 ) )
            for a in sc_atoms:
                ti.tiparm.atoms[a].charge = 0.
                ti.tiparm.parm_data["CHARGE"][a] = 0.

            tiorg = CopyParm( self.org_parm )
            for iat,a in enumerate(tiorg.atoms):
                a.master_index = iat
                a.original_index = [iat]
                a.modified_index = []
            for iat,a in enumerate(ti.tiparm.atoms):
                for i in range(len(a.original_index)):
                    o = a.original_index[i]
                    m = self.req_parm.atoms[o].master_index
                    a.original_index[i] = m
                    tiorg.atoms[m].modified_index.append(iat)
            for a in tiorg.atoms:
                if len(a.modified_index) == 0:
                    a.modified_index = None
                elif len(a.modified_index) == 1:
                    a.modified_index = a.modified_index[0]
            ti.org_parm = tiorg
        
#            self.reqti = ti
#            self.reqti.tiparm.remake_parm()

            for ires in range(len(ti.tiparm.residues)):
                 res = ti.tiparm.residues[ires]
                 if len(res.atoms) == 1:
                    atm = res.atoms[0]
                    if atm.atomic_number == 1:
                       if len(atm.bond_partners) > 0:
                          jres = atm.bond_partners[0].residue.idx
                       else:
                          raise Exception("PMEMD does not support residues consisting solely of an H atom")
                       pres = ti.tiparm.residues[jres]
                       if jres < ires:
                          print "merging :%i into :%i"%(res.idx+1,pres.idx+1)
                          pres.atoms.append( atm )
                          for a in pres.atoms:
                                a.residue = pres
                          res.atoms = []
                       else:
                          print "merging :%i into :%i"%(pres.idx+1,res.idx+1)
                          res.atoms.extend( pres.atoms )
                          for a in res.atoms:
                                a.residue = res
                          pres.atoms = []

            self.reqti = ti
            self.reqti.tiparm.remake_parm()

        
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.reqti.save(basename, subdir=dirname)
                
        

##########################################################
##########################################################
# PARAMETER INTERPOLATED TI PATHWAY
##########################################################
##########################################################
    
        
class pisctisetup(altsctisetup):
    
    ###################################################################################
    ###################################################################################
    def __init__( self,               \
                  parmfile, rstfile,  \
                  molmask1, molmask2, \
                  timask1,  timask2,  \
                  scmask1,  scmask2   ):
        super(pisctisetup,self).__init__(parmfile,rstfile,molmask1,molmask2,timask1,timask2,scmask1,scmask2)
        self.pathway = "pi"
        self.true_dihedral_deletion = False # there's no other sensible option


    def test(self):
        calc_piti_dvdls(self)

    def scale_deq_charges( self, tiobj, lam ):

        for o1 in self.ti_o1_to_o2:
            # the original atom
            o2    = self.ti_o1_to_o2[o1]
            o1atm = tiobj.org_parm.atoms[o1]
            o2atm = tiobj.org_parm.atoms[o2]
            q0 = o1atm.charge
            q1 = o2atm.charge
            q  = q0 + lam*(q1-q0)
            matm = tiobj.tiparm.atoms[ o1atm.modified_index[1] ] # the copied atom
            matm.charge = q
            tiobj.tiparm.parm_data["CHARGE"][matm.idx] = q
        for o1 in self.scsel1:
            o1atm = tiobj.org_parm.atoms[o1]
            q0 = o1atm.charge
            q1 = 0.
            q  = q0 + lam*(q1-q0)
            matm = tiobj.tiparm.atoms[ o1atm.modified_index[1] ]
            matm.charge = q
            tiobj.tiparm.parm_data["CHARGE"][matm.idx] = q

            
        
    def prepare_decharge(self,lams,save=True,basename="pideq"):
        from collections import defaultdict as ddict
        import copy

        print "PISCTI DEQ SETUP"
        sctisetup.prepare_decharge(self,lams,save=False,basename=basename)
        self.deqti.basename=basename
        self.deq_parms = ddict(int)
        
        if True:
            self.reassign_unique_types_to_timask1_copy( self.deqti )
            self.insert_dihe_params( self.deqti, imodregion=1, lam=None )
            self.delete_dihe_params( self.deqti, imodregion=1, lam=None, onlymarkidx=True )
                    
            parm0 = None
            parm1 = None
            for lam in [0,1]:
                print "PISCTI DEQ ENDPOINT LAM %.8f"%(lam)
                tiobj = copy.deepcopy( self.deqti )
                tiobj.tiparm = CopyParm( self.deqti.tiparm )
                tiobj.org_parm = self.deqti.org_parm

                self.copy_params( tiobj, lam=lam )
                self.delete_dihe_params( tiobj, imodregion=1, lam=lam )
                self.insert_dihe_params( tiobj, imodregion=1, lam=lam )
                self.scale_deq_charges( tiobj, lam )
                tiobj.tiparm.remake_parm()
                tiobj.tiparm.strip( tiobj.merged_timask1 )
                #self.deq_parms[lam] = tiobj.tiparm
                if lam == 0:
                    parm0 = tiobj.tiparm
                else:
                    parm1 = tiobj.tiparm
            for lam in lams:
                print "PISCTI DEQ INTERPOLATING LAM %.8f"%(lam)
                if lam == 0.0:
                    self.deq_parms[lam] = parm0
                elif lam == 1.0:
                    self.deq_parms[lam] = parm1
                else:
                    self.deq_parms[lam] = nonsoftcore_interpolation(parm0,parm1,lam)

        if save:
            import os.path
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            for ilam,lam in enumerate(sorted(self.deq_parms)):
                parm = self.deq_parms[lam]
                rebase = os.path.join(dirname,"%s_%.8f"%(basename,lam))
                SaveParm( parm, "%s.parm7"%(rebase) )
                #parmed.tools.writeCoordinates( parm, "%s.rst7"%(rebase) ).execute()
                if ilam == 0:
                    parmed.tools.writeCoordinates( parm, os.path.join(dirname,"%s.rst7"%(basename))).execute()

            
    def prepare_softcore(self,lams,save=True,basename="pisc"):
        from collections import defaultdict as ddict
        print "PICSTI SC  SETUP"
        super(pisctisetup,self).prepare_softcore(lams,save=False,basename=basename)
        self.scti.basename=basename

        self.sc_parms = ddict(int)
        
        if True:
            for ilam,lam in enumerate(lams):
                print "PISCTI SC  INTERPOLATING LAM %.8f"%(lam)
                parm = softcore_interpolation(\
                        self.scti.tiparm, self.scti.merged_timask1, self.scti.merged_timask2, lam )
                self.sc_parms[lam] = parm
        if save:
            import os.path
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            for ilam,lam in enumerate(sorted(self.sc_parms)):
                parm = self.sc_parms[lam]
                rebase = os.path.join(dirname,"%s_%.8f"%(basename,lam))
                SaveParm( parm, "%s.parm7"%(rebase) )
                #parmed.tools.writeCoordinates( parm, "%s.rst7"%(rebase) ).execute()
                if ilam == 0:
                    parmed.tools.writeCoordinates( parm, os.path.join(dirname,"%s.rst7"%(basename))).execute()

                
    def prepare_recharge(self,lams,save=True,basename="pireq"):
        from collections import defaultdict as ddict
        print "PISCTI REQ SETUP"
        super(pisctisetup,self).prepare_recharge(lams,save=False,basename=basename)
        self.reqti.basename=basename

        self.req_parms = ddict(int)
        if False:
            v0ats = self.reqti.merged_timask1
            v1ats = self.reqti.merged_timask2
            req_endstates = []
            for lam in [0,1]:
                parm = CopyParm( self.req.tiparm )
                for v0a,v1a in zip(v0ats,v1ats):
                    q0 = parm.atoms[v0a].charge
                    q1 = parm.atoms[v1a].charge
                    parm.atoms[v1a].charge = q0 + lam*(q1-q0)
                parm.strip( self.reqti.merged_timask1 )
                parm.remake_parm()
                req_endstates.append( parm )
            for lam in lams:
                self.req_parms[lam] = nonsoftcore_interpolation(\
                        req_endstates[0],req_endstates[1],lam)
        elif False:
            v0ats = GetSelectedAtomIndices(self.reqti.tiparm, self.reqti.merged_timask1)
            v1ats = GetSelectedAtomIndices(self.reqti.tiparm, self.reqti.merged_timask2)
            for ilam,lam in enumerate(lams):
                parm = CopyParm( self.reqti.tiparm )
                for v0a,v1a in zip(v0ats,v1ats):
                    q0 = parm.atoms[v0a].charge
                    q1 = parm.atoms[v1a].charge
                    parm.atoms[v1a].charge = q0 + lam*(q1-q0)
                parm.strip( self.reqti.merged_timask1 )
                parm.remake_parm()
                self.req_parms[lam] = parm
        else:
            ats = GetSelectedAtomIndices(self.org_parm, self.scmask2)
            for lam in lams:
                print "PISCTI REQ INTERPOLATING LAM %.8f"%(lam)
                parm = CopyParm( self.org_parm )
                for a in ats:
                    parm.atoms[a].charge *= lam
                    #if lam == 0.0: print AtmStr( parm.atoms[a] )
                parm.strip( self.molmask1 )
                #parm.remake_parm()
                self.req_parms[lam] = parm
            
        if save:
            import os.path
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            for ilam,lam in enumerate(sorted(self.req_parms)):
                parm = self.req_parms[lam]
                rebase = os.path.join(dirname,"%s_%.8f"%(basename,lam))
                SaveParm( parm, "%s.parm7"%(rebase) )
                if ilam == 0:
                    parmed.tools.writeCoordinates( parm, os.path.join(dirname,"%s.rst7"%(basename))).execute()

                
    def execute(self,lams,compinfo,sclams=None):
        #from parmutils import ListToSelection as l2s
        #from parmutils import GetSelectedAtomIndices as s2l
        l2s=ListToSelection
        s2l=GetSelectedAtomIndices

        basename=self.pathway + "deq"
        self.prepare_decharge(lams,save=True,basename=basename)
        o = self.deqti
        crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
        noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
        write_stdti_mdin_script( lams,basename,compinfo,
                                 crgmask,
                                 noshakemask,
                                 o.merged_timask1,o.merged_timask2,
                                 o.merged_scmask1,o.merged_scmask2,
                                 is_piti=True,
                                 is_gas=(o.org_parm.box is None) )

        if len(self.scsel1 + self.scsel2) > 0:
            if sclams is None:
                sclams=lams
            basename=self.pathway + "sc"
            self.prepare_softcore(sclams,save=True,basename=basename)
            o = self.scti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script( sclams,basename,compinfo,
                                     crgmask,
                                     noshakemask,
                                     o.merged_timask1,o.merged_timask2,
                                     o.merged_scmask1,o.merged_scmask2,
                                     is_piti=True,
                                     is_gas=(o.org_parm.box is None) )

        if len(self.scsel2) > 0:
            basename=self.pathway + "req"
            self.prepare_recharge(lams,save=True,basename=basename)
            o = self.reqti
            crgmask = l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script( lams,basename,compinfo,
                                     crgmask,
                                     noshakemask,
                                     o.merged_timask1,o.merged_timask2,
                                     o.merged_scmask1,o.merged_scmask2,
                                     is_piti=True,
                                     is_gas=(o.org_parm.box is None) )
                    









##########################################################
##########################################################
# STANDARD TI PATHWAY
##########################################################
##########################################################
    
class unitisetup(sctisetup):
    
    ###################################################################################
    ###################################################################################
    def __init__( self,               \
                  parmfile, rstfile,  \
                  molmask1, molmask2, \
                  timask1,  timask2,  \
                  scmask1,  scmask2 ):
        #print "scmask1=",scmask1
        self.refit="" #refit
        ti_init(self,parmfile,rstfile,molmask1,molmask2,timask1,timask2,scmask1,scmask2)
        self.true_dihderal_deletion = False # not possible otherwise
        self.prepare_common_atom_maps()
        self.pathway = "uni"

    def execute(self,lams,compinfo,sclams=None):
        #from parmutils import ListToSelection as l2s
        #from parmutils import GetSelectedAtomIndices as s2l
        l2s=ListToSelection
        s2l=GetSelectedAtomIndices

        #if len(self.scsel1 + self.scsel2) > 0:
        if True:
            if sclams is None:
                sclams = lams
            basename=self.pathway + "sc"
            print "prepare_softcore"
            self.prepare_softcore(sclams,save=True,basename=basename)
            o = self.scti
            c1 = s2l( o.tiparm, o.merged_scmask1 )
            c2 = s2l( o.tiparm, o.merged_scmask2 )
            crgmask = "" # l2s( s2l( o.tiparm, o.merged_scmask1 ) + s2l( o.tiparm, o.merged_scmask2 ) )
            noshakemask = l2s( s2l( o.tiparm, o.merged_timask1 ) + s2l( o.tiparm, o.merged_timask2 ) )
            write_stdti_mdin_script(sclams,basename,compinfo,
                                    crgmask,
                                    noshakemask,
                                    o.merged_timask1,o.merged_timask2,
                                    o.merged_scmask1,o.merged_scmask2,
                                    is_gas=(o.org_parm.box is None))

        
                
    ###################################################################################
    ###################################################################################
    def prepare_softcore(self,lams,save=True,basename="unisc"):

        if len(self.molsel2) > 0:
            self.scti = stdtimerge( self.org_parm,
                                    self.molmask1,self.molmask2,
                                    self.timask1,self.timask2,
                                    keep_bondparms_from_mol = 2 )

            self.scti.merged_timask1 = self.scti.remask( self.timask1 )
            self.scti.merged_timask2 = self.scti.remask( self.timask2 )
            self.scti.merged_scmask1 = self.scti.remask( self.scmask1 )
            self.scti.merged_scmask2 = self.scti.remask( self.scmask2 )
            
            
        else:
            self.scti = stdtimerge( self.org_parm,
                                    self.molmask1,self.molmask2,
                                    self.timask1,self.timask2,
                                    keep_bondparms_from_mol = 1 )
            
            self.scti.merged_timask1 = self.scti.remask( self.timask1 )
            self.scti.merged_timask2 = ""
            self.scti.merged_scmask1 = self.scti.remask( self.scmask1 )
            self.scti.merged_scmask2 = ""
            

            
        self.scti.tiparm.remake_parm()
        if save:
            import os
            dirname = os.path.join(basename,"template")
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            self.scti.save( basename, subdir=dirname )
        





            
##########################################################
##########################################################
# END
##########################################################
##########################################################



def maincli():
    import argparse

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""
      """,
      epilog="""
THE COPY UTILITY
   
   To setup a TI calculation you need appropriately crafted parm7 and rst7
   files.  In typical siutations, this means that you have two copies
   of your molecule. In a standard (non-TI) MD simulation, you'd have:

     (let's call this md.parm7 and md.rst7)
      molecule
      nonwater (ions)
      solvent  (waters)

   but the input of this script requires the system to look like:

     (let's call this tigen_input.parm7 tigen_input.rst7)
      molecule (lambda=0)
      molecule (lambda=1)
      nonwater (ions)
      solvent  (waters)

   Given md.parm7/rst7, tigen.py contains a utility to help you prepare
   tigen_input.parm7/rst7.  For example, suppose the molecule was ':1-10',
   then you could run:
 
   tigen.py -p md.parm7 -c md.rst7 --molmask=":1-10" --copy="tigen_input"

   This would create a series of files:
      tigen_input.sh
      tigen_input.lib
      tigen_input.frcmod
      tigen_input.mol0.pdb
      tigen_input.mol1.pdb
      tigen_input.nonwater.pdb
      tigen_input.solvent.pdb

   The tigen_input.sh script will use tleap to construct tigen_input.parm7 and
   tigen_input.rst7.  If you are doing a pKa calculation, then you should edit
   the sh script by adding tleap commands to change the the atomic charges in
   mol1; e.g., to change the N and H1 charges in residue 5 within mol1:
      set mol1.5.N    charge  0.159200
      set mol1.5.H1   charge  0.198400

   If you are replacing a residue with a different kind of residue, then you
   should delete or change the appropriate atoms in tigen_input.mol1.pdb to
   match the residue template found in your Amber OFF (lib) file.  This may
   require you to change the ordering of the atoms.  You should not manually
   ADD atoms to the pdb file, nor should you change the atom coordinates.
   You should only: delete atoms, rename atoms, and/or swap lines.

   The only time that you wouldn't need to prepare tigen_input.parm7/rst7 is
   for absolute solvation energies.
   

EXAMPLE TI TRANSFORMATIONS

   1. Absolute solvation energy of a small molecule
   tigen.py --molmask1=":1" --timask1=":1" --scmask1=":1"

   2. Change charge vector, e.g., pKa calculation
   tigen.py --molmask1=":1-10"  --timask1=":5" \\
            --molmask2=":11-20" --timask2=":15"

   3. Relative free energy from changing a functional group
   tigen.py --molmask1=":1" --timask1=":1" --scmask1=":1@H1" \\
            --molmask2=":2" --timask2=":2" --scmask2=":2@OH,HO"

   3. Deazo substitution (changing a C into a N-H -- we add a H)
   tigen.py --molmask1=":1-10"  --timask1=":5" \\
            --molmask2=":11-20" --timask2=":15" --scmask2=":15@HN"


CALCULATION OF THE NUMBER OF ALLOCATED NODES

The script will choose to use the fewest number of nodes with
the maximum occupancy that is a multiple of the number of lambda
simulations.  The choice is influenced by the options --min-nodes
--max-nodes --cpus-per-node and --nlambda

   1. --nlambda=12 --min-nodes=12 --max-nodes=12 --cpus-per-node=24
        ...will use 12 nodes and 24 cores/node (occupancy=24/24)
        12 nodes*24 cores/node = 288 is a multiple of nlambda

   2. --nlambda=12 --min-nodes=1  --max-nodes=12 --cpus-per-node=24
        ...will use  1 node and 24 cores/node (occupancy=24/24)
        1 nodes*24 cores/node = 24 is a multiple of nlambda

   3. --nlambda=12 --min-nodes=1  --max-nodes=12 --cpus-per-node=28
        ...will use 3 nodes; 28 cores/node (occupancy=(28/28)
        1 nodes*28 cores/node = 28 is NOT a multiple of nlambda
        2 nodes*28 cores/node = 56 is NOT a multiple of nlambda
        3 nodes*28 cores/node = 84 is a multiple of nlambda

   4. --nlambda=12 --min-nodes=1  --max-nodes=2  --cpus-per-node=28
        ...will use 1 node;  24 cores/node (occupancy=(24/28)
        1 nodes*28 cores/node = 28 is NOT a multiple of nlambda
        2 nodes*28 cores/node = 56 is NOT a multiple of nlambda
        1 nodes*27 cores/node = 27 is NOT a multiple of nlambda
        2 nodes*27 cores/node = 54 is NOT a multiple of nlambda
        1 nodes*26 cores/node = 26 is NOT a multiple of nlambda
        2 nodes*26 cores/node = 52 is NOT a multiple of nlambda
        1 nodes*25 cores/node = 25 is NOT a multiple of nlambda
        2 nodes*25 cores/node = 50 is NOT a multiple of nlambda
        1 nodes*24 cores/node = 24 is a multiple of nlambda

As you can see, the decision is made by looping from min-nodes to max-nodes and
asking whether the total number of cores is a multiple of nlambda.  If max-nodes
is reached, then the number of cpus-per-node is decreased by 1 and the loop is
performed again. The default behavior is to allocate nlambda nodes by setting:
   max-nodes = min-nodes = nlambda
If you want to force the use of X nodes, then use the options: 
  --max-nodes=X --min-nodes=X
""")
    


    parser.add_argument \
        ("-n","--nlambda",
         help="number of states [default: 12]",
         type=int,
         default=12,
         required=False )
    
    parser.add_argument \
        ("--nlambda-softcore",
         help="number of states for the softcore stage. If -1, then same as --nlambda [default: -1]",
         type=int,
         default=-1,
         required=False )
    
    parser.add_argument \
        ("-p","--parm",
         help="amber parm7 file",
         type=str,
         required=True )
    
    parser.add_argument \
        ("-c","--crd",
         help="amber restart file",
         type=str,
         required=True )

    parser.add_argument \
        ("--molmask1",
         help="amber mask selecting the first copy of the molecule",
         type=str,
         required=True )

    parser.add_argument \
        ("--molmask2",
         help="amber mask selecting the second copy of the molecule [default: '']",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--timask1",
         help="amber mask selecting the TI-region in the first molecule",
         type=str,
         required=False )

    parser.add_argument \
        ("--timask2",
         help="amber mask selecting the TI-region in the second molecule [default: '']",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--scmask1",
         help="amber mask selecting the softcore atoms in the first TI region [default: '']",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--scmask2",
         help="amber mask selecting the softcore atoms in the second TI region [default: '']",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--stdti",
         help="prepare a standard TI pathway calculations",
         action='store_true' )

    parser.add_argument \
        ("--altti",
         help="prepare the alternate TI pathway calculations",
         action='store_true' )

    parser.add_argument \
        ("--piti",
         help="prepare the parameter-interpolated TI pathway calculations",
         action='store_true' )


    parser.add_argument \
        ("--serti",
         help="prepare the series TI pathway calculations",
         action='store_true' )


    parser.add_argument \
        ("--uniti",
         help="prepare the unified protocol TI calculations",
         action='store_true' )
    
    parser.add_argument \
        ("--test",
         help="calculate DV/DL as a function of lambda using the fixed input crds",
         action='store_true' )

    parser.add_argument \
        ("--refit",
         help="refit the non-softcore TI atom charges within the intermediate end-states",
         action='store_true' )

    parser.add_argument \
        ("--refitq",
         help="refit the non-softcore TI atom charges within the intermediate end-states under the constraint that commonly charged atoms retain their commonality in the fit",
         action='store_true' )

    parser.add_argument \
        ("--refith",
         help="refit the non-softcore TI atom charges within the intermediate end-states under the constraint that commonly charged hydrogens retain their commonality",
         action='store_true' )

    
    
    parser.add_argument \
        ("--copy",
         help="write pdb, sh, lib, and frcmod files for use in tleap",
         type=str,
         default="",
         required=False )
    
    parser.add_argument \
        ("--nstlim",
         help="number of MD steps per restart [default: 1000000]",
         type=int,
         default=1000000,
         required=False )

    parser.add_argument \
        ("--ntpr",
         help="steps per print [default: 20000]",
         type=int,
         default=20000,
         required=False )

    parser.add_argument \
        ("--ntwx",
         help="steps per trajectory write [default: 20000]",
         type=int,
         default=20000,
         required=False )

    parser.add_argument \
        ("--numexchg",
         help="number of replica exchange attempts per simulation [default: 0]",
         type=int,
         default=0,
         required=False )



    parser.add_argument \
        ("--gti",
         help="if present, then set the gti_add_sc=1, gti_cut=1 flag in TEMPLATE.mdin",
         action='store_true',
         default=False,
         required=False )


    parser.add_argument \
        ("--nmropt",
         help="if 1, then add a restraint section to TEMPLATE.mdin [default: 0]",
         type=int,
         default=0,
         required=False )

    
    parser.add_argument \
        ("--cut",
         help="Nonbond cutoff [default: 10.0]",
         type=float,
         default=10.,
         required=False )


    parser.add_argument \
        ("--parallel-cpu",
         help="Run all lambdas at the same time on CPUs using a groupfile",
         action="store_true",
         default=False,
         required=False )

    parser.add_argument \
        ("--serial-cpu",
         help="Run lambdas sequentially on one-or-more CPU nodes",
         action="store_true",
         default=False,
         required=False )

    parser.add_argument \
        ("--parallel-gpu",
         help="Run all lambdas at the same time on GPUs using a groupfile",
         action="store_true",
         default=False,
         required=False )

    parser.add_argument \
        ("--serial-gpu",
         help="Run lambdas sequentially on a single GPU",
         action="store_true",
         default=False,
         required=False )


    parser.add_argument \
        ("--oversubscribe",
         help="Run all lambdas at the same time on all GPUs from 1 node",
         action="store_true",
         default=False,
         required=False )

    
    parser.add_argument \
        ("--cpu-partition",
         help="name of the CPU slurm partition [default: main]",
         type=str,
         default="main",
         required=False )

    parser.add_argument \
        ("--gpu-partition",
         help="name of the GPU slurm partition default: gpu]",
         type=str,
         default="gpu",
         required=False )

    parser.add_argument \
        ("--cpus-per-node",
         help="number of CPU cores per node [default: 24]",
         type=int,
         default=24,
         required=False )

    parser.add_argument \
        ("--gpus-per-node",
         help="number of GPUs per node [default: 4]",
         type=int,
         default=4,
         required=False )


    parser.add_argument \
        ("--qos",
         help="SBATCH qos option, necessary for some clusters",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--account",
         help="SBATCH account option, necessary for some clusters",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--exclude",
         help="SBATCH node exclusion option",
         type=str,
         default="",
         required=False )

    parser.add_argument \
        ("--constraint",
         help="SBATCH node constraint option",
         type=str,
         default="",
         required=False )

    
    parser.add_argument \
        ("--min-nodes",
         help="minimum number of CPU nodes to use (if less than 1, then it is set to nlambda) [default: -1]",
         type=int,
         default=-1,
         required=False )

    parser.add_argument \
        ("--max-nodes",
         help="maximum number of CPU nodes to use (if less than 1, then it is set to nlambda) [default: -1]",
         type=int,
         default=-1,
         required=False )

    parser.add_argument \
        ("--days",
         help="time limit in days [default: 2]",
         type=int,
         default=2,
         required=False )
    
    parser.add_argument \
        ("--launch",
         help="launch command [default: 'srun --mpi=pmi2']",
         type=str,
         default="srun --mpi=pmi2",
         required=False )
    

    parser.add_argument \
        ("--xstream",
         help="overrides gpu options for xstream.stanford.xsede.org and uses a special slurm header for that particular system",
         action='store_true' )


    parser.add_argument \
        ("--steps-per-slurm",
         help="production cycles per slurm file [default: 10]",
         type=int,
         default=10,
         required=False )

    parser.add_argument \
        ("--steps",
         help="total number of production cycles [default: 100]",
         type=int,
         default=100,
         required=False )

    
    

    args = parser.parse_args()

    if args.copy != "":
        if args.stdti:
            raise Exception("can't use --copy with --stdti")
        if args.serti:
            raise Exception("can't use --copy with --serti")
        elif args.altti:
            raise Exception("can't use --copy with --altti")
        elif args.piti:
            raise Exception("can't use --copy with --piti")
        elif args.uniti:
            raise Exception("can't use --copy with --uniti")
        elif len(args.molmask2) > 0:
            raise Exception("can't use --copy with --molmask2")
        elif args.timask1 is not None:
            raise Exception("can't use --copy with --timask1")
        elif len(args.timask2) > 0:
            raise Exception("can't use --copy with --timask2")
        elif len(args.scmask1) > 0:
            raise Exception("can't use --copy with --scmask1")
        elif len(args.scmask2) > 0:
            raise Exception("can't use --copy with --scmask2")
        elif args.test:
            raise Exception("can't use --copy with --test")
        elif len(args.molmask1) == 0:
            raise Exception("--copy requires --molmask1")
        p = OpenParm( args.parm, xyz=args.crd )
        ticopy( p, args.molmask1, base=args.copy )
        exit(0)
    else:
        if args.timask1 is None:
            raise Exception("error: argument --timask1 is required")

    
    cinfo = computer_info()
    cinfo.num_cpu_per_node = args.cpus_per_node
    cinfo.num_gpu_per_node = args.gpus_per_node
    cinfo.min_cpu_nodes = args.min_nodes
    cinfo.max_cpu_nodes = args.max_nodes
    cinfo.num_days  =  args.days
    cinfo.nstlim = args.nstlim
    cinfo.ntpr = args.ntpr
    cinfo.ntwx = args.ntwx
    cinfo.numexchg = args.numexchg
    cinfo.xstream = args.xstream
    cinfo.cpu_partition = args.cpu_partition
    cinfo.gpu_partition = args.gpu_partition
    cinfo.steps_per_slurm = args.steps_per_slurm
    cinfo.steps = args.steps
    cinfo.launch = args.launch
    cinfo.exclude= args.exclude
    cinfo.constraint= args.constraint
    cinfo.nmropt = args.nmropt
    cinfo.gti = args.gti

    cinfo.cut = args.cut

    cinfo.oversub = args.oversubscribe
    cinfo.parallel_cpu=False
    cinfo.parallel_gpu=False
    cinfo.serial_cpu=False
    cinfo.serial_gpu=False
    if args.parallel_cpu:
        cinfo.parallel_cpu=True
    elif args.parallel_gpu:
        cinfo.parallel_gpu=True
    elif args.serial_cpu:
        cinfo.serial_cpu=True
    elif args.serial_gpu:
        cinfo.serial_gpu=True
    else:
        cinfo.parallel_cpu=True
    cinfo.qos=args.qos
    cinfo.account=args.account

    
    if cinfo.xstream:
        cinfo.num_gpu_per_node=16
        cinfo.num_days=2
        cinfo.gpu_partition="normal"

        
    nlam = args.nlambda
    lams = [ i/(nlam-1.) for i in range(nlam) ]
    if args.nlambda_softcore < 0:
        sclams = lams
    else:
        sclams=[ i/(args.nlambda_softcore-1.) for i in range(args.nlambda_softcore) ]
        

    refit=""
    if args.refit:
        refit="Y"
    if args.refith:
        refit="H"
    if args.refitq:
        refit="Q"
        
    if args.stdti:
        ti = sctisetup( args.parm, args.crd,
                        args.molmask1, args.molmask2,
                        args.timask1, args.timask2,
                        args.scmask1, args.scmask2,
                        refit )
        ti.execute(lams,cinfo,sclams=sclams)
        if args.test:
            ti.test()


    if args.uniti:
        ti = unitisetup( args.parm, args.crd,
                         args.molmask1, args.molmask2,
                         args.timask1, args.timask2,
                         args.scmask1, args.scmask2 )
        ti.execute(lams,cinfo,sclams=sclams)


            
    if args.serti:
        ti = sersctisetup( args.parm, args.crd,
                           args.molmask1, args.molmask2,
                           args.timask1, args.timask2,
                           args.scmask1, args.scmask2,
                           refit )
        ti.execute(lams,cinfo,sclams=sclams)
        if args.test:
            ti.test()
            
    if args.altti:
        ti = altsctisetup( args.parm, args.crd,
                           args.molmask1, args.molmask2,
                           args.timask1, args.timask2,
                           args.scmask1, args.scmask2 )
        ti.execute(lams,cinfo,sclams=sclams)
        if args.test:
            ti.test()
            
    if args.piti:
        ti = pisctisetup( args.parm, args.crd,
                          args.molmask1, args.molmask2,
                          args.timask1, args.timask2,
                          args.scmask1, args.scmask2 )
        ti.execute(lams,cinfo,sclams=sclams)
        if args.test:
            ti.test()

if __name__ == "__main__":
    maincli()
    exit(0)

    

if __name__ == "__main__":

    nlam = 12
    lams = [ i/(nlam-1.) for i in range(nlam) ]
    
    compinfo = computer_info()
    compinfo.cpu_partition="main"
    compinfo.gpu_partition="gpu"
    compinfo.num_cpu_per_node=24
    compinfo.num_gpu_per_node=4
    compinfo.nstlim=4000
    compinfo.ntpr=1000

    if True:
        
        parmfile = "lj.parm7"
        #parmfile = "wBZFBNZ.parm7"
    
        if False:
            scti = sctisetup( parmfile, "wBZFBNZ.rst7",
                              ":1", ":2",
                              ":1", ":2",
                              ":1@H3,C3,C2,H2,O1", ":2@H1,H6" )
            scti.execute(lams,compinfo)
            scti.test()
        
        if True:
            altscti = altsctisetup( parmfile, "wBZFBNZ.rst7",
                                    ":1", ":2",
                                    ":1", ":2",
                                    ":1@H3,C3,C2,H2,O1", ":2@H1,H6" )
            altscti.execute(lams,compinfo)
            #altscti.test()

        if True:

            piscti = pisctisetup( parmfile, "wBZFBNZ.rst7",
                                  ":1", ":2",
                                  ":1", ":2",
                                  ":1@H3,C3,C2,H2,O1", ":2@H1,H6" )
            #piscti.prepare_decharge(lams)
            #piscti.prepare_softcore(lams)
            #piscti.prepare_recharge(lams)
            piscti.execute(lams,compinfo)
            #piscti.test()





            
            
    if False:
        
        if False:
            scti = sctisetup( "prot.parm7", "prot.rst7",
                              ":1-3", ":4-6",
                              ":2",   ":5",
                              ":2@HD22,CD1,HD11,HD12,HD13", ":5@HG3,CE,HE2,HE3,NZ,HZ1,HZ2,HZ3" )
            #scti.prepare_decharge(lams)
            #scti.prepare_softcore(lams)
            #scti.prepare_recharge(lams)
            scti.execute(lams,compinfo)
            scti.test()

        if True:
            altscti = altsctisetup( "prot.parm7", "prot.rst7",
                                    ":1-3", ":4-6",
                                    ":2",   ":5",
                                    ":2@HD22,CD1,HD11,HD12,HD13", ":5@HG3,CE,HE2,HE3,NZ,HZ1,HZ2,HZ3" )
            #altscti.prepare_decharge(lams)
            #altscti.prepare_softcore(lams)
            #altscti.prepare_recharge(lams)
            altscti.execute(lams,compinfo)
            altscti.test()

        if True:
            piscti = pisctisetup( "prot.parm7", "prot.rst7",
                                   ":1-3", ":4-6",
                                   ":2",   ":5",
                                   ":2@HD22,CD1,HD11,HD12,HD13", ":5@HG3,CE,HE2,HE3,NZ,HZ1,HZ2,HZ3" )
            #piscti.prepare_decharge(lams)
            #piscti.prepare_softcore(lams)
            #piscti.prepare_recharge(lams)
            piscti.execute(lams,compinfo)
            piscti.test()

            #  deq      21.4947
            #   sc      -0.6814
            #  req      27.5912
            #  sum      48.4045
            
            #  deq      27.3049
            #   sc      -6.4915
            #  req      27.5912
            #  sum      48.4046
            
            #  deq      27.3049
            #   sc      -6.4915
            #  req      27.5912
            #  sum      48.4045


# p0 = OpenParm( "sc_%.8f.parm7"%(0.), xyz="sc.rst7" )
# p1 = OpenParm( "sc_%.8f.parm7"%(1.), xyz="sc.rst7" )
# lam = 0.
# dedparm = "sc.nc.dedparm"
# dvdl = softcore_tigradient(p0,p1,lam,dedparm)
# print dvdl


# ALL      : 27.3049 + 0.9048 = 28.2097
# LJ-ONLY  : -0.9048 + 0.9048 =  0.0000  OK

#!/usr/bin/python
import parmed
import argparse

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
        if abs(param.box[3]-109.471219)<1.e-4 and \
           abs(param.box[4]-109.471219)<1.e-4 and \
           abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param


def GetSelectedAtomIndices(param,maskstr):
    sele = []
    if len(maskstr) > 0:
        newmaskstr = maskstr.replace("@0","!@*")
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""Converts an arbitrary amber mask to an atom selection mask
""" )
    
    parser.add_argument \
        ("-p","--parm",
         help="parm7 file",
         type=str,
         required=True )


    parser.add_argument \
        ("-w","--whole-molecule",
         help="if present, then selecting any portion of a molecule selects the whole molecule",
         action='store_true' )

    
    parser.add_argument \
        ("-r","--byresidue",
         help="if present, then return a list of matching residues",
         action='store_true' )


    parser.add_argument \
        ("-n","--norange",
         help="if present, then list each match individually -- as opposed to using dashes to represent ranges",
         action='store_true' )

    
    parser.add_argument \
        ("-c","--crd",
         help="input restart file",
         type=str,
         default=None,
         required=False )
    
    parser.add_argument \
        ("-m","--mask",
         help="amber mask selection of atoms to polarize",
         type=str,
         required=True )
    
    args = parser.parse_args()
    

    parmfile = args.parm
    crd7file = args.crd
    maskstr  = args.mask
    maskstr = maskstr.replace('\\',"").replace('\\',"").replace('\\',"").replace('\\',"").replace('\\',"")

    if crd7file is not None:
        try:
            parm = OpenParm( parmfile, xyz=crd7file )
        except:
            parm = OpenParm( parmfile )
    else:
        parm = OpenParm( parmfile )


        
    atoms = GetSelectedAtomIndices( parm, maskstr )

    if args.whole_molecule:
       smols=[]
       iat = 0
       for imol,mnat in enumerate(parm.parm_data["ATOMS_PER_MOLECULE"]):
           for i in range(mnat):
               if iat in atoms:
                  smols.append( imol )
               iat += 1
       smols=list(set(smols))
       smols.sort()
       atoms = []
       iat=0
       for imol,mnat in enumerate(parm.parm_data["ATOMS_PER_MOLECULE"]):
           for i in range(mnat):
               if imol in smols:
                  atoms.append(iat)
               iat += 1 

    res = []
    for a in atoms:
        res.append( parm.atoms[a].residue.idx )
    res = list(set(res))
    res.sort()

    if args.byresidue:
        if args.norange:
            s = ":%s"%( ",".join( ["%i"%(x+1) for x in res ] ) )
        else:
            s = ListToSelection( res )
            s = s.replace("@",":")
    else:
        if args.norange:
            s = "@%s"%( ",".join( ["%i"%(x+1) for x in atoms ] ) )
        else:
            s = ListToSelection( atoms )
            
    print s


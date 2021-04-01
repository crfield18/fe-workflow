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
        if abs(param.box[3]-109.471219)<1.e-4 and \
           abs(param.box[4]-109.471219)<1.e-4 and \
           abs(param.box[5]-109.471219)<1.e-4:
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
       if abs(param.box[3]-109.471219)<1.e-4 and \
          abs(param.box[4]-109.471219)<1.e-4 and \
          abs(param.box[5]-109.471219)<1.e-4:
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
    if param.box is not None:
        rst.box = [param.box[0], param.box[1], param.box[2], param.box[3], param.box[4], param.box[5]]
    rst.write( "{}.rst7".format(fname) )


def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]

def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


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
        sele = [ param.atoms[i].idx for i in parmed.amber.mask.AmberMask( param, newmaskstr ).Selected() ]
    return sele


def GetSelectedResidueIndices(param,maskstr):
    a = GetSelectedAtomIndices(param,maskstr)
    b = list(set([ param.atoms[c].residue.idx for c in a ]))
    b.sort()
    return b


def ExtractFrcmod(p,sele,fname,masses=True,bonds=True,angles=True,dihedrals=True,lj=True):
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from collections import defaultdict as ddict
    import re
    
    q = Extract(p,sele)
    satoms = [ a.idx for a in q.atoms ]
    namemap = ddict(str)
    for a in satoms:
        oldname = q.atoms[a].atom_type.name
        newname = oldname
        if not "MSK" in newname:
            newname += "MSK"
        namemap[newname] = oldname
        q.atoms[a].atom_type.name = newname
        q.atoms[a].type = newname

    try:
        self = parmed.amber.parameters.AmberParameterSet.from_structure(q)
    except Exception as e:
        print("AmberParameterSet.from_structure failed because:",e)
        print("Retrying with ParameterSet.from_structure(q,allow_unequal_duplicates=True)")
        self = parmed.amber.parameters.ParameterSet.from_structure(q,allow_unequal_duplicates=True)
    angfact=0.9999995714245039
    
    fh = open(fname,"w")

    fh.write("modified parameters")
    fh.write('\n')
    # Write the atom mass
    fh.write('MASS\n')
    if masses:
        for atom, typ in iteritems(self.atom_types):
            if atom in namemap:
                fh.write('%s%11.8f\n' % (namemap[atom].ljust(6), typ.mass))

    fh.write('\n')
        
    # Write the bonds
    fh.write('BOND\n')
    if bonds:
        cdone = set()
        for (a1, a2), typ in iteritems(self.bond_types):
            typ.k = float("%.8f"%(typ.k))
            delta = 0
            if a1 in namemap and a2 in namemap:
                qq = (namemap[a1],namemap[a2])
                if qq in cdone: continue
                qq = (namemap[a2],namemap[a1])
                if qq in cdone: continue
                cdone.add(qq)
            else:
                continue    
            fh.write('%s-%s   %19.14f  %11.8f\n' %
                     (namemap[a1].ljust(2), namemap[a2].ljust(2), typ.k, typ.req))
    fh.write('\n')

    # Write the angles
    fh.write('ANGLE\n')
    if angles:
        cdone = set()
        for (a1, a2, a3), typ in iteritems(self.angle_types):
            typ.k = float("%.8f"%(typ.k))
            delta = 0.
            if a1 in namemap and a2 in namemap and \
               a3 in namemap:
                qq = (namemap[a1],namemap[a2],namemap[a3])
                if qq in cdone: continue
                qq = (namemap[a3],namemap[a2],namemap[a1])
                if qq in cdone: continue
                cdone.add(qq)
            else:
                continue    
            fh.write('%s-%s-%s   %19.14f  %17.3f\n' %
                     (namemap[a1].ljust(2), namemap[a2].ljust(2), namemap[a3].ljust(2), typ.k+delta,
                      typ.theteq * angfact))
    fh.write('\n')

    # Write the dihedrals
    fh.write('DIHE\n')
    if dihedrals:
        cdone = set()
        for (a1, a2, a3, a4), typ in iteritems(self.dihedral_types):
            isnew = False
            if a1 in namemap and a2 in namemap and \
               a3 in namemap and a4 in namemap:
                qq = (namemap[a1],namemap[a2],namemap[a3],namemap[a4])
                if qq in cdone: continue
                qq = (namemap[a4],namemap[a3],namemap[a2],namemap[a1])
                if qq in cdone: continue
                cdone.add(qq)
                isnew = True
            else:
                continue

            if isinstance(typ, DihedralType) or len(typ) == 1:
                if not isinstance(typ, DihedralType):
                    typ = typ[0]
                    typ.phi_k = float("%.8f"%(typ.phi_k))
                if abs(typ.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                    namemap[a3].ljust(2), namemap[a4].ljust(2), 1, typ.phi_k, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                    namemap[a3].ljust(2), namemap[a4].ljust(2), 1, typ.phi_k, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
            else:
                typ = sorted( typ, key=lambda x: x.per, reverse=False )
                for dtyp in typ[:-1]:
                    dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                    if abs(dtyp.phase-180) < 0.0001:
                        #print "%20.16f"%(180.0/dtyp.phase)
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                      namemap[a3].ljust(2), namemap[a4].ljust(2), 1, dtyp.phi_k,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                    else:
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                      namemap[a3].ljust(2), namemap[a4].ljust(2), 1, dtyp.phi_k,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                dtyp = typ[-1]
                dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                if abs(dtyp.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                    namemap[a3].ljust(2), namemap[a4].ljust(2), 1, dtyp.phi_k,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (namemap[a1].ljust(2), namemap[a2].ljust(2),
                                                    namemap[a3].ljust(2), namemap[a4].ljust(2), 1, dtyp.phi_k,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                    
    fh.write('\n')
    # Write the impropers
    fh.write('IMPROPER\n')
    if dihedrals:
        for (a1, a2, a3, a4), typ in iteritems(self.improper_periodic_types):
            # Make sure wild-cards come at the beginning
            if a2 == 'X':
                assert a4 == 'X', 'Malformed generic improper!'
                a1, a2, a3, a4 = a2, a4, a3, a1
            elif a4 == 'X':
                a1, a2, a3, a4 = a4, a1, a3, a2

            typ.phi_k = float("%.8f"%(typ.phi_k))
            if a1 in namemap and a2 in namemap and \
               a3 in namemap and a4 in namemap:

                if abs(typ.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %20.14f %13.3f %5.1f\n' %
                             (namemap[a1].ljust(2), namemap[a2].ljust(2), namemap[a3].ljust(2), namemap[a4].ljust(2),
                              typ.phi_k, typ.phase * angfact, typ.per))
                else:
                    fh.write('%s-%s-%s-%s %20.14f %13.8f %5.1f\n' %
                             (namemap[a1].ljust(2), namemap[a2].ljust(2), namemap[a3].ljust(2), namemap[a4].ljust(2),
                              typ.phi_k, typ.phase * angfact, typ.per))

                
    fh.write('\n')
    
    # Write the LJ terms
    fh.write('NONB\n')
    if lj:
        if True:
            for atom, typ in iteritems(self.atom_types):
                #typ.rmin = float("%.8f"%(typ.rmin))
                typ.epsilon = float("%.9f"%(typ.epsilon))
                if atom in namemap:
                    fh.write('%-3s  %12.8f %18.9f\n' %
                             (namemap[atom].ljust(2), typ.rmin, typ.epsilon))

    fh.write('\n')
    
    # Write the NBFIX terms
    if lj:
        if self.nbfix_types:
            fh.write('LJEDIT\n')
            for (a1, a2), (eps, rmin) in iteritems(self.nbfix_types):
                if a1 in namemap and a2 in namemap:
                    fh.write('%s %s %13.8f %13.8f %13.8f %13.8f\n' %
                             (namemap[a1].ljust(2), namemap[a2].ljust(2), eps, rmin/2,
                              eps, rmin/2))


                    
if __name__ == "__main__":
    
    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
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
                        help="Amber selection mask. Only residues within this mask will be renamed",
                        type=str,
                        default="@*",
                        required=False )

    parser.add_argument("-n","--name",
                        help="Rename all selected residues to this name",
                        type=str,
                        required=False)

    parser.add_argument("-e1","--extract1",
                        help="Amber selection mask. Only the parameters from the residues within this mask will be extracted",
                        type=str,
                        required=False)

    parser.add_argument("-e2","--extract2",
                        help="Amber selection mask. Only the parameters from the residues within this mask will be extracted",
                        type=str,
                        required=False)

    parser.add_argument("-o","--output",
                        help="Output file basename if -n/--name is specified",
                        type=str,
                        required=False)


    parser.add_argument("-o1","--output1",
                        help="Output file basename if -e1/--extract1 is specified",
                        type=str,
                        required=False)

    parser.add_argument("-o2","--output2",
                        help="Output file basename if -e2/--extract2 is specified",
                        type=str,
                        required=False)

    parser.add_argument("-f","--frcmod",
                        help="Output frcmod file",
                        type=str,
                        required=False)

    parser.add_argument("-l","--lib",
                        help="Output lib file",
                        type=str,
                        required=False)



    args = parser.parse_args()

    p    = OpenParm(args.parm,xyz=args.crd)

    sres = GetSelectedResidueIndices(p,args.mask)
    resmask = ":%s"%( ",".join( [ "%i"%(res+1) for res in sres ] ) )
    if args.name is not None:
        if len(args.name) > 0:
            name = args.name
            if len(name) > 3:
                name = name[0:3]
            for i in parmed.amber.mask.AmberMask( p, args.mask ).Selected():
                p.atoms[i].residue.name = name
            SaveParmRst(p, args.output)
            parmed.formats.Mol2File.write(p, "{}.mol2".format(args.output))
            parmed.formats.PDBFile.write(p,"{}.pdb".format(args.output))
            fh = open("%s.seq"%(args.output),"w")
            seq = GetResSeq( p )
            seqchunks=[]
            for chunk in divide_chunks(seq,10):
                seqchunks.append( " ".join(chunk) )
                seqstr = "\n".join(seqchunks)
            fh.write(seqstr)

            ExtractFrcmod(p,"@*","{}.frcmod".format(args.output))
            parmed.tools.writeOFF(p,"{}.lib".format(args.output)).execute()


    if args.extract1 is not None:
        if len(args.extract1) > 0:
            p1 = CopyParm (p)
            sres1 = GetSelectedResidueIndices(p1,args.extract1)
            resmask1 = ":%s"%( ",".join( [ "%i"%(res+1) for res in sres1 ] ) )
            q1    = Extract(p1,resmask1)
            parmed.formats.Mol2File.write(q1, "{}.mol2".format(args.output1))
            parmed.formats.PDBFile.write(q1,"{}.pdb".format(args.output1))

            fh = open("%s.seq"%(args.output1),"w")
            seq = GetResSeq( q1 )
            seqchunks=[]
            for chunk in divide_chunks(seq,10):
                seqchunks.append( " ".join(chunk) )
                seqstr = "\n".join(seqchunks)
            fh.write(seqstr)

            ExtractFrcmod(q1,"@*","{}.frcmod".format(args.output1))
            parmed.tools.writeOFF(q1,"{}.lib".format(args.output1)).execute()


    if args.extract2 is not None:
        if len(args.extract2) > 0:
            p2       = CopyParm (p)
            sres2    = GetSelectedResidueIndices(p2,args.extract2)
            resmask2 = ":%s"%( ",".join( [ "%i"%(res+1) for res in sres2 ] ) )
            q2    = Extract(p2,resmask2)
            parmed.formats.Mol2File.write(q2, "{}.mol2".format(args.output2))
            parmed.formats.PDBFile.write(q2,"{}.pdb".format(args.output2))

            fh = open("%s.seq"%(args.output2),"w")
            seq = GetResSeq( q2 )
            seqchunks=[]
            for chunk in divide_chunks(seq,10):
                seqchunks.append( " ".join(chunk) )
                seqstr = "\n".join(seqchunks)
            fh.write(seqstr)

            ExtractFrcmod(q2,"@*","{}.frcmod".format(args.output2))
            parmed.tools.writeOFF(q2,"{}.lib".format(args.output2)).execute()
 

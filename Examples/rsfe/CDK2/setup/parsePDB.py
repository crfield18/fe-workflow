#!/usr/bin/env python3


def CopyPDB( pdb ):
    import copy
    p = copy.copy( pdb )
    p.coordinates = copy.copy( pdb.coordinates )
    return p


def Strip( pdb, mask ):
    p = CopyPDB( pdb )
    p.strip( "%s"%(mask) )
    return p


def Extract( pdb, mask ):
    return Strip( pdb, "!(%s)"%(mask) )

def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]

def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


def Joincom( mol1, mol2, nonligand ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    elif nonligand.coordinates.shape[0] == 0:
        pq=nonligand
    else:
        pq = mol1 + mol2 + nonligand
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates, nonligand.coordinates) )
        return pq

def Joinaq( mol1, mol2 ):
    import numpy
    import copy
    if mol1.coordinates.shape[0] == 0:
        pq=mol1
    elif mol2.coordinates.shape[0] == 0:
        pq=mol2
    else:
        pq = mol1 + mol2
        pq.coordinates = numpy.concatenate( (mol1.coordinates, mol2.coordinates) )
        return pq


if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    from sys import exit

    parser = argparse.ArgumentParser     ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Takes as input - a PDB file, an atommap file, two ligand masks \"orglig\" and \"tarlig\", name of the two ligand mol2 files as \"orgmol2\" and \"tarmol2\", and output basename. Splits PDB into water.pdb, nonwater.pdb, mol1.pdb containing \"orglig\", and mol2.pdb containing a new residue \"tarlig\" which consists of a subset of \"orglig\" atoms that are mapped as per provided atommap file. Merges the different pdbs into a single pdb and generates the scmasks and timasks" )

    parser.add_argument("-p1","--pdb1",
                        help="PDB file 1",
                        type=str,
                        required=True)

    parser.add_argument("-p2","--pdb2",
                        help="PDB file 2",
                        type=str,
                        required=True)

    parser.add_argument("-m","--map",
                        help="Atom mapfile",
                        type=str,
                        required=True )

    parser.add_argument("-l1","--lig1",
                        help="name of ligand 1 in PDB file 1" ,
                        type=str,
                        required=False)

    parser.add_argument("-l2","--lig2",
                        help="name of ligand 2 in PDB file 2" ,
                        type=str,
                        required=False)

    parser.add_argument("-m1","--lig1mol2",
                        help="name of mol2 file of ligand 1" ,
                        type=str,
                        required=False)

    parser.add_argument("-m2","--lig2mol2",
                        help="name of mol2 file of ligand 2" ,
                        type=str,
                        required=False)

    parser.add_argument("-s","--system",
                        help="sys either complex (com) or aqueous (aq)" ,
                        type=str,
                        required=False)

    parser.add_argument("-o","--output",
                        help="output basename of merged files, sc and ti masks" ,
                        type=str,
                        required=False)


    args = parser.parse_args()

    f1 = open ("{}".format(args.map),'r')
    lines = f1.readlines()

    mappedmol1atomnames=[]
    mappedmol2atomnames=[]
    for line in lines:
        atoms = re.split('=>|\n',line.replace(" ",""))
        mappedmol1atomnames.append(atoms[0])
        mappedmol2atomnames.append(atoms[1])

    if ".mol2" in args.pdb1:
        q = parmed.load_file(args.pdb1)
        q.save("{}.pdb".format(args.pdb1.split('.mol2')[0]),overwrite=True)
        pdb1 = parmed.load_file("{}.pdb".format(args.pdb1.split('.mol2')[0]))
    else:
        pdb1 = parmed.load_file(args.pdb1)

    if ".mol2" in args.pdb2:
        q = parmed.load_file(args.pdb2)
        q.save("{}.pdb".format(args.pdb2.split('.mol2')[0]),overwrite=True)
        pdb2 = parmed.load_file("{}.pdb".format(args.pdb2.split('.mol2')[0]))
    else:
        pdb2 = parmed.load_file(args.pdb2)


#### manipulate pdb1 
    if args.system == "com":
        pdb1copy1 = CopyPDB (pdb1)
        p1mol1 = Extract(pdb1copy1, ":{}".format(args.lig1))
        p1mol2 = CopyPDB (p1mol1)
        pdb1copy2 = CopyPDB (pdb1)
        p1nl = Strip(pdb1copy2, ":{}".format(args.lig1))
    else:
        p1mol1 = CopyPDB (pdb1)
        p1mol2 = CopyPDB (pdb1)

    p1mol2.residues[0].name = "{}".format(args.lig2)

    p1mol1atomnames=[]
    for a in p1mol1.atoms:
        p1mol1atomnames.append(a.name)

    atomstodel=set (p1mol1atomnames) ^ set (mappedmol1atomnames)

    if len(atomstodel) > 0:
    	atoms=','.join([str(i) for i in atomstodel])
    	atoms="".join(('@',atoms))
    	p1mol2 = Strip(p1mol2,atoms)

    p1mol1tomol2map = []
    for j in p1mol2.atoms:
        count=0
        for i in mappedmol1atomnames:
            if ( j.name == i ):
                p1mol1tomol2map.append(count)
            count+=1

    count=0
    for j in p1mol2.atoms:
        j.name = mappedmol2atomnames[p1mol1tomol2map[count]]
        count+=1
###################

#### manipulate pdb2
    if args.system == "com":
        pdb2copy1 = CopyPDB (pdb2)
        p2mol1 = Extract(pdb2copy1, ":{}".format(args.lig2))
        p2mol2 = CopyPDB (p2mol1)
        pdb2copy2 = CopyPDB (pdb2)
        p2nl = Strip(pdb2copy2, ":{}".format(args.lig2))
    else:
        p2mol1 = CopyPDB (pdb2)
        p2mol2 = CopyPDB (pdb2)

    p2mol2.residues[0].name = "{}".format(args.lig1)

    p2mol1atomnames=[]
    for a in p2mol1.atoms:
        p2mol1atomnames.append(a.name)

    atomstodel=set (p2mol1atomnames) ^ set (mappedmol2atomnames)

    if len(atomstodel) > 0:
    	atoms=','.join([str(i) for i in atomstodel])
    	atoms="".join(('@',atoms))
    	p2mol2 = Strip(p2mol2,atoms)

    p2mol1tomol2map = []
    for j in p2mol2.atoms:
        count=0
        for i in mappedmol2atomnames:
            if ( j.name == i ):
                p2mol1tomol2map.append(count)
            count+=1

    count=0
    for j in p2mol2.atoms:
        j.name = mappedmol1atomnames[p2mol1tomol2map[count]]
        count+=1
###################


    #parmed.tools.writeCoordinates(mol1,"mol1.pdb").execute()
    #parmed.tools.writeCoordinates(mol2,"mol2.pdb").execute()
    #parmed.tools.writeCoordinates(nl,"nonligand.pdb").execute()


    if args.system == "com":
        merge1 = Joincom(p1mol1,p1mol2,p1nl)
        merge2 = Joincom(p2mol2,p2mol1,p2nl)
    elif args.system == "aq":
        merge1 = Joinaq(p1mol1,p1mol2)
        merge2 = Joinaq(p2mol2,p2mol1)

    parmed.tools.writeCoordinates(merge1,"merged1_{}.pdb".format(args.output)).execute()
    parmed.tools.writeCoordinates(merge2,"merged2_{}.pdb".format(args.output)).execute()

    p1fh = open("merged1_%s.seq"%(args.output),"w")
    seq = GetResSeq( merge1 )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    p1fh.write(seqstr)
    p1fh.close()

    p2fh = open("merged2_%s.seq"%(args.output),"w")
    seq = GetResSeq( merge2 )
    seqchunks=[]
    for chunk in divide_chunks(seq,10):
        seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
    p2fh.write(seqstr)
    p2fh.close()


    m1names = []
    m2names = []
    m1 = parmed.load_file(args.lig1mol2)
    m2 = parmed.load_file(args.lig2mol2)
    for j in m1.atoms:
        m1names.append(j.name)

    for j in m2.atoms:
        m2names.append(j.name)

    scatoms1 = set (m1names) ^ set (mappedmol1atomnames)
    scatoms2 = set (m2names) ^ set (mappedmol2atomnames)

    if len(scatoms1) > 0:
        scmask1 = ','.join([str(i) for i in scatoms1])
        scmask1 = "".join((":{}@".format(args.lig1),scmask1))
    else:
        scmask1 = '""'

    if len(scatoms2) > 0:
        scmask2 = ','.join([str(i) for i in scatoms2])
        scmask2 = "".join((":{}@".format(args.lig2),scmask2))
    else:
        scmask2 = '""'

    timask1 = ":{}".format(args.lig1)
    timask2 = ":{}".format(args.lig2)

    fh1 = open("{}.scmask1".format(args.output), "w")
    fh1.write(scmask1)
    fh1.close

    fh2 = open("{}.scmask2".format(args.output), "w")
    fh2.write(scmask2)
    fh2.close

    fh3 = open("{}.timask1".format(args.output), "w")
    fh3.write(timask1)
    fh3.close

    fh4 = open("{}.timask2".format(args.output), "w")
    fh4.write(timask2)
    fh4.close



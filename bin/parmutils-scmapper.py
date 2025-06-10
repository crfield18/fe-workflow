#!/usr/bin/env python3
import parmed
from io import StringIO
from rdkit import rdBase
import rdkit.Chem


def mcss(mol2str_1, mol2str_2, maxtime=60, isotope_map=None, selec=''):
    """
    Maximum common substructure search via RDKit/fmcs.

    :param mol2str_1: first MOL2 string
    :type mol2str_1: string
    :param mol2str_2: second MOL2 string
    :type mol2str_2: string
    :param maxtime: timeout for fmcs in seconds
    :type maxtime: float
    :param isotope_map: explicit user atom mapping
    :type isotope_map: dict
    :param selec: selection method for multiple MCS
    :type selec: string
    :raises: SetupError
    :returns: index map
    :rtype: dict
    """

    from collections import OrderedDict, defaultdict
    #import openbabel as ob

    import sys



    _fmcs_imp = 'c++'                       # 'python' or 'c++'
    if _fmcs_imp == 'c++':
        from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare

        # RDKit 2015.03.1 FMCS C++ implentation, seems not to be exactly the
        # same implementation e.g. SMARTS string more specific
        # NOTE: some different parameters and order!
        #       matchChiralTag = False not implemented before 2015.03.1
        _params = dict(maximizeBonds = False, threshold = 1.0,
                       verbose = False, matchValences = False,
                       ringMatchesRingOnly = True, completeRingsOnly = True,
                       bondCompare = BondCompare.CompareAny)
    else:
        from rdkit.Chem.MCS import FindMCS

        # defaults are minNumAtoms = 2, maximize = 'bonds',
        # atomCompare = 'elements', bondCompare = 'bondtypes',
        # matchValences = False, ringMatchesRingOnly = False,
        # completeRingsOnly = False, timeout = None, threshold = None
        #
        # when completeRingsOnly = True also ringMatchesRingOnly = True
        #
        # if number of atoms in match < minNumAtoms then smarts will be None
        # completeRingsOnly = True disallows partial rings
        # MCS algorithm is exhaustive so use timeout to limit time
        _params = dict(minNumAtoms = 2, maximize = 'atoms', atomCompare = 'elements',
                       bondCompare = 'any', matchValences = False,
                       ringMatchesRingOnly = True, completeRingsOnly = True,
                       threshold = None)



    logger = sys.stderr
    # disable warning about no explicit hydrogens
    rdBase.DisableLog('rdApp.warning')
    mol1 = rdkit.Chem.MolFromMol2Block(mol2str_1, sanitize = False,
                                       removeHs = False)
    mol2 = rdkit.Chem.MolFromMol2Block(mol2str_2, sanitize = False,
                                       removeHs = False)
    rdBase.EnableLog('rdApp.warning')

    _params.update(timeout = int(maxtime) )

    # FIXME: test c++ implementation
    if isotope_map:
        if _fmcs_imp == 'c++':
            _params.update(atomCompare = AtomCompare.CompareIsotopes)
        else:
            _params.update(atomCompare = 'isotopes')

        max_idx1 = mol1.GetNumAtoms()
        max_idx2 = mol2.GetNumAtoms()

        icnt = 0

        # NOTE: would it make sense to have multiple atoms of a molecule tagged
        #       as the same isotope?
        for idx1, idx2 in isotope_map.items():
            icnt += 1

            logger.write('Mapping atom index %i to %i' % (idx1, idx2) )

            if idx2 < 1:
                atom1 = mol1.GetAtomWithIdx(idx1-1)
                atom1.SetIsotope(-1 * idx2 + icnt)
            else:
                if idx1 > max_idx1 or idx2 > max_idx2 or idx1 < 0 or idx2 < 0:
                    logger.write('Error: indices out of bounds (%i, %i)' %
                                 (max_idx1, max_idx2) )
                    raise Exception('Mapping indices out of bounds '
                                    '(%i, %i)' % (max_idx1, max_idx2) )

                # FIXME: guard against non-existing indices
                atom1 = mol1.GetAtomWithIdx(idx1-1)
                atom1.SetIsotope(icnt)

                atom2 = mol2.GetAtomWithIdx(idx2-1)
                atom2.SetIsotope(icnt)
    else:
        if _fmcs_imp == 'c++':
            _params.update(atomCompare = AtomCompare.CompareAny)
        else:
            _params.update(atomCompare = 'any')

        n_chiral1 = len(rdkit.Chem.FindMolChiralCenters(mol1) )
        n_chiral2 = len(rdkit.Chem.FindMolChiralCenters(mol2) )

        if n_chiral1 > 0:
            logger.write('Warning: state 0 has %i chiral center%s. Check if '
                         'configurations are inverted!'
                         % (n_chiral1, 's' if n_chiral1 > 1 else '') )
        if n_chiral2 > 0:
            logger.write('Warning: state 1 has %i chiral center%s. Check if '
                         'configurations are inverted!'
                         % (n_chiral2, 's' if n_chiral2 > 1 else '') )


    mcs = FindMCS( (mol1, mol2), **_params)

    if _fmcs_imp == 'c++':
        smarts = mcs.smartsString
        completed = not mcs.canceled
    else:
        smarts = mcs.smarts
        completed = mcs.completed

    logger.write('Running RDKit/fmcs (%s implementation) with arguments:\n%s\n' %
                 (_fmcs_imp,
                  ', '.join(['%s=%s' % (k,v) for k,v in _params.items()] ) ) )

    if not smarts:
        raise Exception('No MCSS match could be found')

    if not completed:
        logger.write('Warning: MCSS timed out after %.2fs' % maxtime)

    p = rdkit.Chem.MolFromSmarts(smarts)

    #conv = ob.OBConversion()
    #conv.SetInAndOutFormats('mol2', 'mol2')
    ## NOTE: this relies on a modified Openbabel MOL2 writer
    #conv.AddOption('r', ob.OBConversion.OUTOPTIONS)  # do not append resnum
    #obmol1 = ob.OBMol()
    #errlev = ob.obErrorLog.GetOutputLevel()
    #ob.obErrorLog.SetOutputLevel(0)
    #conv.ReadString(obmol1, mol2str_1)
    #ob.obErrorLog.SetOutputLevel(errlev)

    # NOTE: experimental!
    if selec == 'spatially-closest':
        m1 = mol1.GetSubstructMatches(p, uniquify=False, maxMatches=100, useChirality=False)
        m2 = mol2.GetSubstructMatches(p, uniquify=False, maxMatches=100, useChirality=False)

        logger.write('Applying spatially-closest algorithm (%s, %s matches)\n' %
                     (len(m1), len(m2) ) )


        # FIXME: is it possible that the smaller one has more then one matches
        #        when uniquify=True?
        if len(m1) < len(m2):
            m1, m2 = m2, m1
            conf1 = mol2.GetConformer()
            conf2 = mol1.GetConformer()
            swapped = True
        else:
            conf1 = mol1.GetConformer()
            conf2 = mol2.GetConformer()
            swapped = False

        # for x in range(m1)...
        #     match1 = m1[x]
        #     for y in range(m2)...
        #         match2 = m2[y]
        #         find sumd
        #         save x, y with smallest sumd
        mind = 999999.0
        minxy = [-1,-1]
        #DEVTHRESHOLD = 0.2**2# That's 0.2 Angstrom
        for x in range(len(m1)):
            match1 = m1[x]
            for y in range(len(m2)):
                match2 = m2[y]
                sumd = 0.0
                for i, idx1 in enumerate(match1):
                    pos1 = conf1.GetAtomPosition(idx1)
                    idx2 = match2[i]
                    pos2 = conf2.GetAtomPosition(idx2)
                    d2 = (pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 +\
                         (pos1.z - pos2.z)**2
                    sumd += d2
                #print (x,y,sumd,mind)
                if sumd < mind:
                    mind = sumd
                    minxy = [x,y]
        if swapped:
            mapping = dict(zip(m2[minxy[1]], m1[minxy[0]]))
        else:
            # JM bug 11/17 ? 
            #mapping = dict(zip(m1[minxy[1]], m2[minxy[0]]))
             mapping = dict(zip(m1[minxy[0]], m2[minxy[1]]))
        #print (mapping)

        #dist_sum = []
        #for match1 in m1:
        #    dists = []
        #
        #    for i, idx1 in enumerate(match1):
        #        pos1 = conf1.GetAtomPosition(idx1)
        #        for match2 in m2:
        #            idx2 = match2[i]
        #            pos2 = conf2.GetAtomPosition(idx2)
        #
        #            d = math.sqrt((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 +
        #                          (pos1.z - pos2.z)**2)
        #            #print pos1.x, pos2.x
        #            #print idx1, idx2, d
        #            sumd += d
        #            dists.append(d)
        #    dist_sum.append(sum(dists))
        #    print sum(dists), m1
        #min_idx = dist_sum.index(min(dist_sum) )
        #print (min_idx)
        ## FIXME: assume there is only one match for 2nd molecule
        #if swapped:
        #    mapping = dict(zip(m2[0], m1[min_idx]) )
        #else:
        #    mapping = dict(zip(m1[min_idx], m2[0]) )
        #
        #m1 = m1[min_idx]
        #m2 = m2[0]
    else:
        m1 = mol1.GetSubstructMatch(p)
        m2 = mol2.GetSubstructMatch(p)

        mapping = dict(zip(m1, m2) )

    # FIXME: we may have to reconsider this and understand when rings have
    #        to be assumed "broken"
    #
    # delete atoms from mapping that are also part of an map-external ring
    if False:
        ring_info1 = mol1.GetRingInfo()
        ring_info2 = mol2.GetRingInfo()

        rings1 = ring_info1.AtomRings()
        rings2 = ring_info2.AtomRings()

        map1 = set(mapping.keys())
        map2 = set(mapping.values())

        for ring in rings1:
            if not set(ring).issubset(map1):
                for idx in map1:
                    if idx in ring and ring_info1.NumAtomRings(idx) == 1:
                        del(mapping[idx])

        delete_values = []
        for ring in rings2:
            if not set(ring).issubset(map2):
                for idx in map2:
                    if idx in ring and ring_info2.NumAtomRings(idx) == 1:
                        delete_values.append(idx)

        mapping = {k: v for k, v in mapping.items() if v not in delete_values}


    # delete_atoms = []

    # for atom in ob.OBMolAtomIter(obmol1):
    #     idx = atom.GetIdx() - 1

    #     if idx not in mapping:
    #         delete_atoms.append(atom)

    # obmol1.BeginModify()

    # for idx in delete_atoms:
    #     obmol1.DeleteAtom(idx)

    # obmol1.EndModify()

#    conv.WriteFile(obmol1, const.MCS_MOL_FILE)

#    with open(const.MCS_MAP_FILE, 'wb') as pkl:
#        pickle.dump(mapping.keys(), pkl, 0)
#        pickle.dump(mapping.values(), pkl, 0)

    return mapping



def MutateMap(mol2file1,mol2file2):
    import parmed
    import copy
    
    if isinstance(mol2file1,str):
        mol1 = parmed.load_file(mol2file1,structure=False)
    elif isinstance(mol2file1,parmed.modeller.residue.ResidueTemplate):
        mol1 = copy.deepcopy(mol2file1)
    else:
        raise TypeError("mol2file1 must be a string or a parmed residue template")

        
    for a in mol1.atoms:
        for elem,num in parmed.periodic_table.AtomicNum.items():
            if num == a.atomic_number:
                a.type = elem
                break
    mol2str_1 = StringIO()
    mol1.save(mol2str_1,format="MOL2")
    mol2str_1 = mol2str_1.getvalue()

    print(mol2file2)
    if isinstance(mol2file2,str):
        mol2 = parmed.load_file(mol2file2,structure=False)
    elif isinstance(mol2file2,parmed.modeller.residue.ResidueTemplate):
        mol2 = copy.deepcopy(mol2file2)
    else:
        raise Exception("mol2file2 must be a string or a parmed residue template")
    print(mol2)
    for a in mol2.atoms:
        for elem,num in parmed.periodic_table.AtomicNum.items():
            if num == a.atomic_number:
                a.type = elem
                break
    mol2str_2 = StringIO()
    mol2.save(mol2str_2,format="MOL2")
    mol2str_2 = mol2str_2.getvalue()

    
    map1to2 = mcss(mol2str_1, mol2str_2, maxtime=60, isotope_map=None, selec='')

    map2to1 = dict()
    for k in map1to2:
        map2to1[ map1to2[k] ] = k

    return map1to2,map2to1


def AtomsBeyondBondDriver(p,keepatom,delatom,s=[]):
    if delatom in s:
        return s
    else:
        s.append(delatom)
    for a in p.atoms[delatom].bond_partners:
        if a.idx == keepatom:
            pass
        elif a.idx in s:
            pass
        else:
            ss = AtomsBeyondBondDriver(p,delatom,a.idx,s=s)
            for b in ss:
                if b not in s:
                    s.append(b)
    return s

def AtomsBeyondBond(p,keepatom,delatom):
    return AtomsBeyondBondDriver(p,keepatom,delatom,s=[])

def FindIdxByName(p,name):
    idx=None
    for a in p.atoms:
        if a.name == name:
            idx = a.idx
            break
    return idx

def NameMapToIdxMap(p1,p2,nmap):
    from collections import defaultdict as ddict
    i1toi2 = ddict(int)
    for m1 in nmap:
        m2 = nmap[m1]
        #print(m1,FindIdxByName(p1,m1))
        #print(m2,FindIdxByName(p2,m2))
        i1toi2[ FindIdxByName(p1,m1) ] = FindIdxByName(p2,m2)
    return i1toi2


def PartitionAcrossAtom(p,idx):
    from collections import defaultdict as ddict
    parts=ddict(list)
    for b in p.atoms[idx].bond_partners:
        parts[b.idx] = AtomsBeyondBond(p,idx,b.idx)
    return parts

def ChooseSCAtoms(p,m1to2,idx):
    from collections import defaultdict as ddict
    parts = PartitionAcrossAtom(p,idx)
    cnts = ddict(int)
    maxcnt = -1
    for b in parts:
        commoncnt = 0
        for a in parts[b]:
            if a in m1to2:
                commoncnt += 1
        cnts[b] = commoncnt
        #print(cnts[b],parts[b])
        if commoncnt > maxcnt:
            maxcnt = commoncnt
    scparts = []
    for b in parts:
        if cnts[b] < maxcnt:
            scparts.extend(parts[b])
    return scparts


def ModifiedMCSSModel1(mol1file,mol2file,mapfile):
    import parmed
    import sys
    from collections import defaultdict as ddict

    mol1 = parmed.load_file(mol1file,structure=True)
    mol2 = parmed.load_file(mol2file,structure=True)

    m1tom2,m2tom1 = MutateMap(mol1file,mol2file)
    i1toi2 = m1tom2 #NameMapToIdxMap(mol1,mol2,m1tom2)
    i2toi1 = m2tom1 #NameMapToIdxMap(mol2,mol1,m2tom1)

#    for a in mol1.atoms:
#        sys.stdout.write("%4s :"%(a.name))
#        for b in a.bond_partners:
#            sys.stdout.write(" %4s"%(b.name))
#        sys.stdout.write("\n")
#    sys.stdout.write("\n")

    sc1 = []
    for a in mol1.atoms:
        if a.idx not in i1toi2:
            sc1.append( a.idx )
    sc2 = []
    for a in mol2.atoms:
        if a.idx not in i2toi1:
            sc2.append( a.idx )

            
    for i1 in i1toi2:
        i2 = i1toi2[i1]
        a1 = mol1.atoms[i1]
        a2 = mol2.atoms[i2]
        hybrid1 = len( a1.bond_partners )
        hybrid2 = len( a2.bond_partners )
        if hybrid1 != hybrid2:
            sc1.extend( ChooseSCAtoms(mol1,m1tom2,i1) )
            sc2.extend( ChooseSCAtoms(mol2,m2tom1,i2) )

    for i1 in i1toi2:
        i2 = i1toi2[i1]
        z1 = mol1.atoms[i1].atomic_number
        z2 = mol2.atoms[i2].atomic_number
        if (z1 == 1 or z2 == 1) and z1 + z2 > 2:
            if z1 == 1:
                base1 = mol1.atoms[i1].bond_partners[0].idx
                base2 = i1toi2[base1]
            else:
                base2 = mol2.atoms[i2].bond_partners[0].idx
                base1 = i2toi1[base2]
            sc1.extend( AtomsBeyondBond(mol1,base1,i1) )
            sc2.extend( AtomsBeyondBond(mol2,base2,i2) )

            
    sc1 = list(set(sc1))
    sc2 = list(set(sc2))
    #print([mol1.atoms[a].name for a in sc1])
    #print([mol2.atoms[a].name for a in sc2])
        
    #ats = AtomsBeyondBond(mol1,FindIdxByName(mol1,"C1"),FindIdxByName(mol1,"C3"))

    for sc in sc1:
        if sc in i1toi2:
            del i1toi2[sc]
    for sc in sc2:
        if sc in i2toi1:
            del i2toi1[sc]

    if mapfile is not None: 
        fh=open(mapfile,"w")
        for i1 in i1toi2:
            i2 = i1toi2[i1]
            n1 = mol1.atoms[i1].name
            n2 = mol2.atoms[i2].name
            fh.write("%4s => %4s\n"%(n1,n2))
        fh.close()
        
    return i1toi2,i2toi1




def ModifiedMCSSModel2(mol1file,mol2file,mapfile):
    import parmed
    import sys
    from collections import defaultdict as ddict

    mol1 = parmed.load_file(mol1file,structure=True)
    mol2 = parmed.load_file(mol2file,structure=True)

    m1tom2,m2tom1 = MutateMap(mol1file,mol2file)
    i1toi2 = m1tom2 #NameMapToIdxMap(mol1,mol2,m1tom2)
    i2toi1 = m2tom1 #NameMapToIdxMap(mol2,mol1,m2tom1)

#    for a in mol1.atoms:
#        sys.stdout.write("%4s :"%(a.name))
#        for b in a.bond_partners:
#            sys.stdout.write(" %4s"%(b.name))
#        sys.stdout.write("\n")
#    sys.stdout.write("\n")

    sc1 = []
    for a in mol1.atoms:
        if a.idx not in i1toi2:
            sc1.append( a.idx )
    sc2 = []
    for a in mol2.atoms:
        if a.idx not in i2toi1:
            sc2.append( a.idx )

            
    for i1 in i1toi2:
        i2 = i1toi2[i1]
        a1 = mol1.atoms[i1]
        a2 = mol2.atoms[i2]
        hybrid1 = len( a1.bond_partners )
        hybrid2 = len( a2.bond_partners )
        z1 = a1.atomic_number
        z2 = a2.atomic_number
        if z1 != z2:
            sc1.extend( ChooseSCAtoms(mol1,m1tom2,i1) )
            sc2.extend( ChooseSCAtoms(mol2,m2tom1,i2) )
            sc1.append( i1 )
            sc2.append( i2 )
        elif (hybrid1 != hybrid2):
            sc1.extend( ChooseSCAtoms(mol1,m1tom2,i1) )
            sc2.extend( ChooseSCAtoms(mol2,m2tom1,i2) )
        
            
    for i1 in i1toi2:
        i2 = i1toi2[i1]
        z1 = mol1.atoms[i1].atomic_number
        z2 = mol2.atoms[i2].atomic_number
        if (z1 == 1 or z2 == 1) and z1 + z2 > 2:
            if z1 == 1:
                base1 = mol1.atoms[i1].bond_partners[0].idx
                base2 = i1toi2[base1]
            else:
                base2 = mol2.atoms[i2].bond_partners[0].idx
                base1 = i2toi1[base2]
            sc1.extend( AtomsBeyondBond(mol1,base1,i1) )
            sc2.extend( AtomsBeyondBond(mol2,base2,i2) )

            
    sc1 = list(set(sc1))
    sc2 = list(set(sc2))
    #print([mol1.atoms[a].name for a in sc1])
    #print([mol2.atoms[a].name for a in sc2])
        
    #ats = AtomsBeyondBond(mol1,FindIdxByName(mol1,"C1"),FindIdxByName(mol1,"C3"))

    for sc in sc1:
        if sc in i1toi2:
            del i1toi2[sc]
    for sc in sc2:
        if sc in i2toi1:
            del i2toi1[sc]

    if mapfile is not None: 
        fh=open(mapfile,"w")
        for i1 in i1toi2:
            i2 = i1toi2[i1]
            n1 = mol1.atoms[i1].name
            n2 = mol2.atoms[i2].name
            fh.write("%4s => %4s\n"%(n1,n2))
        fh.close()
        
    return i1toi2,i2toi1



def MCSSModel(mol1file,mol2file,mapfile):
    import parmed
    import sys
    from collections import defaultdict as ddict

    mol1 = parmed.load_file(mol1file,structure=True)
    mol2 = parmed.load_file(mol2file,structure=True)

    i1toi2,i2toi1 = MutateMap(mol1file,mol2file)

    if mapfile is not None:
        fh=open(mapfile,"w")
        for i1 in i1toi2:
            i2 = i1toi2[i1]
            n1 = mol1.atoms[i1].name
            n2 = mol2.atoms[i2].name
            fh.write("%4s => %4s\n"%(n1,n2))
        fh.close()

    return i1toi2,i2toi1



def GraphMCSS(method,graph):
    from collections import defaultdict as ddict
    data = ddict(list)
    pairs = ddict( lambda: ddict( str ) )
    tree = []
    gin = open(graph,"r")
    for line in gin:
        mols = line.strip().split()
        if len(mols) > 1:
            tree.append(mols)
    gin.close()
        
    for mols in tree:
        mol1file = mols[0]
        for mol2file in mols[1:]:
            m1to2,m2to1 = method(mol1file,mol2file)
            pairs[mol1file][mol2file] = m1to2,m2to1
            
            if mol1file not in data:
                data[mol1file] = [ a for a in m1to2 ]
            else:
                data[mol1file] = [ a for a in data[mol1file] if a in m1to2 ]
            if mol2file not in data:
                data[mol2file] = [ a for a in m2to1 ]
            else:
                data[mol2file] = [ a for a in data[mol2file] if a in m2to1 ]

    molfiles = ddict(str)
    for mol in data:
        molfiles[mol] = parmed.load_file(mol,structure=True)
                
    for mols in tree:
        mol1file = mols[0]
        for mol2file in mols[1:]:
            m1to2,m2to1 = pairs[mol1file][mol2file]
            t = ddict(str)
            for a in data[mol1file]:
                t[a] = m1to2[a]
            m1to2 = t

            pname = "%s~%s.map.txt"%( mol1file.replace(".mol2",""), mol2file.replace(".mol2",""))
            fout = open(pname,"w")
            for m1 in m1to2:
                m2 = m1to2[m1]
                n1 = molfiles[mol1file].atoms[m1].name
                n2 = molfiles[mol2file].atoms[m2].name
                fout.write("%4s => %4s\n"%(n1,n2))
            fout.close()

            
    
if __name__ == "__main__":

    import argparse
    from collections import defaultdict as ddict

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""Writes a alchemical transformation map file.
Mapping method 0: MCSS
Mapping method 1: Extends MCSS to include chains from atoms changing hybridization into SC region
Mapping method 2: Like (1), but also includes atoms that change atomic number (and chains connected to them)""")

    parser.add_argument \
        ("-t","--method",
         help="Mapping method. 0=MCSS, 1=MCSS+dHyb, 2=MCSS+Hyb+dZ",
         type=int,
         default=2,
         required=False )

    parser.add_argument \
        ("-o","--map",
         help="Name of the output map file",
         type=str,
         required=False )

    parser.add_argument \
        ("-a","--lig1",
         help="mol2 file for ligand 1",
         type=str,
         default=None,
         required=False )
    
    parser.add_argument \
        ("-b","--lig2",
         help="mol2 file for ligand 2",
         type=str,
         default=None,
         required=False )

    parser.add_argument \
        ("-g", "--graph",
         help="graph of mol2 filenames. Each row is a node and each column is an edge",
         type=str,
         default=None,
         required=False)
    
    
    args = parser.parse_args()
    print(args.method)
    print(args.lig1,args.lig2,args.map,args.graph)

    if args.graph is not None:
        if args.lig1 is not None or args.lig2 is not None or args.map is not None:
            raise Exception("If --graph is used, then --lig1, --lig2, and --map must all be used")
    else:
        if args.lig1 is None or args.lig2 is None or args.map is None:
            raise Exception("If --graph is not used, then none of --lig1, --lig2, nor --map can be used")

    if args.graph is None:
        if args.method == 0:
            MCSSModel(args.lig1,args.lig2,args.map)
        elif args.method == 1:
            ModifiedMCSSModel1(args.lig1,args.lig2,args.map)
        elif args.method == 2:
            ModifiedMCSSModel2(args.lig1,args.lig2,args.map)
    else:
        if args.method == 0:
            method = lambda lig1,lig2 : MCSSModel(lig1,lig2,None)
        elif args.method == 1:
            method = lambda lig1,lig2 : ModifiedMCSSModel1(lig1,lig2,None)
        elif args.method == 2:
            method = lambda lig1,lig2 : ModifiedMCSSModel2(lig1,lig2,None)
            
        GraphMCSS(method,args.graph)

    
    

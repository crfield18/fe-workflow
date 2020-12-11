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

        
    for a in mol1.atoms:
        for elem,num in parmed.periodic_table.AtomicNum.items():
            if num == a.atomic_number:
                a.type = elem
    mol2str_1 = StringIO()
    mol1.save(mol2str_1,format="MOL2")
    mol2str_1 = mol2str_1.getvalue()

    
    if isinstance(mol2file2,str):
        mol2 = parmed.load_file(mol2file2,structure=False)
    elif isinstance(mol2file2,parmed.modeller.residue.ResidueTemplate):
    #else:
        mol2 = copy.deepcopy(mol2file2)
    
    for a in mol2.atoms:
        for elem,num in parmed.periodic_table.AtomicNum.items():
            if num == a.atomic_number:
                a.type = elem
    mol2str_2 = StringIO()
    mol2.save(mol2str_2,format="MOL2")
    mol2str_2 = mol2str_2.getvalue()

    
    map1to2 = mcss(mol2str_1, mol2str_2, maxtime=60, isotope_map=None, selec='')

    map2to1 = dict()
    for k in map1to2:
        map2to1[ map1to2[k] ] = k

    return map1to2,map2to1



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

def GetSelectedAtomIndices(param,maskstr):
    #param = parmed.load_file(parmfile)
    #mask = parmed.amber.mask.AmberMask( param, maskstr )
    #aidxs = mask.Selected()
    #for aidx in aidxs:
    #    atom = param.atoms[aidx]
    #    res  = atom.residue
    sele = []
    if len(maskstr) > 0:
        sele = [ param.atoms[i].idx for i in parmed.amber.mask.AmberMask( param, maskstr ).Selected() ]
    return sele

def GetSelectedResidueIndices(param,maskstr):
    a = GetSelectedAtomIndices(param,maskstr)
    b = list(set([ param.atoms[c].residue.idx for c in a ]))
    b.sort()
    return b

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
    sele = ""
    if len(sarr) > 0:
        sele = ",".join(sarr)
    return sele

def AtomListToSelection(p,atomlist):
    from collections import defaultdict as ddict
    reslist = ddict( list )
    haslist = ddict( list )
    atomlist = list(set(atomlist))
    for iat in atomlist:
        ires = p.atoms[iat].residue.idx
        if not ires in reslist:
            reslist[ires] = [ a.idx for a in p.residues[ires].atoms ]
            haslist[ires] = []
        reslist[ires].remove(iat)
        haslist[ires].append(iat)

    ress = []
    atms = []
    for ires in reslist:
        if len(reslist[ires]) == 0:
            ress.append(ires)
        else:
            for iatm in haslist[ires]:
                atms.append(iatm)
    rsel = ListToSelection(ress)
    asel = ListToSelection(atms)
    sel = ""
    if len(rsel) > 0:
        sel = ":"+rsel
        if len(asel) > 0:
            sel += "|"
    if len(asel) > 0:
        sel += "@"+asel
    return sel


def GetMoleculeMask( p, target_sele ):
    import copy
    ires = GetSelectedResidueIndices(p,target_sele)[0]

    apermol = copy.deepcopy(p.parm_data["ATOMS_PER_MOLECULE"])
    for i in range(1,len(apermol)):
        apermol[i] += apermol[i-1]

    first_atom = p.residues[ires].atoms[0].idx
    imol = 0
    for i in range(len(apermol)):
        if first_atom < apermol[i]:
            imol = i
            break
    if imol == 0:
        a1 = 0
        a2 = apermol[0]-1
    else:
        a1 = apermol[imol-1]
        a2 = apermol[imol]-1
            
    molmask = "@%i-%i"%(a1+1,a2+1)
    return molmask



def GetResSeq( parm ):
    rtc=parmed.modeller.residue.ResidueTemplateContainer.from_structure( parm )
    return [r.name for r in rtc]




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
        
        nfile = open(native_frcmod,"w")
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
        


def divide_chunks_generator(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

def divide_chunks(l,n):
    return list(divide_chunks_generator(l,n))


def Mutate(parmfile,rstfile,timask1,mol22,lib2,frcmod2,base,args=None,mapdict=None):
    p = OpenParm( parmfile, xyz=rstfile )

    #
    # convert the timask1 into a mol2/structure instance
    #
    # p1 is a parameter file containing only the timask1
    p1 = Extract( p, timask1 )
    if len(p1.atoms) == 0:
        print("Mask \"%s\" does not match any atoms"%(timask1))
        exit(1)
    # write the timask1 residue to a mol2 file (internally, a string)
    mol21str = StringIO()
    parmed.formats.Mol2File.write(p1,mol21str)
    # re-read the string as a mol2 file
    m21 = parmed.formats.Mol2File.parse( StringIO( mol21str.getvalue() ), structure=False )

    
    tisel1 = GetSelectedAtomIndices(p,timask1)
    if len(tisel1) != len(m21.atoms):
        print("Error: timask \"%s\" contains %i atoms, but %s has %i atoms"%(timask1,len(tisel1),mol21,len(m21.atoms)))
        exit(1)

    err=False
    for i in range(len(m21.atoms)):
        if m21.atoms[i].residue.name != p.atoms[tisel1[i]].residue.name:
            print("Error: Residue name mismatch mol2 @%s (%s) and parm :%i@%s (%s)"%(m21.atoms[i].name,m21.atoms[i].residue.name,p.atoms[tisel1[i]].residue.idx,p.atoms[tisel1[i]].name,p.atoms[tisel1[i]].residue.name))
            err=True
        if m21.atoms[i].name != p.atoms[tisel1[i]].name:
            print("Error: Atom name mismatch mol2 @%s (%s) and parm :%i@%s (%s)"%(m21.atoms[i].name,m21.atoms[i].residue.name,p.atoms[tisel1[i]].residue.idx,p.atoms[tisel1[i]].name,p.atoms[tisel1[i]].residue.name))
            err=True

    if err:
        exit(1)
            
    m22 = parmed.load_file(mol22,structure=False)

    if mapdict is not None:
        map1to2 = dict()

        for on in mapdict:
            mn = mapdict[on]

            k1=-1
            for i,a in enumerate(m21.atoms):
                if a.name == on:
                    k1 = i
            if k1 >= 0:
                k2=-1
                for i,a in enumerate(m22.atoms):
                    if a.name == mn:
                        k2 = i
                if k2 >= 0:
                    map1to2[k1]=k2
                else:
                    raise Exception("Found %s in target mask, but failed to find %s in mol2"%(on,mn))
            else:
                raise Exception("Failed to find %s in target mask"%(on))
  
        
        map2to1 = dict()
        for k in map1to2:
            map2to1[ map1to2[k] ] = k

    else:
        map1to2,map2to1 = MutateMap(m21,mol22)

        fh = open("%s.map.txt"%(base),"w")
        for k1 in map1to2:
            k2 = map1to2[k1]
            fh.write("%-4s  =>  %-4s\n"%(m21.atoms[k1].name,m22.atoms[k2].name))
        fh.close()
        
        
        
    err=False
    rn1 = list(set( [ a.residue.name for a in m21.atoms ] ))
    rn2 = list(set( [ a.residue.name for a in m22.atoms ] ))
    for r2 in rn2:
        if r2 in rn1:
            print("Error: residue %s found in %s and %s"%(r2,mol21,mol22))
            err=True
    if err:
        exit(1)

    sc1=[]
    for a in range(len(m21.atoms)):
        if a not in map1to2:
            sc1.append(a)
    sc2=[]
    for a in range(len(m22.atoms)):
        if a not in map2to1:
            sc2.append(a)
        

    molsele = GetSelectedAtomIndices(p,GetMoleculeMask(p,timask1))
    tisele  = GetSelectedAtomIndices(p,timask1)
    old_mol1_first = molsele[0]
    old_mol1_last  = molsele[-1]
    old_ti1_first  = tisele[0]
    old_ti1_last   = tisele[-1]
    for i in range(len(tisele)):
        if tisele[i] != old_ti1_first + i:
            print("Error: TI selection is not a contiguous region")
            exit(1)
    
    new_mol1_first = 0
    new_mol1_last  = (new_mol1_first+len(molsele)-1)
    new_mol2_first = len(molsele)
    new_mol2_last  = (new_mol2_first+len(molsele)-1) - len(m21.atoms) + len(m22.atoms)

    new_ti1_first = old_ti1_first-old_mol1_first+new_mol1_first
    new_ti1_last  = new_ti1_first + len(m21.atoms) - 1
    new_ti2_first = old_ti1_first-old_mol1_first+new_mol2_first
    new_ti2_last  = new_ti2_first + len(m22.atoms) - 1

    scmask1 = ListToSelection( [ new_ti1_first + i for i in sc1 ] )
    if len(scmask1) > 0:
        scmask1 = "@" + scmask1
    scmask2 = ListToSelection( [ new_ti2_first + i for i in sc2 ] )
    if len(scmask2) > 0:
        scmask2 = "@" + scmask2

    option_str = "--molmask1=\"@%i-%i\" --timask1=\"@%i-%i\" --scmask1=\"%s\" --molmask2=\"@%i-%i\" --timask2=\"@%i-%i\" --scmask2=\"%s\""%(new_mol1_first+1,new_mol1_last+1,new_ti1_first+1,new_ti1_last+1,scmask1,new_mol2_first+1,new_mol2_last+1,new_ti2_first+1,new_ti2_last+1,scmask2)


    m21 = parmed.formats.Mol2File.parse( StringIO( mol21str.getvalue() ), structure=True )
    m22 = parmed.load_file(mol22,structure=True)


    
    for a in range(len(m21.atoms)):
        m21.atoms[a].xx = p1.atoms[a].xx
        m21.atoms[a].xy = p1.atoms[a].xy
        m21.atoms[a].xz = p1.atoms[a].xz
    
    parmed.tools.actions.writeCoordinates( m21, "%s.mol1.pdb"%(base) ).execute()

    c1 = []
    for a in range(len(m21.atoms)):
        c1.extend( [m21.atoms[a].xx,m21.atoms[a].xy,m21.atoms[a].xz] )

    c2 = []
    for a in range(len(m22.atoms)):
        c2.extend( [m22.atoms[a].xx,m22.atoms[a].xy,m22.atoms[a].xz] )

    c1to2 = [None] * len(m21.atoms)
    for a in range(len(m21.atoms)):
        if a in map1to2:
            c1to2[a] = map1to2[a]

    c2wts = [0]*len(m22.atoms)
    c2to1 = [None] * len(m22.atoms)
    for a in range(len(m22.atoms)):
        if a in map2to1:
            b = map2to1[a]
            c2to1[a] = b
            c2wts[a] = max(m22.atoms[a].atomic_number,m21.atoms[b].atomic_number)
            if m22.atoms[a].atomic_number != m21.atoms[b].atomic_number:
                c2wts[a] = c2wts[a]*4



    #
    #c2 = RmsOverlay(c2,c1,mol2refmap=c2to1,molw=c2wts)
    #
    
    for a in range(len(m22.atoms)):
        if a in map2to1:
            b = map2to1[a]
            for k in range(3):
                c2[k+a*3] = c1[k+b*3]

    for a in range(len(m22.atoms)):
        m22.atoms[a].xx = c2[0+a*3]
        m22.atoms[a].xy = c2[1+a*3]
        m22.atoms[a].xz = c2[2+a*3]


        
    mol0parm = Extract(p,"@"+ListToSelection(molsele))
    parmed.tools.actions.writeCoordinates(mol0parm,"%s.mol0.pdb"%(base)).execute()
    
    mol1parm = CopyParm(mol0parm)
    if len(sc1) > 0:
        mol1parm.strip( "@" + ListToSelection( [ i + old_ti1_first - old_mol1_first for i in sc1 ] ) )

    ii=-1
    for i in range(len(m22.atoms)):
        j = c2to1[i]
        if j is None:
            continue
        else:
            ii += 1
        a = mol1parm.atoms[old_ti1_first-old_mol1_first + ii]
        a.name = m22.atoms[i].name
        a.type = m22.atoms[i].type
        a.charge = m22.atoms[i].charge
        a.atomic_number = m22.atoms[i].atomic_number
        a.residue.name = m22.atoms[i].residue.name
        a.xx = c2[0+i*3]
        a.xy = c2[1+i*3]
        a.xz = c2[2+i*3]

    parmed.tools.actions.writeCoordinates(mol1parm,"%s.mol1.pdb"%(base)).execute()

    
    parmed.tools.writeOFF( p, "%s.lib"%(base) ).execute()
    WriteFrcmod( p, "%s.frcmod"%(base), uniqueparams=False )

    
    res2=m22.atoms[0].residue.name
    if lib2 is None:
        lib2 = "%s.lib"%(res2)
    if frcmod2 is None:
        frcmod2 = "%s.frcmod"%(res2)
    
    fh = open("%s.sh"%(base),"w")
    fh.write("""#!/bin/bash

if [ ! -e "%s" ]; then
   parmchk2 -i %s -o %s -f mol2
fi
"""%(frcmod2,mol22,frcmod2))

    fh.write("""
if [ -e "%s" ]; then
if [ "$(grep -l ATTN %s)" != "" ]; then
   echo "The frcmod file "%s" contains unset parameters"
   echo "You need to supply a properly parameterized frcmod file to timutate"
   echo "or you'll need to modify %s.sh"
   exit 1  
fi
fi
"""%(frcmod2,frcmod2,frcmod2,base))

    source="leaprc.gaff"
    if args is not None:
        source = args.source
    water="leaprc.water.tip4pew"
    if args is not None:
        if args.tip3p:
            water="leaprc.water.tip3p"
    
    fh.write("""
if [ ! -e "%s" ]; then
   cat <<EOF > %s.lib.inp
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source %s
source %s
%s = loadmol2 %s
saveoff %s %s
quit
EOF
   tleap -s -f %s.lib.inp
   rm leap.log
   rm %s.lib.inp
   sed -i 's/ZTQ/%s/g' %s
fi
"""%(lib2,base,source,water,"ZTQ",mol22,"ZTQ",lib2,base,base,res2,lib2))

    mol0seq = GetResSeq( mol0parm )
    mol1seq = GetResSeq( mol1parm )

    mol0chunks=[]
    for chunk in divide_chunks(mol0seq,10):
        mol0chunks.append( " ".join(chunk) )
    mol0str = "\n".join(mol0chunks)

    mol1chunks=[]
    for chunk in divide_chunks(mol1seq,10):
        mol1chunks.append( " ".join(chunk) )
    mol1str = "\n".join(mol1chunks)

    
#     fh.write("""
# cat << 'EOF' > %s.sh.cmds
# source leaprc.protein.ff14SB
# source leaprc.RNA.OL3
# source %s
# source %s
# loadOff %s.lib
# loadAmberParams %s.frcmod
# loadOff %s
# loadAmberParams %s

# mol0 = loadPdbUsingSeq %s.mol0.pdb { %s }
# mol1 = loadPdbUsingSeq %s.mol1.pdb { %s }

# """%(base,source,water,base,base,lib2,frcmod2,base,mol0str,base,mol1str))



    fh.write("""
cat << 'EOF' > %s.sh.cmds
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source %s
source %s
loadOff %s.lib
loadAmberParams %s.frcmod
loadOff %s
loadAmberParams %s

mol0 = loadPdb %s.mol0.pdb
mol1 = loadPdb %s.mol1.pdb 

"""%(base,source,water,base,base,lib2,frcmod2,base,base))


    

    
    tparm = CopyParm( p )
    tparm.strip( "(%s)|:WAT"%( AtomListToSelection(p,molsele) ) )
    nonwat_len = len(tparm.atoms)
    if nonwat_len > 0:
        parmed.tools.writeCoordinates(tparm, "%s.nonwater.pdb"%(base)).execute()
        seq = GetResSeq( tparm )
        seqchunks=[]
        for chunk in divide_chunks(seq,10):
            seqchunks.append( " ".join(chunk) )
        seqstr = "\n".join(seqchunks)
        fh.write("nonwater = loadPdbUsingSeq %s.nonwater.pdb { %s }\n\n"%(base,seqstr))
        #fh.write("nonwater = loadPdb %s.nonwater.pdb\n\n"%(base))
         

    tparm = CopyParm( p )
    tparm.strip( "(%s)|(!:WAT)"%( AtomListToSelection(p,molsele) ) )
    wat_len = len(tparm.atoms)
    if wat_len > 0:
        parmed.tools.writeCoordinates(tparm, "%s.solvent.pdb"%(base)).execute()
        fh.write("solvent = loadPdb %s.solvent.pdb\n\n"%(base))

    
    fh.write("\n\nx = combine { mol0 mol1")
    if nonwat_len > 0:
        fh.write(" nonwater")
    if wat_len > 0:
        fh.write(" solvent")
    fh.write(" }\n\n")
    
    if p.box is not None:
        fh.write("setbox x centers\n")

    incdir=""
    if args is not None:
        if len(args.include) > 0:
            incdir="-I %s"%(args.include)
        
    fh.write("""
saveAmberParm x %s.parm7 %s.rst7
quit
EOF

tleap -s %s -f %s.sh.cmds | grep -v "+---" | grep -v "+Currently" > %s.sh.out
cat %s.sh.out

"""%(base,base,incdir,base,base,base))

    fh.write("if [ ! -e %s.parm7 -o ! -e %s.rst7 ]; then echo \"Failed to make parameter file\"; exit 1; fi\n"%(base,base))

    if p.box is not None:
        fh.write("""
# Set the box dimensions
ChBox -X %.12f -Y %.12f -Z %.12f -al %.12f -bt %.12f -gm %.12f -c %s.rst7 -o %s.rst7.tmp; mv %s.rst7.tmp %s.rst7
"""%(p.box[0],p.box[1],p.box[2],p.box[3],p.box[4],p.box[5],base,base,base,base))

        if abs(p.box[-1]-109.471219000000) < 1.e-4:
            fh.write("""
# Reset the ifbox flag
sed -i 's|0       0       1|0       0       2|' %s.parm7

"""%(base))

            

    fh.write("\n\n\n################################\n\n\n")

    if args is None:
        fh.write("# ")
    if args.altti:
        fh.write("parmutils-tigen.py --altti -p %s.parm7 -c %s.rst7 \\\n"%(base,base))
    elif args.piti:
        fh.write("parmutils-tigen.py --piti -p %s.parm7 -c %s.rst7 \\\n"%(base,base))
    elif args.uniti:
        fh.write("parmutils-tigen.py --uniti -p %s.parm7 -c %s.rst7 \\\n"%(base,base))
    elif args.serti:
        fh.write("parmutils-tigen.py --serti -p %s.parm7 -c %s.rst7 \\\n"%(base,base))
    else:
        fh.write("parmutils-tigen.py --stdti -p %s.parm7 -c %s.rst7 \\\n"%(base,base))
    if args is None:
        fh.write("# ")
    fh.write("%s \\\n"%(option_str))
    if args is not None:
        fh.write("--nlambda=%i \\\n"%(args.nlambda))
        fh.write("--nlambda-softcore=%i \\\n"%(args.nlambda_softcore))

        if args.parallel_cpu:
            fh.write("--parallel-cpu \\\n")
        elif args.serial_cpu:
            fh.write("--serial-cpu \\\n")
        elif args.parallel_gpu:
            fh.write("--parallel-gpu \\\n")
        elif args.serial_gpu:
            fh.write("--serial-gpu \\\n")

        if args.oversubscribe:
            fh.write("--oversubscribe \\\n")
            
        if args.refitq:
            fh.write("--refitq \\\n")
        elif args.refith:
            fh.write("--refith \\\n")
        elif args.refit:
            fh.write("--refit \\\n")

        if args.gti:
            fh.write("--gti \\\n")
        fh.write("--nmropt=%i \\\n"%(args.nmropt))
        fh.write("--cut=%.2f \\\n"%(args.cut))

            
        fh.write("--nstlim=%i \\\n"%(args.nstlim))
        fh.write("--numexchg=%i \\\n"%(args.numexchg))
        fh.write("--ntpr=%i \\\n"%(args.ntpr))
        fh.write("--ntwx=%i \\\n"%(args.ntwx))
        fh.write("--cpu-partition=\"%s\" \\\n"%(args.cpu_partition))
        fh.write("--cpus-per-node=%i \\\n"%(args.cpus_per_node))
        fh.write("--gpu-partition=\"%s\" \\\n"%(args.gpu_partition))
        fh.write("--gpus-per-node=%i \\\n"%(args.gpus_per_node))
        fh.write("--qos=\"%s\" \\\n"%(args.qos))
        fh.write("--account=\"%s\" \\\n"%(args.account))
        fh.write("--exclude=\"%s\" \\\n"%(args.exclude.replace("\"","")))
        fh.write("--constraint=\"%s\" \\\n"%(args.constraint.replace("\"","")))
        fh.write("--min-nodes=%i \\\n"%(args.nodes))
        fh.write("--max-nodes=%i \\\n"%(args.nodes))
        fh.write("--days=%i \\\n"%(args.days))
        fh.write("--launch=\"%s\" \\\n"%(args.launch))
        fh.write("--steps=%i \\\n"%(args.steps))
        fh.write("--steps-per-slurm=%i \\\n"%(args.steps_per_slurm))

    fh.write("\n\n")
    fh.write("#POST\n\n")
    fh.close()

    
    import os
    import stat
    st = os.stat("%s.sh"%(base))
    os.chmod("%s.sh"%(base), st.st_mode | stat.S_IEXEC)

    # num2elem = {}
    # for elem,num in parmed.periodic_table.AtomicNum.items():
    #     num2elem[num]=elem
        
    # fh = file("c1.xyz","w")
    # fh.write("%i\n\n"%(len(m21.atoms)))
    # for i,a in enumerate(m21.atoms):
    #     fh.write("%4s %12.5f %12.5f %12.5f\n"%(num2elem[a.atomic_number],c1[0+i*3],c1[1+i*3],c1[2+i*3]))
    # fh.close()
 
    # fh = file("c2.xyz","w")
    # fh.write("%i\n\n"%(len(m22.atoms)))
    # for i,a in enumerate(m22.atoms):
    #     fh.write("%4s %12.5f %12.5f %12.5f\n"%(num2elem[a.atomic_number],c2[0+i*3],c2[1+i*3],c2[2+i*3]))
    # fh.close()

    

    
    

def RmsOverlay(mol,ref,mol2refmap=None,molw=None):
    import numpy
    import copy
    
    nm = len(mol)//3
    nr = len(ref)//3
    
    if molw == None:
        molw = [1.] * nm
    if len(molw) != nm:
        raise Exception("molw size mismatch %i (expected %i)"%(len(molw),nm))
    if mol2refmap is None:
        mol2refmap = []
        for i in range(nm):
            if i < nr:
                mol2refmap.append(i)
            else:
                mol2refmap.append(None)
                molw[i] = 0
    if len(mol2refmap) != nm:
        raise Exception("mol2refmap size mismatch %i (expected %i)"%(len(mol2refmap),nm))

    
    totw = 0
    for i in range(nm):
        if mol2refmap[i] is not None:
            totw += molw[i]
            
    for i in range(nm):
        molw[i] = molw[i] / float(totw)

        
    molcom = [0,0,0]
    for i in range(nm):
        if mol2refmap[i] is not None:
            for k in range(3):
                molcom[k] += molw[i]*mol[k+i*3]
      
    refcom = [0,0,0]
    for i in range(nm):
        j = mol2refmap[i]
        if j is not None:
            for k in range(3):
                refcom[k] += molw[i]*ref[k+j*3]

    sumsq = 0
    fm = [0]*9
    for iat in range(nm):
        jat = mol2refmap[iat]
        if jat is not None:
            for i in range(3):
                c1 = mol[i+iat*3]-molcom[i]
                c2 = ref[i+jat*3]-refcom[i]
                sumsq += molw[iat] * ( c1*c1+c2*c2 )
                for j in range(3):
                    c1 = mol[j+iat*3]-molcom[j]
                    fm[j+i*3] += molw[iat] * c1 * c2

    #fm = [ i for i in range(9) ]
    fmm = numpy.array(fm).reshape(3,3)
    #print fmm
    U,S,VT = numpy.linalg.svd(fmm)
    detU = numpy.linalg.det(U)
    detV = numpy.linalg.det(VT)
    dot = S.sum()
    if detU*detV < 0:
        jmin=0
        if S[1] < S[jmin]:
            jmin=1
        if S[2] < S[jmin]:
            jmin=2
        dot -= 2*S[jmin]
        U[jmin,:] = -U[jmin,:]

    rmsPD = numpy.sqrt( max( 0, sumsq - 2*dot ) )
    rot = numpy.matmul( U, VT ) #numpy.transpose(VT), numpy.transpose(U) )

    ocrd = numpy.array( mol ).reshape(nm,3)
    molcom = numpy.array(molcom)
    refcom = numpy.array(refcom)
    rcrd = [0.]*(3*nm)
    #print molcom
    #print refcom
    for i in range(nm):
        #ocrd[i,:] = ocrd[i,:]-molcom+refcom
        ocrd[i,:] = numpy.dot( rot, ocrd[i,:]-molcom ) + refcom
        for k in range(3):
            rcrd[k+i*3] = ocrd[i,k]
    return rcrd

    
if __name__ == "__main__":

    import argparse
    from collections import defaultdict as ddict

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""
      """,
      epilog="""

""")


    parser.add_argument \
        ("--map",
         help="if present, then the argument is a file that maps the target atom names to the mutated atom names using the format oldname => newname with one entry-per-line",
         type=str,
         default="",
         required=False )
    
    
    parser.add_argument \
        ("-p","--parm",
         help="amber parm7 file of state 0",
         type=str,
         required=True )
    
    parser.add_argument \
        ("-c","--crd",
         help="amber restart file of state 0",
         type=str,
         required=True )

      
    parser.add_argument \
        ("--target",
         help="amber mask of the residue to be mutated, e.g., \":1\" for residue 1",
         type=str,
         required=True )

    parser.add_argument \
        ("--mol2",
         help="mol2 file of the state 1 residue",
         type=str,
         required=True )

    parser.add_argument \
        ("--lib",
         help="(optional) parm off/lib file for the state 1 residue. If not specified, a lib file will be generated from the mol2 file using tleap",
         type=str,
         default=None,
         required=False )

    parser.add_argument \
        ("--frcmod",
         help="(optional) frcmod file for the state 1 residue. If not specified, a frcmod file will be generated using parmchk2",
         type=str,
         default=None,
         required=False )

    ####

    parser.add_argument \
        ("--stdti",
         help="prepare a standard TI pathway calculations",
         action='store_true' )

    parser.add_argument \
        ("--uniti",
         help="prepare a unified protocol TI pathway calculation",
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
         help="prepare a series deletion/addition pathway calculation",
         action='store_true' )

    parser.add_argument \
        ("-b","--base",
         help="base filename for system preparation [default: ticopy]",
         type=str,
         default="ticopy",
         required=False )
    
    parser.add_argument \
        ("-n","--nlambda",
         help="number of states [default: 12]",
         type=int,
         default=12,
         required=False )
    
    parser.add_argument \
        ("--nlambda-softcore",
         help="number of states for softcore. If -1, then same as --nlambda [default: -1]",
         type=int,
         default=-1,
         required=False )


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
         help="if present, then set the gti_add_sc=2 flag in TEMPLATE.mdin",
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
         default=10.0,
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
        ("--oversubscribe",
         help="Run all lambdas at the same time on all GPUs from 1 node",
         action="store_true",
         default=False,
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
        ("--nodes",
         help="number of CPU nodes to use (if less than 1, then it is set to nlambda) [default: -1]",
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
        ("--steps-per-slurm",
         help="production cycles per slurm file [default: 1]",
         type=int,
         default=1,
         required=False )

    parser.add_argument \
        ("--steps",
         help="total number of production cycles [default: 1]",
         type=int,
         default=1,
         required=False )


    parser.add_argument \
        ("--tip3p",
         help="use tip3p rather than tip4p/ew",
         action="store_true",
         default=False,
         required=False )

    
    parser.add_argument \
        ("-I","--include",
         help="tleap include directory",
         type=str,
         default="",
         required=False )

    
    parser.add_argument \
        ("-S","--source",
         help="leaprc to source [default: leaprc.gaff]",
         type=str,
         default="leaprc.gaff",
         required=False )


    args = parser.parse_args()


    mapdict=None
    if len(args.map) > 0:
        mapdict = ddict(str)
        fh = open(args.map,"r")
        for line in fh:
            cols = line.strip().split()
            if len(cols) == 3:
                mapdict[cols[0]]=cols[2]
                #print("%s => %s"%(cols))
            #print("%3i => %3i"%(i+1,map1to2[i]+1))

        
    Mutate(args.parm,args.crd,
           args.target,
           args.mol2,
           args.lib,
           args.frcmod,
           args.base,
           args=args,
           mapdict=mapdict)
        


    
    #print MutateMap(mol1,"mol2_dir/hexane.mol2")
    
    # scmask1_list = []
    # for a in range(len(mol1.atoms)):
    #     if a not in map1to2:
    #         scmask1_list.append(a)

    # scmask2_list=[]
    # for a in range(len(mol2.atoms)):
    #     if a not in map2to1:
    #         scmask2_list.append(a)

    # print "scmask1 serial "," ".join(["%i"%(i+1) for i in scmask1_list ])
    # print "scmask2 serial "," ".join(["%i"%(i+1) for i in scmask2_list ])

    # for a in range(len(mol2.atoms)):
    #     if a in map2to1:
    #         b = map2to1[a]
    #         print "%4s => %4s"%(mol2.atoms[a].name,mol1.atoms[b].name)

            
    

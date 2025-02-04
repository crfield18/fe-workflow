#!/usr/bin/env python3
import numpy as np
import rdkit
import parmed
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, rdDetermineBonds
import io
from glob import glob


def show_atom_number(mol, label, offset):
    for atom in mol.GetAtoms():
        atom.SetProp(label, str(atom.GetIdx()+1+offset))
    return mol


def GetSelectedAtomIndices(param,maskstr):
    import parmed
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


class Ligand:
    def __init__(self, input, parmfile, rstfile, image_loc, edge_name, showidxs, size):
        """ This is a class to generate an image of a ligand with softcore atoms highlighted

        Parameters
        ----------
        input : str
            Path to the input mdin file [mdin]
        parmfile : str
            Path to the parameter file [parm7]
        rstfile : str
            Path to the restart file [rst7]
        image_loc : str
            Path to the directory where the images will be saved
        edge_name : str
            Name of the edge (e.g. 0~1)
        showidxs : bool
            If True, show the atom indexes in the image
        size : tuple of int
            The image size in pixels.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If the input, parmfile, or rstfile does not exist

        
        """
        self.input = Path(input)
        self.parmfile = Path(parmfile)
        self.rstfile = Path(rstfile)    
        self.softcore = []
        self.timask = []
        self.image_loc = Path(image_loc)
        self.edge_name = edge_name
        self.showidxs = showidxs
        self.size = size
        if not self.input.exists():
            raise FileNotFoundError(f"{self.input} does not exist")
        if not self.parmfile.exists():
            raise FileNotFoundError(f"{self.parmfile} does not exist")
        if not self.rstfile.exists():
            raise FileNotFoundError(f"{self.rstfile} does not exist")

    def run(self):
        """ This function runs the class and generates the images """
        self.parm = self._read_parm()
        self._find_softcore()
        ti_only, sc_only, ti_all = [], [], []
        for i in range(len(self.softcore)):
            ti_only.append(self._find_difference(self.softcore[i], self.timask[i]))
            sc_only.append(self._parse_range(self.softcore[i]))
            ti_all.append(self._parse_range(self.timask[i]))

        self.ti_only = ti_only
        self.sc_only = sc_only
        self.ti_all = ti_all
        self._write_2d_structure()
        return

    def _read_parm(self):
        """ This function reads the parameter and restart files and returns a parameter file object"""
        parm = None
        try:
            parm = parmed.load_file(str(self.parmfile))
            parm.load_rst7(str(self.rstfile))
        except Exception as e:
            raise RuntimeError(f"Failed to read parameter or restart file: {e}")
        return parm
    
    def _write_2d_structure(self):
        """ This function writes the 2D structure of the ligand with the softcore atoms highlighted 
        
        Parameters
        ----------
        parm : parmed.Structure
            A parmed structure object
            
        Returns
        -------
        None
        
        TODO
        ----
        This function should be modified to check charges on the ligands. 
        
        """
        
        
        ligand_resnames = [res.name for res in self.parm.residues[0:2]]
        fake_pdb = FakeFile()
        self.parm.write_pdb(fake_pdb)
        for i, mol in enumerate(ligand_resnames):
            mol_content = fake_pdb.getvalue(mask=[mol])
            rdmol = Chem.MolFromPDBBlock(mol_content, removeHs=False)

            ##for iat,atom in enumerate(self.parm.residues[i].atoms):
            ##    rdmol.GetAtomWithIdx(iat).SetFormalCharge(round(atom.charge))
            
            q1 = round(sum([atom.charge for atom in self.parm.residues[i].atoms]))

            rdDetermineBonds.DetermineBonds(rdmol,charge=q1)
            
            #q2 = Chem.GetFormalCharge(rdmol)
            #print(rdmol,
            #      "charge=",q1,q2,
            #      "natoms",len(rdmol.GetAtoms()),len(self.parm.residues[i].atoms))
            Chem.rdDepictor.SetPreferCoordGen(True)

            Chem.rdDepictor.Compute2DCoords(rdmol)
            sc = [idx - min(self.ti_all[i]) for idx in self.sc_only[i]]

            #img = Draw.MolToImage(rdmol, highlightAtoms=sc, kekulize=True, wedgeBonds=True,size=size)
            #node_names = self.edge_name.split('~')
            #img.save(f"{self.image_loc}/{self.edge_name}_{node_names[i]}.png")
            
            img = Draw.rdMolDraw2D.MolDraw2DCairo(*self.size)
            opts = img.drawOptions()
            opts.baseFontSize = 0.46
            opts.useMolBlockWedging = True
            opts.singleColourWedgeBonds = True
            opts.highlightBondWidthMultiplier = 12

            img.DrawMolecule(rdmol,highlightAtoms=sc)
            img.FinishDrawing()
            node_names = self.edge_name.split('~')
            with open(f"{self.image_loc}/{self.edge_name}_{node_names[i]}.png",'wb') as f:
                f.write(img.GetDrawingText())

              
            if self.showidxs:

                
                offset = self.parm.residues[i].atoms[0].idx
                #show_atom_number(rdmol,'molAtomMapNumber')
                show_atom_number(rdmol,'atomLabel',offset)
                img = Draw.rdMolDraw2D.MolDraw2DCairo(*self.size)
                opts = img.drawOptions()
                opts.baseFontSize = 0.42
                opts.useMolBlockWedging = True
                opts.singleColourWedgeBonds = True
                opts.highlightBondWidthMultiplier = 12

                img.DrawMolecule(rdmol,highlightAtoms=sc)
                img.FinishDrawing()
                node_names = self.edge_name.split('~')
                #print(f"{self.edge_name}_{node_names[i]} : {offset}")
                with open(f"{self.image_loc}/{self.edge_name}_{node_names[i]}_idxs.png",'wb') as f:
                    f.write(img.GetDrawingText())
            
        return 

    def _find_softcore(self):
        """ This function finds the softcore atoms from the mdin file """
        with self.input.open() as f:
            for line in f:
                if "scmask" in line:
                    line = line.split('!')[0].strip()
                    key,value = line.split('=')
                    value=value.strip()
                    self.softcore.append(value)
                if "timask" in line:
                    line = line.split('!')[0].strip()
                    key,value = line.split('=')
                    value=value.strip()
                    self.timask.append(value)

    def _parse_range(self, mask):
        value=mask.strip()
        if value[0] == "\"" and value[1] == "\"":
            value = value[1:-1]
        if value[0] == "'":
            value = value[1:]
            if value[-1] == "'":
                value = value[:-1]
        atomids = [idx+1 for idx in GetSelectedAtomIndices(self.parm,value)]
        return atomids
        
                    
    # def _parse_range(self, range_atoms):
    #     """ This function parses the range atoms and returns a set of atom ids
        
    #     This is meant to take a range of atoms and return a set of atom ids. For example,
    #     if the input is '1-3,5,7-9', the output will be {1, 2, 3, 5, 7, 8, 9}
        
    #     Parameters
    #     ----------
    #     range_atoms : str
    #         A string of range atoms
        
    #     Returns
    #     -------
    #     set
    #         A set of atom ids

    #     """
    #     numbers = set()
    #     for r in range_atoms.split(','):
    #         if '-' not in r:
    #             numbers.add(int(r))
    #             continue
    #         start, end = map(int, r.split('-'))
    #         numbers.update(range(start, end + 1))
    #     return numbers
    
    def _find_difference(self, list1, list2):
        """ This function finds the difference between two sets """
        set1 = set()
        set2 = set()
        set1.update(self._parse_range(list1))
        set2.update(self._parse_range(list2))
        return sorted(set2-set1)
    
class FakeFile:
    def __init__(self):
        """ This is a fake file class to capture the output of parmed.write_pdb """
        self.content = ""

    def write(self, text):
        """ This function writes the text to the content 
        
        Parameters
        ----------
        text : str
            The text to write to the content

        """
        self.content += text

    def getvalue(self, mask=None):
        """ This function returns the content of the file 
        
        Parameters
        ----------
        mask : list
            A list of strings to filter the content
            
        Returns
        -------
        str
            The content of the file

        """
        content = []
        if mask == None:
            return self.content
        for line in self.content.split('\n'):
            for elem in mask:
                if elem in line:
                    content.append(line)
        content = '\n'.join(content)
        return content
    
class GenerateImages:
    def __init__(self, system='1Y27', image_dir="analysis/images", sub_dir="aq", showidxs=False,size=(-1,-1)):
        """ This class generates images of the ligands with the softcore atoms highlighted 
        
        Parameters
        ----------
        system : str
            The system to analyze
        image_dir : str
            The directory to save the images
        
        Returns
        -------
        None

        """
        self.system = system
        self.find_edges()
        self.image_dir = Path(image_dir)
        self.sub_dir = sub_dir
        self.showidxs = showidxs
        self.size = size
        if not self.image_dir.exists():
            self.image_dir.mkdir(parents=True, exist_ok=True)
        self.generate_images()
        return
    
    def find_edges(self):
        """ This function finds the edges to analyze """
        self.edges = []
        for file in glob(f"{self.system}/unified/run/*~*"):
            edge = Path(file)
            if edge.is_dir():
                self.edges.append(edge)
        return
    def generate_images(self):
        """ This function generates the images """
        for edge in self.edges:
            edge_name = str(edge).split('/')[-1]
            print(f"Working on {edge_name}")
            ligand = Ligand(edge / f"{self.sub_dir}/inputs/0.00000000_ti.mdin", 
                            edge / f"{self.sub_dir}/unisc.parm7", 
                            edge / f"{self.sub_dir}/stateA.rst7", 
                            image_loc = self.image_dir,
                            edge_name=edge_name,
                            showidxs=self.showidxs,
                            size=self.size)
            ligand.run()
        
        # This section generates an html file to view the images in a browser
        # This is useful for quick visualization of the images. 
        if len(self.edges) > 0:
            import xml.etree.ElementTree as ET
            import html as HTML
            import xml.dom.minidom as md

            html = ET.Element('html', attrib={'lang':'en'})
            head = ET.SubElement(html,'head')
            ET.SubElement(head,'title').text = self.system
            ET.SubElement(head,'meta',
                          attrib={ 'http-equiv': "content-type",
                                   'content': "text/html; charset=utf-8" })

            body = ET.SubElement(html,'body')
            table = ET.SubElement(body,'table',attrib={'style': 'border: 1px solid black; border-collapse:collapse'})
            tr = ET.SubElement(table,'tr')
            th = ET.SubElement(tr,'th',attrib={'colspan': '2'})
            th.text = self.system
            for iedge,edge in enumerate(self.edges):
                edge_name = str(edge).split('/')[-1]
                lig1,lig2 = edge_name.split('~')
                
                tr = ET.SubElement(table,'tr')
                th = ET.SubElement(tr,'th',attrib={'colspan': '2'})
                th.text = edge_name
                tr = ET.SubElement(table,'tr')
                th = ET.SubElement(tr,'th')
                th.text = lig1
                th = ET.SubElement(tr,'th')
                th.text = lig2
                
                tr  = ET.SubElement(table,'tr',attrib={'style': 'border-bottom: 1px solid black;'})
                td  = ET.SubElement(tr,'td')
                img = ET.SubElement(td,'img',attrib={'src': "%s_%s.png"%(edge_name,lig1)})
                td  = ET.SubElement(tr,'td')
                img = ET.SubElement(td,'img',attrib={'src': "%s_%s.png"%(edge_name,lig2)})
                
            fh = open( str(self.image_dir / "index.html"), "w" )
            fh.write( HTML.unescape( ET.tostring(html,encoding="unicode",
                                                 method='html') ) )
            fh.close()

        
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sys", type=str, required=True)
    parser.add_argument("--image_dir", type=str, default="results/imgdir")
    parser.add_argument("--single", type=str)
    parser.add_argument("--sub_dir", type=str, default="aq")
    parser.add_argument("--showidxs", action='store_true')
    parser.add_argument("--width", type=int, default=400)
    parser.add_argument("--height", type=int, default=-1)
    
    args = parser.parse_args()

    size = (args.width,args.height)

    if args.single is not None:
        ligand = Ligand(input=f"{args.sys}/unified/run/{args.single}/{args.sub_dir}/inputs/0.00000000_ti.mdin", 
                        parmfile=f"{args.sys}/unified/run/{args.single}/{args.sub_dir}/unisc.parm7", 
                        rstfile=f"{args.sys}/unified/run/{args.single}/{args.sub_dir}/stateA.rst7", 
                        image_loc="analysis/images", 
                        edge_name=f"{args.single}",
                        showidxs=args.showidxs,
                        size=size)
        ligand.run()
    else:
        Images = GenerateImages(system=args.sys, image_dir=args.image_dir, sub_dir=args.sub_dir,showidxs=args.showidxs,size=size)

import json

import numpy as np


from pathlib import Path

class ABFESetup:
    """ A class for setting up ABFE simulations using python. """
    def __init__(self):
        self.mdins = []

    def add_mdin(self, pymdin=None):
        if pymdin is None:
            raise ValueError("No PyMDin object provided.")
        else:
            self.mdins.append(pymdin)
            print(f"Added PyMDin object with {len(pymdin.blocks)} blocks.")

class PyAmber:
    """ A class for creating Amber input files.
    
    'PyAmber' is a class that can be used to create Amber input files. It is designed to be used in a Python script to generate input files for Amber simulations. The class is designed
    to be used in conjunction with the 'PyAmberBlock' class, which is used to create individual blocks in the input file.
    
    Attributes:
    -----------
    blocks : list
        A list of 'PyAmberBlock' objects that make up the input file.
    
    """
    def __init__(self):
        self.blocks = []

    def write(self, filename="equil_file.mdin"):
        """ Writes the PyAmber object to an Amber input file.
        
        Parameters:
        -----------
        filename : str
            The name of the file to which the input will be written. Default is "equil_file.mdin".
            
        Raises:
        -------
        ValueError
            If no blocks have been added to the PyAmber object.
        
        """
        if len(self.blocks) == 0:
            raise ValueError("No blocks have been added to the PyAmber object.")
        with open(filename, 'w') as f:
            for block in self.blocks:
                block_lines = block.write_lines()
                for line in block_lines:
                    f.write(line + "\n")
                f.write("\n")
                
    def add_block(self, block):
        """ Adds a block to the PyAmber object.
        
        Parameters:
        -----------
        block : PyAmberBlock
            The block to add to the PyAmber object.
        
        """
        self.blocks.append(block)
    
    def print_blocks(self):
        """ Prints the blocks in the PyAmber object.
        
        """
        print(f"PyAmber object with {len(self.blocks)} blocks.")
        for i, block in enumerate(self.blocks):
            print(f"Block {i}: {block.name}")
    
    def reorder_blocks(self, block_order):
        """ Reorders the blocks in the PyAmber object.
        
        Parameters:
        -----------
        block_order : list
            A list of integers indicating the new order of the blocks. For example, if block_order = [1, 0, 2], the second block in the PyAmber object will be moved to the first position, the first block will be moved to the second position, and the third block will remain in the third position.
        
        """
        new_blocks = []
        for i in block_order:
            new_blocks.append(self.blocks[i])
        self.blocks = new_blocks

class PyAmberBlock:
    def __init__(self, name="cntrl", default_cntrl=True):
        self.name = name
        self.data =  {}
        if self.name == "cntrl" and default_cntrl:
            self.set_default_cntrl()

        

    def write_lines(self):
        """ Writes the block as a list of lines which are returned.
        
        Returns:
            list: a list of strings, each of which is a line in the block.
        
        """
        primary_group = ["imin", "nstlim", "dt", "irest", "ntx", "ntxo", "ntc", "ntf", "ntwx", "ntpr", "cut", "iwrap"]
        ensemble_group = ["ntb", "ntp", "tempi", "temp0", "ntt", "gamma_ln", "tautp", "barostat", "ig"]

        lines = []
        lines.append(f"&{self.name}")
        flag=False
        for key, value in self.data.items():
            if key in primary_group:
                lines.append(f"{key:<20} = {value}")
                flag=True
        if flag:
            lines.append("\n")
        flag=False
        for key, value in self.data.items():
            if key in ensemble_group:
                lines.append(f"{key:<20} = {value}")
                flag=True
        if flag:
            lines.append("\n")
        flag=False
        for key, value in self.data.items():
            if key not in primary_group and key not in ensemble_group:
                lines.append(f"{key:<20} = {value}")

        lines.append("/")
        return lines
    
    def add_keypair(self, key, value):
        """ Adds a key-value pair to the block.
        
        Parameters:
        -----------
        key : str
            The key to add to the block.
        value : str
            The value to add to the block.
        
        """
        self.data[key] = value

    def print_keys(self):
        """ Prints the keys in the block.
        
        """
        primary_group = ["imin", "nstlim", "dt", "irest", "ntx", "ntxo", "ntc", "ntf", "ntwx", "ntpr", "cut", "iwrap"]
        ensemble_group = ["ntb", "ntp", "tempi", "temp0", "ntt", "gamma_ln", "tautp", "barostat", "ig"]
        for key in self.data.keys():
            if key in primary_group:
                print(f"{key} = {self.data[key]}")
        print(f"")
        for key in self.data.keys():
            if key in ensemble_group:
                print(f"{key} = {self.data[key]}")
        print(f"")
        for key in self.data.keys():
            if key not in primary_group and key not in ensemble_group:
                print(f"{key} = {self.data[key]}")

    def add_from_json(self, json_data):
        """ Adds key-value pairs to the block from a JSON block

        Parameters:
        -----------
        json_data : dict
            A dictionary containing key-value pairs to add to the block.

        """
        for key in json_data.keys():
            self.data[key] = json_data[key]

    def add_from_json_file(self, json_file):
        """ Adds key-value pairs to the block from a JSON file.

        Parameters:
        -----------
        json_file : str
            The name of the JSON file to read.

        """
        json_file = Path(json_file)
        if not json_file.exists():
            raise FileNotFoundError(f"File {json_file} not found.")
        with open(json_file, 'r') as f:
            json_data = json.load(f)
        self.add_from_json(json_data)

    def set_default_cntrl(self):
        """ Sets default data for the block.

        """
        default_data = {
            # Primary group
            "imin":     0,
            "nstlim":   1000,
            "dt":       0.001,
            "irest":    0,
            "ntx":      1,
            "ntxo":     1,
            "ntc":      2,
            "ntf":      1,
            "ntwx":     1000,
            "ntpr":     1000,
            "cut":      12.0,
            "iwrap":    1,
            # Ensemble group
            "ntb":      1,
            "ntp":      1,
            "tempi":    298.15,
            "temp0":    298.15,
            "ntt":      3,
            "gamma_ln": 1.0,
            "tautp":    2,
            "barostat": 2,
            # Other
            "ig":       -1,
        }
        for key in default_data.keys():
            self.data[key] = default_data[key]
        return
    
def MinRecipe():
    cntrl_recipe = {
        "imin": 1,
        "maxcyc": 5000,
        "ntmin": 2,
        "ntx": 1,
        "ntxo": 1,
        "ntpr": 100,
        "cut": 10,
        
        "ntr":1,
        "restraint_mask": "'!:WAT,Cl-,K+,Na+ & !@H='",
        "restraint_wt": 5.0,

        "ifsc": 1,
        "icfe": 1,
    }
    block = PyAmberBlock(name="cntrl")
    MDIN = PyAmber()
    block.add_from_json(cntrl_recipe)
    MDIN.add_block(block)
    return MDIN


if __name__ == "__main__":
    TestBlock = PyAmberBlock(name="cntrl")
    TestBlock.add_keypair("imin", 0)
    TestBlock.add_keypair("irest", 0)
    example_json = {"ntx": 1, "ntb": 1}
    TestBlock.add_from_json(example_json)
    TestBlock.add_from_json_file("test.json")

    Test2 = PyAmberBlock(name="wt")
    Test2.add_keypair("TYPE", "'TEMP0'")
    Test2.add_keypair("ISTEP1", 0)
    Test2.add_keypair("ISTEP2", 2000)
    Test2.add_keypair("VALUE1", 0.0)
    Test2.add_keypair("VALUE2", 298.15)

    Test3 = PyAmberBlock(name="wt")
    Test3.add_keypair("TYPE", "'TEMP0'")
    Test3.add_keypair("ISTEP1", 2000)
    Test3.add_keypair("ISTEP2", 4000)
    Test3.add_keypair("VALUE1", 298.15)
    Test3.add_keypair("VALUE2", 298.15)

    Test4 = PyAmberBlock(name="wt")
    Test4.add_keypair("TYPE", "'END'")



    MDIN = PyAmber()
    MDIN.add_block(TestBlock)
    MDIN.add_block(Test2)
    MDIN.add_block(Test3)
    MDIN.add_block(Test4)

    MDIN.write(filename="test.mdin")
    MDIN.print_blocks()
    MDIN.reorder_blocks([0, 2, 1, 3])
    MDIN.print_blocks()
    print("****ABFE****")

    NewABFE = ABFESetup()
    NewABFE.add_mdin(MDIN)

    min_recipe = MinRecipe()
    min_recipe.write(filename="min.mdin")
    NewABFE.add_mdin(min_recipe)
    print("****ABFE****")
    print(NewABFE.mdins)
    print("****ABFE****")

    

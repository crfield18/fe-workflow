#!/usr/bin/env python3
"""
This is a short script that writes a fewf map file for a given edge. This is useful in combination with mapinspect=1, which pauses before building the system. You can use this in combination with the fewf-image-writer code to generate images
of the system and work out problems with the soft core atoms. (Note - imagewriter requires you to move to the second stage of setup, which is less than ideal, but maybe we can fix that in the future.)
"""

def write_map_file(edgename, atomlist):
    with open(f"{edgename}.map.txt", 'w') as f:
        for atom in atomlist:
            f.write(f"  {atom} =>   {atom}\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Write a fewf map file")
    parser.add_argument("--edgename", required=True, help="The name of the edge.")
    parser.add_argument("--atomlist", nargs='+', help="The list of atoms that should be in the COMMON core.")
    args = parser.parse_args()

    write_map_file(args.edgename, args.atomlist)

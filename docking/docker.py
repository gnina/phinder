
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 00:35:41 2019
Program name: docker.py
@author: samcho
This literally docks the fragmet files into the receptor file.
This will be the first program you run when testing on something new. It only needs to be
run once. All the receptors have been docked already, so you should be find for a while.
The outputs ([fragment name]Docked.sdf) are used in the PharmacoFind set of programs.
NOTE: This hasn't been automated yet, and only works on one receptor at a time.
It shouldn't be a problem for you for a while, but in the
meantime, I'll be working on writing an automated version, and I'll send it to you when
it's ready. 
If you need an automated version, and I haven't sent it yet,
send me an email at Samthefluffysheep@gmail.com and I'll
write a script up within 3 days. :)
"""
from __future__ import print_function  # ensures print function compatibility with Python3
import argparse
# safely deal with file paths
import rdkit
import math
import os
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as scip
import numpy as np
import json
# library for loading the image
from rdkit import Chem, ForceField
import subprocess
import fileinput

# library for viewing the image

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------

"""
Puts it all together.
"""
def main():
    parser = argparse.ArgumentParser(description='docker.py takes a receptor-ligand pair and performs fragment docking with a given probe',formatter_class=argparse.RawDescriptionHelpFormatter)
    #Dock and include guanidine
    #Also, it looks like imidazole is the one that's trying to be positive.
    parser.add_argument("-r","--receptorFile",type=str, default="receptor.pdb",help="When docking a single Ligand-Receptor Pair, Receptor File. Takes .pdb, and might take other formats.")
    parser.add_argument("-l","--ligandFile",type=str, default="crystal_ligand.mol2",help="When docking a single Ligand-Receptor Pair, Ligand File. Takes .mol2, and might take other formats.")
    parser.add_argument("-p","--probeFiles",type=list, default=['benzene.sdf','isopropylamine.sdf', 'acetamide.sdf', 'isopropanol.sdf', 'imidazole.sdf', 'acetate.sdf', 'isobutane.sdf', 'water.sdf', 'guanidine.sdf'], help="List of probe files. Format: .sdf . (Default probes: benzene, isopropylamine, acetamide, isopropanol, imidazole, acetate, isobutane, water, guanidine)")
    parser.add_argument("-n","--numMode",type=int, default=100, help="Number of dockings to do for each probe. (Default: 100)")
    parser.add_argument("-w","--writeDirectory", type=str,help="When docking a single Ligand-Receptor Pair, Filename to write to. (Default: './PharmacoFindWorkplace/specifiedreceptor/')")
    parser.add_argument('-f', '--inputFile', type=str, help='Can read a csv file containing the locations for your ligands, receptors, and desired output directory in the column format. If expects that the file pairs the correct crystal ligand (.mol2 file) with the corresponding receptor (.pdb file). The format should look as follows (ligand, receptor, outputdirectory). All locations should be relative paths from the ./docker directory')
    working_dir = os.getcwd()
    args = parser.parse_args()
    if args.inputFile == None:
        probes = args.probeFiles
        recept = args.receptorFile
        out_dir = args.writeDirectory

        for probe in probes:     
            #For some reason, this bit doesn't do anything. 
            subprocess.run(["gnina", "-r", recept, "-l", probe, "-o", args.writeDirectory + "/" + probe[:-4]+"Docked.sdf", "--num_modes", str(args.numMode), "--autobox_ligand", args.ligandFile])
            print(probe)
        print(probes)
    else:
        file = np.genfromtxt(args.inputFile, dtype=str, delimiter=',')
        probes = args.probeFiles
        ligands, receptors, outs = file[:, 0], file[:, 1], file[:, 2]
        
        for probe in probes:
            for i,ligand in enumerate(receptors):
                subprocess.run(['gnina', "-r", str(receptors[i]) , "-l", str(probe), "-o", outs[i] + "/" + str(probe[:-4])+"Docked.sdf", "--num_modes", str(args.numMode), "--autobox_ligand", str(ligands[i])])
                print(probe, receptors[i], ligands[i], outs[i])
        print(probes)
if __name__ == '__main__':
    main()

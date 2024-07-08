#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
# IMPORTS
#-------------------------------------------------------------------------------

from __future__ import print_function  # ensures print function compatibility with Python3
import argparse
from sys import argv
# safely deal with file paths
import rdkit
import math
import os
#from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as scip
import numpy as np
from rdkit import Chem, ForceField
import subprocess
import fileinput
import json, sys
import scipy
from scipy import spatial
import io
import pandas as pd
import pickle

#-------------------------------------------------------------------------------
# FRAGMENTS
#-------------------------------------------------------------------------------

kinds=["Aromatic", "PositiveIon", "NegativeIon", "HydrogenDonor", "HydrogenAcceptor", "Hydrophobic"]

fragments = {
    'benzene': """ 
     RDKit          3D 
 
 12 12  0  0  0  0  0  0  0  0999 V2000 
   -1.1992    0.6647   -0.0436 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0117    1.3962   -0.0484 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.1860    0.7271   -0.0046 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.1942   -0.6703    0.0438 C   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0214   -1.3860    0.0482 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.1914   -0.7258    0.0045 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -2.1298    1.1810   -0.0775 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0041    2.4774   -0.0854 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1025    1.3011   -0.0086 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1297   -1.1818    0.0775 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0067   -2.4830    0.0858 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -2.1127   -1.3007    0.0084 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  2  0 
  2  3  1  0 
  3  4  2  0 
  4  5  1  0 
  5  6  2  0 
  6  1  1  0 
  1  7  1  0 
  2  8  1  0 
  3  9  1  0 
  4 10  1  0 
  5 11  1  0 
  6 12  1  0 
M  END 
""",
    "isopropylamine": """6363 
  -OEChem-07022016443D 
 
 13 12  0     0  0  0  0  0  0999 V2000 
   -0.1164    1.3888   -0.1236 N   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0019    0.0143    0.3611 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.3114   -0.6051   -0.1193 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.1969   -0.7981   -0.1181 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0012    0.0332    1.4566 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.4143   -1.6288    0.2562 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.3627   -0.6418   -1.2133 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1726   -0.0313    0.2406 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -2.1356   -0.3601    0.2390 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.1440   -1.8242    0.2614 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.2409   -0.8467   -1.2121 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.1088    1.3978   -1.1431 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.6993    1.9261    0.1680 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  1  0  0  0  0 
  1 12  1  0  0  0  0 
  1 13  1  0  0  0  0 
  2  3  1  0  0  0  0 
  2  4  1  0  0  0  0 
  2  5  1  0  0  0  0 
  3  6  1  0  0  0  0 
  3  7  1  0  0  0  0 
  3  8  1  0  0  0  0 
  4  9  1  0  0  0  0 
  4 10  1  0  0  0  0 
  4 11  1  0  0  0  0 
M  END 
> <PUBCHEM_COMPOUND_CID> 
6363 
 
> <PUBCHEM_CONFORMER_RMSD> 
0.4 
 
> <PUBCHEM_CONFORMER_DIVERSEORDER> 
1 
 
> <PUBCHEM_MMFF94_PARTIAL_CHARGES> 
4 
1 -0.99 
12 0.36 
13 0.36 
2 0.27 
 
> <PUBCHEM_EFFECTIVE_ROTOR_COUNT> 
0 
 
> <PUBCHEM_PHARMACOPHORE_FEATURES> 
3 
1 1 cation 
1 1 donor 
3 2 3 4 hydrophobe 
 
> <PUBCHEM_HEAVY_ATOM_COUNT> 
4 
 
> <PUBCHEM_ATOM_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ATOM_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ISOTOPIC_ATOM_COUNT> 
0 
 
> <PUBCHEM_COMPONENT_COUNT> 
1 
 
> <PUBCHEM_CACTVS_TAUTO_COUNT> 
1 
 
> <PUBCHEM_CONFORMER_ID> 
000018DB00000001 
 
> <PUBCHEM_MMFF94_ENERGY> 
-5.6468 
 
> <PUBCHEM_FEATURE_SELFOVERLAP> 
15.223 
 
> <PUBCHEM_SHAPE_FINGERPRINT> 
139733 1 8790891778173935434 
20096714 4 17183903008921445365 
21015797 1 9798534824276916602 
5943 1 11932023695004990610 
 
> <PUBCHEM_SHAPE_MULTIPOLES> 
77.34 
1.44 
1.24 
0.64 
0.04 
0.44 
0.02 
-0.55 
-0.09 
-0.03 
-0.09 
0 
-0.02 
0 
 
> <PUBCHEM_SHAPE_SELFOVERLAP> 
125.392 
 
> <PUBCHEM_SHAPE_VOLUME> 
55.2 
 
> <PUBCHEM_COORDINATE_TYPE> 
2 
5 
10 
 
$$$$""",
    "acetamide": """178 
  -OEChem-02222113093D 
 
  9  8  0     0  0  0  0  0  0999 V2000 
    0.6073   -1.1657    0.0018 O   0  0  0  0  0  0  0  0  0  0  0  0 
    0.7522    1.1306    0.0015 N   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.4275    0.0969    0.0013 C   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0680   -0.0618   -0.0046 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.8768   -0.6906   -0.6101 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.7917    0.0140    1.0288 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.7274    1.0655   -0.4084 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.7667    1.1492    0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.2815    2.0298    0.0081 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  4  2  0  0  0  0 
  2  4  1  0  0  0  0 
  2  8  1  0  0  0  0 
  2  9  1  0  0  0  0 
  3  4  1  0  0  0  0 
  3  5  1  0  0  0  0 
  3  6  1  0  0  0  0 
  3  7  1  0  0  0  0 
M  END 
> <PUBCHEM_COMPOUND_CID> 
178 
 
> <PUBCHEM_CONFORMER_RMSD> 
0.4 
 
> <PUBCHEM_CONFORMER_DIVERSEORDER> 
1 
 
> <PUBCHEM_MMFF94_PARTIAL_CHARGES> 
6 
1 -0.57 
2 -0.8 
3 0.06 
4 0.57 
8 0.37 
9 0.37 
 
> <PUBCHEM_EFFECTIVE_ROTOR_COUNT> 
0 
 
> <PUBCHEM_PHARMACOPHORE_FEATURES> 
2 
1 1 acceptor 
1 2 donor 
 
> <PUBCHEM_HEAVY_ATOM_COUNT> 
4 
 
> <PUBCHEM_ATOM_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ATOM_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ISOTOPIC_ATOM_COUNT> 
0 
 
> <PUBCHEM_COMPONENT_COUNT> 
1 
 
> <PUBCHEM_CACTVS_TAUTO_COUNT> 
2 
 
> <PUBCHEM_CONFORMER_ID> 
000000B200000001 
 
> <PUBCHEM_MMFF94_ENERGY> 
4.2976 
 
> <PUBCHEM_FEATURE_SELFOVERLAP> 
10.148 
 
> <PUBCHEM_SHAPE_FINGERPRINT> 
139733 1 9293842481710923235 
20096714 4 18268432510784477049 
21015797 1 8790889574902944993 
5943 1 8603650939448069255 
 
> <PUBCHEM_SHAPE_MULTIPOLES> 
71.47 
1.35 
1.13 
0.57 
0.57 
0.03 
0 
-0.1 
0 
-0.4 
0 
0.04 
0 
0 
 
> <PUBCHEM_SHAPE_SELFOVERLAP> 
122.421 
 
> <PUBCHEM_SHAPE_VOLUME> 
48.5 
 
> <PUBCHEM_COORDINATE_TYPE> 
2 
5 
10 
 
$$$$""",
    "isopropanol": """3776 
  -OEChem-07022016383D 
 
 12 11  0     0  0  0  0  0  0999 V2000 
   -0.0004    1.3572   -0.1242 O   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0000    0.0177    0.3601 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.2599   -0.6878   -0.1179 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.2603   -0.6871   -0.1179 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0001    0.0646    1.4540 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.3079   -1.7139    0.2590 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -2.1507   -0.1484    0.2213 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.3061   -0.7138   -1.2122 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.3089   -1.7132    0.2590 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1508   -0.1471    0.2213 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.3066   -0.7131   -1.2122 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0006    1.3242   -1.0961 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  1  0  0  0  0 
  1 12  1  0  0  0  0 
  2  3  1  0  0  0  0 
  2  4  1  0  0  0  0 
  2  5  1  0  0  0  0 
  3  6  1  0  0  0  0 
  3  7  1  0  0  0  0 
  3  8  1  0  0  0  0 
  4  9  1  0  0  0  0 
  4 10  1  0  0  0  0 
  4 11  1  0  0  0  0 
M  END 
> <PUBCHEM_COMPOUND_CID> 
3776 
 
> <PUBCHEM_CONFORMER_RMSD> 
0.4 
 
> <PUBCHEM_CONFORMER_DIVERSEORDER> 
1 
 
> <PUBCHEM_MMFF94_PARTIAL_CHARGES> 
3 
1 -0.68 
12 0.4 
2 0.28 
 
> <PUBCHEM_EFFECTIVE_ROTOR_COUNT> 
0 
 
> <PUBCHEM_PHARMACOPHORE_FEATURES> 
3 
1 1 acceptor 
1 1 donor 
3 2 3 4 hydrophobe 
 
> <PUBCHEM_HEAVY_ATOM_COUNT> 
4 
 
> <PUBCHEM_ATOM_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ATOM_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ISOTOPIC_ATOM_COUNT> 
0 
 
> <PUBCHEM_COMPONENT_COUNT> 
1 
 
> <PUBCHEM_CACTVS_TAUTO_COUNT> 
1 
 
> <PUBCHEM_CONFORMER_ID> 
00000EC000000001 
 
> <PUBCHEM_MMFF94_ENERGY> 
1.1698 
 
> <PUBCHEM_FEATURE_SELFOVERLAP> 
15.223 
 
> <PUBCHEM_SHAPE_FINGERPRINT> 
139733 1 8790885181088344531 
20096714 4 17400358344043462849 
21015797 1 9870305445774720250 
5943 1 13812955797099969682 
 
> <PUBCHEM_SHAPE_MULTIPOLES> 
76.45 
1.45 
1.19 
0.64 
0 
0.38 
0.02 
-0.54 
-0.09 
0 
-0.08 
0 
-0.03 
0 
 
> <PUBCHEM_SHAPE_SELFOVERLAP> 
124.806 
 
> <PUBCHEM_SHAPE_VOLUME> 
54.3 
 
> <PUBCHEM_COORDINATE_TYPE> 
2 
5 
10 
 
$$$$""",
    "imidazole": """795 
  -OEChem-07022016443D 
 
  9  9  0     0  0  0  0  0  0999 V2000 
    0.7451   -0.8866    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.1524    0.2675    0.0002 N   0  0  0  0  0  0  0  0  0  0  0  0 
    1.1083    0.4307    0.0001 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.6213   -0.9350   -0.0001 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0796    1.1234   -0.0002 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.3716   -1.6804    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1427    0.7424    0.0003 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.1627   -1.8702   -0.0002 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.2167    2.1949   -0.0003 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  3  1  0  0  0  0 
  1  4  1  0  0  0  0 
  1  6  1  0  0  0  0 
  2  4  2  0  0  0  0 
  2  5  1  0  0  0  0 
  3  5  2  0  0  0  0 
  3  7  1  0  0  0  0 
  4  8  1  0  0  0  0 
  5  9  1  0  0  0  0 
M  END 
> <PUBCHEM_COMPOUND_CID> 
795 
 
> <PUBCHEM_CONFORMER_RMSD> 
0.4 
 
> <PUBCHEM_CONFORMER_DIVERSEORDER> 
1 
 
> <PUBCHEM_MMFF94_PARTIAL_CHARGES> 
9 
1 0.03 
2 -0.57 
3 -0.3 
4 0.04 
5 0.08 
6 0.27 
7 0.15 
8 0.15 
9 0.15 
 
> <PUBCHEM_EFFECTIVE_ROTOR_COUNT> 
0 
 
> <PUBCHEM_PHARMACOPHORE_FEATURES> 
3 
1 1 donor 
3 1 2 4 cation 
5 1 2 3 4 5 rings 
 
> <PUBCHEM_HEAVY_ATOM_COUNT> 
5 
 
> <PUBCHEM_ATOM_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ATOM_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ISOTOPIC_ATOM_COUNT> 
0 
 
> <PUBCHEM_COMPONENT_COUNT> 
1 
 
> <PUBCHEM_CACTVS_TAUTO_COUNT> 
1 
 
> <PUBCHEM_CONFORMER_ID> 
0000031B00000001 
 
> <PUBCHEM_MMFF94_ENERGY> 
0.5083 
 
> <PUBCHEM_FEATURE_SELFOVERLAP> 
15.223 
 
> <PUBCHEM_SHAPE_FINGERPRINT> 
20096714 4 18338802341446370673 
21015797 1 9223230745331489536 
21040471 1 18410575076189396737 
5943 1 13961677281797930628 
 
> <PUBCHEM_SHAPE_MULTIPOLES> 
92.94 
1.3 
1.2 
0.58 
0 
0.02 
0 
-0.02 
0 
0.02 
0 
0.01 
0 
0 
 
> <PUBCHEM_SHAPE_SELFOVERLAP> 
189.355 
 
> <PUBCHEM_SHAPE_VOLUME> 
56.1 
 
> <PUBCHEM_COORDINATE_TYPE> 
2 
5 
10 
 
$$$$""",
    "acetate": """175 
  -OEChem-07022016433D 
 
  7  6  0     0  0  0  0  0  0999 V2000 
   -0.6906   -1.1195    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0 
   -0.6023    1.1668    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 
    1.4054   -0.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.1125    0.0025    0.0001 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.7934    0.4532   -0.8908 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.7688   -1.0820    0.0004 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.7937    0.4539    0.8902 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  4  1  0  0  0  0 
  2  4  2  0  0  0  0 
  3  4  1  0  0  0  0 
  3  5  1  0  0  0  0 
  3  6  1  0  0  0  0 
  3  7  1  0  0  0  0 
M  CHG  1   1  -1 
M  END 
> <PUBCHEM_COMPOUND_CID> 
175 
 
> <PUBCHEM_CONFORMER_RMSD> 
0.4 
 
> <PUBCHEM_CONFORMER_DIVERSEORDER> 
1 
 
> <PUBCHEM_MMFF94_PARTIAL_CHARGES> 
4 
1 -0.9 
2 -0.9 
3 -0.11 
4 0.91 
 
> <PUBCHEM_EFFECTIVE_ROTOR_COUNT> 
0 
 
> <PUBCHEM_PHARMACOPHORE_FEATURES> 
3 
1 1 acceptor 
1 2 acceptor 
3 1 2 4 anion 
 
> <PUBCHEM_HEAVY_ATOM_COUNT> 
4 
 
> <PUBCHEM_ATOM_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ATOM_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_DEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_BOND_UDEF_STEREO_COUNT> 
0 
 
> <PUBCHEM_ISOTOPIC_ATOM_COUNT> 
0 
 
> <PUBCHEM_COMPONENT_COUNT> 
1 
 
> <PUBCHEM_CACTVS_TAUTO_COUNT> 
1 
 
> <PUBCHEM_CONFORMER_ID> 
000000AF00000001 
 
> <PUBCHEM_MMFF94_ENERGY> 
1.0491 
 
> <PUBCHEM_FEATURE_SELFOVERLAP> 
15.277 
 
> <PUBCHEM_SHAPE_FINGERPRINT> 
139733 1 9222419305866417282 
20096714 4 18050288068225412184 
21015797 1 9727635016030451018 
5943 1 11571732774193885094 
 
> <PUBCHEM_SHAPE_MULTIPOLES> 
70.58 
1.31 
1.11 
0.56 
0.57 
0.01 
0 
-0.01 
0 
-0.37 
0 
0.04 
0 
0 
> <PUBCHEM_SHAPE_SELFOVERLAP> 
122.282 
 
> <PUBCHEM_SHAPE_VOLUME> 
47.5 
 
> <PUBCHEM_COORDINATE_TYPE> 
2 
5 
10 
 
$$$$""",
    "isobutane": """ 
     RDKit          3D 

 14 13  0  0  0  0  0  0  0  0999 V2000 
   -1.2723   -0.6848   -0.0019 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0032    0.0307    0.3220 C   0  0  0  0  0  0  0  0  0  0  0  0 
    1.2364   -0.7867    0.0340 C   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0604    1.3947   -0.3210 C   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.0665   -1.4688   -0.7509 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -1.7660   -1.0964    0.9229 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -2.0243    0.0455   -0.4084 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.0022    0.2094    1.4299 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.4066   -1.5153    0.8695 H   0  0  0  0  0  0  0  0  0  0  0  0 
    1.1969   -1.3222   -0.9144 H   0  0  0  0  0  0  0  0  0  0  0  0 
    2.1011   -0.0722    0.1098 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.8179    1.9658    0.0846 H   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.0133    1.3717   -1.4151 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.9597    1.9286    0.0390 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  1  0 
  2  3  1  0 
  2  4  1  0 
  1  5  1  0 
  1  6  1  0 
  1  7  1  0 
  2  8  1  0 
  3  9  1  0 
  3 10  1  0 
  3 11  1  0 
  4 12  1  0 
  4 13  1  0 
  4 14  1  0 
M  END""",
    "water": """ 
     RDKit          3D 
 
  3  2  0  0  0  0  0  0  0  0999 V2000 
   -0.0026    0.3674   -0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 
   -0.8225   -0.1846   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 
    0.8251   -0.1828    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  1  0 
  1  3  1  0 
M  END""",
    "guanidine": """3520
  -OEChem-02222113143D
  
  9  8  0     0  0  0  0  0  0999 V2000
   -0.7372   -1.1253    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6615    1.1693    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3504   -0.0414    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0483   -0.0026    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3272   -2.0530   -0.0003 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7495   -1.0635   -0.0009 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1948    2.0697    0.0004 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6757    1.1725    0.0009 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7069    0.9207    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  2  4  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  3  4  2  0  0  0  0
  3  9  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
3520
> <PUBCHEM_CONFORMER_RMSD>
0.4
> <PUBCHEM_MMFF94_ENERGY>
9.3667
$$$$"""
}


#-------------------------------------------------------------------------------
# FEATURE TYPES (SMILES)
#-------------------------------------------------------------------------------


AROMATIC= ["a1aaaaa1", "a1aaaa1"]
HYDROGEN_DONOR=["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
"[#8!H0&!$([OH][C,S,P]=O)]",
"[#16!H0]"]
HYDROGEN_ACCEPTOR=["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
"[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"]
POSITIVE_ION =["[+,+2,+3,+4]",
"[$(CC)](=N)N",
"[$([NH2]C(N)(=N))]", #originally [$(C(N)(N)=N)]]
"[$(n1cc[nH]c1)]"]
NEGATIVE_ION = ["[-,-2,-3,-4]",
"C(=O)[O-,OH,OX1]",
"[$([S,P](=O)[O-,OH,OX1])]",
"c1[nH1]nnn1",
"c1nn[nH1]n1",
"C(=O)N[OH1,O-,OX1]",
"C(=O)N[OH1,O-]",
"CO(=N[OH1,O-])",
"[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]"]
HYDROPHOBIC=[
    "a1aaaaa1",
"a1aaaa1",
"[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
"[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
"*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
"[C&r3]1~[C&r3]~[C&r3]1",
"[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
"[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
"[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
"[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
"[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
"[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
"[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
"[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
"[$([S]~[#6])&!$(S~[!#6])]",
    ]
    
#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------


def check_docked(dockedList, directory):
    dockedCount=len(dockedList) #total number of docked files that should be in the folder
    count=0
    dirFiles=os.listdir(directory) #turn the folder contents into a list 
    for x in dockedList: #check if each docked file is in the folder
        if x in dirFiles:
            count += 1
        else:
            count +=0
    if count>=dockedCount: #if each docked file is in the folder it will return true 
        return True
    else:
        return False

#Finds the center coordinates of each instance of pharmacophore by identifying each docked fragment that meets the requirement for a given pharmacophore type
#The "pharmacophores" generated by this function have not undergone any filtering or processing and are literally just one docked fragment
#Outputs: [[coordinates], minAffinity, cnnAffinity, cnnScore] 
def getCenters(fileName, smartList): #fileName is each docked file sdf (like benzenDocked.sdf), smartList (AKA SMART) is the type of pharmacophore like Aromatic 
    print(fileName)
    sdfFileData = rdkit.Chem.rdmolfiles.SDMolSupplier(fileName)
    #Coordlist is a list of coordinates (so a list of multiple [x,y,z] coordinates) of the center
    #of each present instance of the specified pharmacophore type
    centerList=[]
    energyList=[]
    cnnenergyList=[]
    cnnscoreList=[]
    for mol in sdfFileData:
        for smart in smartList:
            minAffinity=0
            patt = Chem.MolFromSmarts(smart)
            if(mol.HasSubstructMatch(patt)):
                #substructIndexList is a list of all instances that a given pharmacophore kind overlaps with mol (the docked fragment within the sdf)
                #substructIndexList=[substructIndex1,substructIndex2,...]
                #Each substructIndex is a list of atoms in each instance of overlap
                substructIndexList=mol.GetSubstructMatches(patt)
                for substructIndex in substructIndexList:
                    molConform=mol.GetConformers()
                    minAffinity=float(mol.GetProp('minimizedAffinity'))
                    cnnAffinity=float(mol.GetProp('CNNaffinity'))
                    cnnScore=float(mol.GetProp('CNNscore'))
                    
                    #It's at 0 because there is only one conformer for the molecule. 
                    #I turned it into conformer because there are functions that only work with
                    #Conformer objects, but this is still essentially the same thing as mol.
                    #Sort of like how Tomatoes are legally a vegetable, allowing them to go
                    #through customs as a vegetable, even if technically they are fruit.
                    atomPositions=[]
                    for atomIndex in substructIndex:
                        atomPositions.append(molConform[0].GetAtomPosition((atomIndex)))
                        #print(molConform[0].GetPositions())
                    x=0
                    y=0
                    z=0
                    for atomPosition in atomPositions:
                        x=x+atomPosition[0]
                        y=y+atomPosition[1]
                        z=z+atomPosition[2]
                        numAtoms=len(atomPositions)
                        
                    x=x/numAtoms
                    y=y/numAtoms
                    z=z/numAtoms
                    #coordList is a list of points, where each point represents an instsance of overlap between mol and feature type structure
                    #the actual coordinates are an average between the atoms in the overlapping structure
                    centerList.append([x,y,z])
                    energyList.append(minAffinity)
                    cnnenergyList.append(cnnAffinity)
                    cnnscoreList.append(cnnScore)
    return centerList,energyList,cnnenergyList,cnnscoreList

#Return a list of indices of each element in the cluster ordered by cluster number 
def cluster_indices(cluster_assignments):
    n = cluster_assignments.max() #n=number of clusters
    indices = []
    for cluster_number in range(1, n + 1):
        indices.append(np.where(cluster_assignments == cluster_number)[0].tolist())
    return indices

def calcClusterEnergy(indexList, energyList, cnnenergyList, cnnscoreList, fileList): 
    aveEnergyList=[] #average minimized affinity
    minEnergyList=[] #minimum min affinity 
    maxEnergyList=[] #maximimum min affinity 
    avecnnEnergyList=[] #average CNN affinity
    mincnnEnergyList=[] #minimum CNN affinity
    maxcnnEnergyList=[] #max CNN affinity
    avecnnScoreList=[] #average CNN score 
    mincnnScoreList=[] #minimum CNN score
    maxcnnScoreList=[] #maximum CNN score 
    fragtypeList=[] #list of fragment types that make up a cluster 
    clustersizeList=[] #list of cluster sizes
    for cluster in indexList:
        aveEnergy=0.0
        avecnnEnergy=0.0
        avecnnScore=0.0
        fragtypelistcluster=[]
        for item in cluster:
            aveEnergy+=energyList[item] 
            try:
                minEnergy
            except NameError:
                minEnergy=energyList[item]
            else:
                if energyList[item] < minEnergy:
                    minEnergy=energyList[item]
            try:
                maxEnergy
            except NameError:
                maxEnergy=energyList[item]
            else:
                if energyList[item] > maxEnergy:
                    maxEnergy=energyList[item]
            avecnnEnergy+=cnnenergyList[item]
            try:
                mincnnEnergy
            except NameError:
                mincnnEnergy=cnnenergyList[item]
            else:
                if cnnenergyList[item] < mincnnEnergy:
                    mincnnEnergy=cnnenergyList[item]
            try:
                maxcnnEnergy
            except NameError:
                maxcnnEnergy=cnnenergyList[item]
            else:
                if cnnenergyList[item] > maxcnnEnergy:
                    maxcnnEnergy=cnnenergyList[item]
            avecnnScore+=cnnscoreList[item]
            try:
                mincnnScore
            except NameError:
                mincnnScore=cnnscoreList[item]
            else:
                if cnnscoreList[item] < mincnnScore:
                    mincnnScore=cnnscoreList[item]
            try:
                maxcnnScore
            except NameError:
                maxcnnScore=cnnscoreList[item]
            else:
                if cnnscoreList[item] > maxcnnScore:
                    maxcnnScore=cnnscoreList[item]
            fragname=fileList[item] #clean up the fragname and remove "Docked.sdf"
            fragname=fragname[0:-10]
            fragtypelistcluster.append(fragname)     
        aveEnergy=aveEnergy/len(cluster) #average minimized affinity
        aveEnergyList.append(aveEnergy) 
        minEnergyList.append(minEnergy)
        maxEnergyList.append(maxEnergy)     
        avecnnEnergy=avecnnEnergy/len(cluster) #average CNN affinity 
        avecnnEnergyList.append(avecnnEnergy)
        mincnnEnergyList.append(mincnnEnergy)
        maxcnnEnergyList.append(maxcnnEnergy)
        avecnnScore=avecnnScore/len(cluster) #average CNN score 
        avecnnScoreList.append(avecnnScore)  
        mincnnScoreList.append(mincnnScore)
        maxcnnScoreList.append(maxcnnScore)
        fragtypeList.append(fragtypelistcluster) #fragment types 
        clustersizeList.append(len(cluster))
    return aveEnergyList, minEnergyList, maxEnergyList, avecnnEnergyList, mincnnEnergyList, maxcnnEnergyList, avecnnScoreList, mincnnScoreList, maxcnnScoreList, fragtypeList, clustersizeList

def unSplitter(li1, li2, li3, li4, li5, li6, li7, li8, li9, li10, li11, li12):
    comList=[]
    x=0
    while x<len(li2):
        comList.append([li1[x], li2[x], li3[x], li4[x], li5[x], li6[x], li7[x], li8[x], li9[x], li10[x], li11[x], li12[x]])
        x+=1
    return comList

#Sorts first by cluster size, then by increasing order of energy.
def sorter(combinedList):
    brick=sorted(combinedList, reverse=True, key=lambda cluster: (len(cluster[0]), -1*cluster[1]))
    return brick


#Converts the raw cluster into something readable, and sorts the results, and 
#Outputs a list [[[coordinatesOfCenterOfCluster], radius, avgerage minimized affinity, minimum min affinity, maximum min affinity, average cnn affinity, min cnn affinity, max cnn affinity, average cnn score, min cnn score, max cnn score, cluster size, guanidine count, water count, isobute count, ...]
def sortCenter(coordinates, cluster, energy, cnnEnergy, cnnScore, file):
    clustEnergy, clustEnergymin, clustEnergymax, clustcnnEnergy, clustcnnEnergymin, clustcnnEnergymax, clustcnnScore, clustcnnScoremin, clustcnnScoremax, fragList, sizeList=calcClusterEnergy(cluster_indices(cluster), energy, cnnEnergy, cnnScore, file)
    li=unSplitter(cluster_indices(cluster),clustEnergy, clustEnergymin, clustEnergymax, clustcnnEnergy, clustcnnEnergymin, clustcnnEnergymax, clustcnnScore, clustcnnScoremin, clustcnnScoremax, fragList, sizeList)
    sortedLi=sorter(li)
    returnList =[]
    for clust in sortedLi:
        guanidine=0
        water=0
        isobutane=0
        acetate=0
        imidazole=0
        isopropanol=0
        acetamide=0
        isopropylamine=0
        benzene=0
        clusterIndex=clust[0] #clusterIndex=a list of all of the indices in the cluster
        clusterCoordinates = [coordinates[i] for i in clusterIndex ]
        cent=clusterCenter(clusterCoordinates)
        radi=radius(cent, clusterCoordinates)
        typefragList=clust[10]
        for frag in typefragList:
            if frag == "guanidine":
                guanidine+=1
            if frag == "water":
                water+=1
            if frag == "isobutane":
                isobutane+=1 
            if frag == "acetate":
                acetate+=1
            if frag == "imidazole":
                imidazole+=1
            if frag == "isopropanol":
                isopropanol+=1
            if frag == "acetamide":
                acetamide+=1
            if frag == "isopropylamine":
                isopropylamine+=1
            if frag == "benzene":
                benzene+=1
        returnList.append([cent,radi, clust[1],clust[2],clust[3],clust[4],clust[5],clust[6],clust[7],clust[8],clust[9],len(clusterIndex),guanidine,water,isobutane,acetate,imidazole,isopropanol,acetamide,isopropylamine,benzene])
    return returnList

def clusterCenter(coordList):
    x=0.0
    y=0.0
    z=0.0
    for coordinate in coordList:
        x+=coordinate[0]
        y+=coordinate[1]
        z+=coordinate[2]
    nu=len(coordList)
    centerCoord=[x/nu,y/nu,z/nu]
    return centerCoord

def radius(center, coordinates): #use raidus of 1 A for cut off from center for feature labeling, >2 is false example, 1-2 is intermediate emit from training 
    #Distance between them all. Find the farthest distance, that's the radius.
    radiList=[]
    for coord in coordinates:
        radiList.append(math.sqrt((coord[0] - center[0])**2 + (coord[1] - center[1])**2 + (coord[2] - center[2])**2))
    radi=max(radiList)
    return radi

def jsoniload():
    jason = json.load(open('pharmit.json'))
    
    #[(print(pt['name'],pt['x'],pt['y'],pt['z'])) for pt in jason['points'] if pt['enabled']] 
    arra=[]
    [(arra.append([pt['name'],pt['x'],pt['y'],pt['z'], pt['radius']]))for pt in jason['points'] if pt['enabled']]
    #print(arra)

def jsonify(pharmaType,lia):
    
    li=[]
    for sphere in lia:
        x=sphere[0][0]
        y=sphere[0][1]
        z=sphere[0][2]
        radiu=sphere[1]
        enabled= True
        #sPharmit={"name": pharmaType, "x":x, "y":y, "z":z, "radius":radiu, "enabled":enabled}
        sPharmit={"name": pharmaType, "x":x, "y":y, "z":z, "radius":radiu, "enabled":enabled}
        li.append(sPharmit)
    return li

def allJsonify(typeName,thing):
    pointList=[]
    pointList.extend(jsonify(typeName,thing))
    pointDict={"points": pointList}
    #print(pointDict)
    return pointDict

#USE THIS BIT IF WANT TO FILTER OUT THE DATA
#Aromatic ->Higher threshold for number of instances
def purge(bigList, thresholds):
    purgedList=[]
    numPoints=thresholds.get("numPoint")
    energ=thresholds.get("energy")
    for ite in bigList:
        #If number of points is greater than numPoint and energy less than energ (set based on pharma kind)
        if(ite[11]>=numPoints and ite[2]<=energ):
            purgedList.append(ite)
    return purgedList

#returns a list [[pharmakind, [x,y,z], radius, avg. min affinity, min min affinity, max min affinity, avg cnn affinity, min cnn affinity, max cnn affinity, avg cnn score, min cnn score, max cnn score, cluster size, fragment counts for each type],...]
def all_kinds_list(spheres_list,kind):
    li=[]
    for sphere in spheres_list:
        li2=[]
        li2.append(kind)
        for item in sphere:
            li2.append(item)
        li.append(li2)
    return li

#returns a dataframe of all spheres of all types 
def make_df(final_spheres):
    li0,li1,li2,li3,li4,li5,li6,li7,li8,li9,li10,li11,li12,li13,li14,li15,li16,li17,li18,li19,li20,li21,li22,li23=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for sphere in final_spheres:
        li0.append(sphere[0]) #kind
        li1.append(sphere[1][0]) #x coord
        li2.append(sphere[1][1]) #y coord
        li3.append(sphere[1][2]) #z coord
        li4.append(sphere[2]) #radius
        li5.append(sphere[3]) #ave min affinity
        li6.append(sphere[4]) #min min affinity
        li7.append(sphere[5]) #max min affinity
        li8.append(sphere[6]) #ave cnn affinity
        li9.append(sphere[7]) #min cnn affinity
        li10.append(sphere[8]) #max cnn affinity
        li11.append(sphere[9]) #ave cnn score
        li12.append(sphere[10]) #min cnn score
        li13.append(sphere[11]) #max cnn score 
        li14.append(sphere[12]) #size
        li15.append(sphere[13]) #guanidine
        li16.append(sphere[14]) #water
        li17.append(sphere[15]) #isobutane
        li18.append(sphere[16]) #acetate
        li19.append(sphere[17]) #imidazole
        li20.append(sphere[18]) #isopropanol
        li21.append(sphere[19]) #acetamide 
        li22.append(sphere[20]) #isopropylamine 
        li23.append(sphere[21]) #benzene
    zipped=list(zip(li1,li2,li3,li0,li4,li5,li6,li7,li8,li9,li10,li11,li12,li13,li14,li15,li16,li17,li18,li19,li20,li21,li22,li23))
    df=pd.DataFrame(zipped,columns=['x','y','z','pharmaKind','radius','ave minAffinity','min minAffinity','max minAffinity','ave cnnAffinity','min cnnAffinity',
        'max cnnAffinity','ave cnnScore','min cnnScore','max cnnScore','size','guanidineCount','waterCount','isobutaneCount','acetateCount','imidazoleCount',
        'isopropanolCount','acetamideCount','isopropylamineCount','benzenCount'])
    return df

#rank features with machine learning models and outputs in dataframe
models=["LinearRegression","LogisticRegression","NeuralNetwork"]
def rank_ml(features): #input is a dataframe of all features of all kinds 
    all_kinds_ranked=[]
    for kind in kinds:
        df=features[features.pharmaKind==kind].copy() #sort for just a given pharmaKind 
        D=df.drop(["x","y","z","pharmaKind"],axis=1).to_numpy() #remove unncessary columns for ranking 
        for m in models:
            filename="/net/pulsar/home/koes/ron33/pdb-bind-refined/"+m+"Final_"+kind+".sav"
            model=pickle.load(open(filename,'rb')) #load pre-trained model
            if hasattr(model,"predict_proba"):
                p=model.predict_proba(D)[:,1]
            else:
                p=model.predict(D)
            df[m]=p #create a new column for predictions for each model type 
        all_kinds_ranked.append(df)
    all_kinds_ranked_df=pd.concat(all_kinds_ranked)
    return all_kinds_ranked_df

def rank_ml_thresholded(df,n,model): #n is the number of features to be returned in the final json
    top1=df[model].max() #top prediction with given model
    index=df[df[model]==top1].index #index of top prediction
    T=df[df[model]==top1] #T will become the dataframe of the top x 
    df=df.drop(index=index, axis=1) #drop the top row so we can select from the others
    while len(T)<n:
        topn=df[model].max()
        row=df[df[model]==topn] #row with the top prediction
        pk=row.pharmaKind.item() #pharmaKind, index, x, y, and z info from the row
        index=row.index
        x1=row.x.item()
        y1=row.y.item()
        z1=row.z.item()
        if len(T[T.pharmaKind==pk])==0: #no feature of this pharmaKind yet in the top 10
            T=pd.concat([T,row],ignore_index=True) #add to list of top 10
            df=df.drop(index=index,axis=0) #remove from overall dataframe
        else:
            df2=T[T.pharmaKind==pk] #dataframe of all of the existing features of this type in the topn 
            overlap=0 
            for r in range(len(df2)):
                x2=df2.iloc[r].x.item()
                y2=df2.iloc[r].y.item()
                z2=df2.iloc[r].z.item()
                dist=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                if dist<=2:
                    overlap+=1
            if overlap!=0: #this indicates that there is one or more feature of the same pharmakind within less than a 2A distnace already in the top 10
                df=df.drop(index=index,axis=0) #remove from overall dataframe and don't add to top 10  
            else:
                T=pd.concat([T,row],ignore_index=True) #add to list of top N
                df=df.drop(index=index,axis=0) #remove from overall dataframe
    return(T) #T is a dataframe of the Top N with no features of the same kind within 2 A
#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------


"""
Puts it all together.
"""
def main():

    fileList = ['benzeneDocked.sdf','isopropylamineDocked.sdf', 'acetamideDocked.sdf', 'isopropanolDocked.sdf', 'imidazoleDocked.sdf', 'acetateDocked.sdf', 'isobutaneDocked.sdf', 'waterDocked.sdf', 'guanidineDocked.sdf']
    
##

    parser = argparse.ArgumentParser(description='Phinder generates pharmacophore features from receptor structures using fragment docking via GNINA',formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    ### Receptor and Ligand Parameters
    parser.add_argument("-r","--receptorFile",type=str, default="receptor.pdb",help="When docking a single Ligand-Receptor Pair, Receptor File. Takes .pdb, and might take other formats.")
    parser.add_argument("-p", "--pocket", choices=["ligand","coords"], default="ligand", help='Choices are "ligand" and "coords." If pocket="ligand", the user must provide a ligand file by which the pocket will be defined. If pocket="coords", the user must provide coordinates.')
    parser.add_argument("-l","--ligandFile",type=str, required=("ligand" in argv),default="crystal_ligand.mol2",help="When docking a single Ligand-Receptor Pair, Ligand File. Takes .mol2, and might take other formats.")
    parser.add_argument("--cords",required=("coords" in argv), type=list, help="Input a list of coordinates formated as [X coordinate of the center, Y coordinate of the center, Z coordinate of the center, size in the X dimension, size in the Y dimension, size in the Z dimension]. Sizes are in Angstroms")
    
    ### Docking Parameter
    parser.add_argument("-n","--numMode",type=int, default=100, help="Number of dockings to do for each probe. (Default: 100)")
    
    ### Directory to Write
    parser.add_argument("-w","--writeDirectory", type=str, default=os.getcwd(),help="When docking a single Ligand-Receptor Pair, Filename to write to. (Default is the Current Directory)")

    ### Machine Learning Ranking of Features 
    parser.add_argument("-m", "--machineLearning",type=bool, default=False, help="If True, machine learning ranking of features will be performed") #specificies if machine learning-based ranking of features should be performed 


    args = parser.parse_args()
    fragmentPrefilter = [-2, -2, -2, -2, -2, -2, -2, -1, -2]

    recept = args.receptorFile
    out_dir = args.writeDirectory



#-------------------------------------------------------------------------------
# RUN GNINA
#-------------------------------------------------------------------------------

    if check_docked(fileList,out_dir)==False: #run gnina only if the docked sdf files are not already in the output folder 
        for key in fragments:
            f = open(out_dir+"/tempFragment.sdf", "w")
            f.write(fragments[key])
            f.close()
            if args.pocket=="ligand": #use autobox to identify the pocket
                subprocess.run(["gnina", "-r", recept, "-l", out_dir+"/tempFragment.sdf", "-o", out_dir + "/" + key+"Docked.sdf", "--num_modes", str(args.numMode), "--autobox_ligand", args.ligandFile])
            if args.pocket=="coords": #user must manually define the pocket
                subprocess.run(["gnina", "-r", recept, "-l", out_dir+"/tempFragment.sdf", "-o", out_dir + "/" + key+"Docked.sdf", "--num_modes", str(args.numMode), "--center_x", args.pocket[0], "--center_y", args.pocket[1], 
                    "--center_z", args.pocket[2], "--size_x", args.pocket[3], "--size_y", args.pocket[4], "--size_z", args.pocket[5]])
            os.remove(out_dir+'/tempFragment.sdf')
    else:
        pass 
#-------------------------------------------------------------------------------
# Generate Pharmacophore
#-------------------------------------------------------------------------------
    all_kinds=[]
    receptor_name=recept #parse out pdb id of receptor
    if "pdb" in receptor_name:
        receptor_name=receptor_name.replace(".pdb","")
        if "/" in receptor_name:
            receptor_name=receptor_name.split("/")
            receptor_name=receptor_name[-1]
    elif "mol2" in receptor_name:
        receptor_name=receptor_name.replace(".mol2","")
        if "/" in receptor_name:
            receptor_name=receptor_name.split("/")
            receptor_name=receptor_name[-1]
    else:
        raise Exception("Receptor file has unknwon format. Receptor name can only be extracted from .pdb and .mol2 files")   

    for kind in kinds:
        print(kind)
        if(kind=="Hydrophobic"):
            SMART=HYDROPHOBIC
            nuPo= 10
            maxEne= -2
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
        elif(kind=="Aromatic"):
            SMART= AROMATIC
            nuPo= 10
            maxEne= -2
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
        elif(kind=="HydrogenDonor"):
            SMART=HYDROGEN_DONOR
            nuPo= 10
            maxEne= -2
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
        elif(kind=="HydrogenAcceptor"):
            SMART=HYDROGEN_ACCEPTOR
            nuPo= 10
            maxEne= -2
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
        elif(kind=="NegativeIon"):
            SMART=NEGATIVE_ION
            nuPo= 10
            maxEne= -2
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
        elif(kind=="PositiveIon"):
            SMART=POSITIVE_ION
            nuPo=6
            maxEne= -2
            clusterThreshold = 1
            THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}

        allSpheres=[]
        thingCenters, thingEnergy, thingcnnEnergy, thingcnnScore, thingFile=[],[],[],[],[] #thingFile is a list of the fragments were each center came fromhbg nb
        filteredCenters, filteredEnergy, filteredcnnEnergy, filteredcnnScore, filteredFile = [],[],[],[],[]

        for aia in fileList: #fileList is the list of docked fragment sdf files 
            thingCentersTem, thingEnergyTem, thingcnnEnergyTem, thingcnnScoreTem=getCenters(os.path.join(out_dir,aia), SMART)
            thingCenters.extend(thingCentersTem)
            thingEnergy.extend(thingEnergyTem)
            thingcnnEnergy.extend(thingcnnEnergyTem)
            thingcnnScore.extend(thingcnnScoreTem)
            for i, k in enumerate(thingEnergyTem):
                thingFile.append(aia)

        print(thingCenters)
            
        #below is used to prefilter fragments
        #thingCluster=[]
        #if((len(thingCenters)>0)):
            # Defining a prefiltering of Affinities by Feature Type
            #things = [thingCenters, thingEnergy, thingcnnEnergy, thingcnnScore, thingFile] #is this necessary? 
            # Docked fragment prefiltering 
            # Default Criteria: Distance; Default Method: Complete (As defined in the args section above)
            #fragmentPrefilt = fragmentPrefilter
            #for i, thing in enumerate(thingCenters):
                #fileIndex = fileList.index(thingFile[i])
                #fragPrefilter = fragmentPrefilt[fileIndex]
                #if thingEnergy[i] < fragPrefilter:
                    #filteredCenters.append(thingCenters.pop(i))
                    #filteredEnergy.append(thingEnergy.pop(i))
                    #filteredcnnEnergy.append(thingcnnEnergy.pop(i))
                    #filteredcnnScore.append(thingcnnScore.pop(i))
                    #filteredFile.append(thingFile.pop(i)) #need to do this step to match fragment with cluster 

            #filteredCenters = np.array(filteredCenters)
            #filteredEnergy = np.array(filteredEnergy)
            #filteredcnnEnergy = np.array(filteredcnnEnergy)
            #filteredcnnScore = np.array(filteredcnnScore)
            #filteredFile = np.array(filteredFile)
                
        #cluster features 
        #hierarchechal clustering settings
        #settings are currently the same across all features types  
        criteria="distance"
        method="complete"
        clusterThreshold=1

        if len(thingCenters) != 0:
            thingCluster=scip.fclusterdata(thingCenters, t = clusterThreshold, criterion = criteria, method = method)
            #thingCluster is a list of the same length as thingCenters were each number is the assigned cluster for each item in thingCenters
            #for example thingCluster=[5,4,...], then the first center is in cluster 5, the second one is in cluster 4, etc
            thingCenterspheres=sortCenter(thingCenters, thingCluster, thingEnergy, thingcnnEnergy, thingcnnScore, thingFile) #if had done fragment prefiltering, would have to change these inputs
            allSpheres.extend(thingCenterspheres)
            #filter clusters based on energy and cluster size cutoffs for each pharma kind
            #allSpheres=purge(allSpheres, THRESHOLD)


        allSphereskind=all_kinds_list(allSpheres,kind) 
        for item in allSphereskind:
            all_kinds.append(item)

    features_df=make_df(all_kinds) #dataframe of all features

    if args.machineLearning==False:
        if os.path.exists(out_dir+'/'+receptor_name+'.csv') == False: #make a csv of all features only if it is not already present 
            features_df.to_csv(out_dir+'/'+receptor_name+'.csv', index=False)
        if os.path.exists(out_dir+'/'+'GeneratedPharma.json') == False:
            json_list=[]
            for index, row in features_df.iterrows():
                enabled= True
                s={'name': row["pharmaKind"],'x': row["x"],'y':row["y"],'z':row["z"],'radius':1,'enabled':enabled} #I set radius for all json features to 1. will still retain ML feature of varying radius
                json_list.append(s)  
            pointDict={"points": json_list}
            with open(out_dir+'/GeneratedPharma.json','w') as f:
                jason=json.dumps(pointDict,indent=4)
                f.write(jason)
                f.close()
        else:
            pass

    if args.machineLearning==True:
        features_df=rank_ml(features_df)
        features_df=features_df.sort_values(by=["NeuralNetwork"],ascending=False)
        if os.path.exists(out_dir+'/'+receptor_name+'.csv') == False: #make a csv of all features only if it is not already present 
            features_df.to_csv(out_dir+'/'+receptor_name+'.csv', index=False)
        features_df.to_csv(out_dir+'/'+receptor_name+'.csv', index=False)
        #this would make a json of the Top N features
        topN=10
        top_features_df=rank_ml_thresholded(features_df, topN, "NeuralNetwork")
        #this makes a json of all features
        if os.path.exists(out_dir+'/GeneratedPharma.json') == False:
            json_list=[]
            for index, row in features_df.iterrows():
                enabled= True
                #if machine learning is set to true, features will be ranked based on neural network scores (although this could be changed) and it will also print out the score 
                s={'name': row["pharmaKind"],'x': row["x"],'y':row["y"],'z':row["z"],'radius':1,'enabled':enabled, "rank":row["NeuralNetwork"]} #I set radius for all json features to 1. will still retain ML feature of varying radius
                json_list.append(s)  
            pointDict={"points": json_list}
            with open(out_dir+'/GeneratedPharma.json','w') as f:
                jason=json.dumps(pointDict,indent=4)
                f.write(jason)
                f.close()
        if os.path.exists(out_dir+'/Top'+str(topN)+'.json') == False:
            json_list=[]
            for index, row in top_features_df.iterrows():
                enabled= True
                #if machine learning is set to true, features will be ranked based on neural network scores (although this could be changed) and it will also print out the score 
                s={'name': row["pharmaKind"],'x': row["x"],'y':row["y"],'z':row["z"],'radius':1,'enabled':enabled, "rank":row["NeuralNetwork"]} #I set radius for all json features to 1. will still retain ML feature of varying radius
                json_list.append(s)  
            pointDict={"points": json_list}
            with open(out_dir+'/Top'+str(topN)+'.json','w') as f:
                jason=json.dumps(pointDict,indent=4)
                f.write(jason)
                f.close()

if __name__ == '__main__':
    main()

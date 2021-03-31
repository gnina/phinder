#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
# IMPORTS
#-------------------------------------------------------------------------------

from __future__ import print_function  # ensures print function compatibility with Python3
import argparse
# safely deal with file paths
import rdkit
import math
import os
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as scip
import numpy as np
from rdkit import Chem, ForceField
import subprocess
import fileinput
import json, sys
import scipy
from scipy import spatial
import io
#-------------------------------------------------------------------------------
# FRAGMENTS
#-------------------------------------------------------------------------------

dockedFiles = ["benzeneDocked.sdf", "acetateDocked.sdf", "acetamideDocked.sdf", "isobutaneDocked.sdf", "isopropanolDocked.sdf", "isopropylamine.sdf", "imidazoleDocked.sdf", "guanadineDocked.sdf", "waterDocked.sdf"]

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
"[$(C(N)(N)=N)]",
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
pharmakinds = ["HydrogenDonor","HydrogenAcceptor","PositiveIon","NegativeIon","Hydrophobic","Aromatic"]

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------


def check_docked(dockedList, directory):
    for x in dockedList:
        allItems = len(dockedList)
        count = 0
        if any(directory) == x:
            count += 1
        else:
            count +=0
    if allItems >= count:
        return True
    else:
        return False

def compute_matches(ref, test, kind):
    '''Return number of matches and not matched points for given kind'''
    refxyz = [(pt['x'],pt['y'],pt['z']) for pt in ref if pt['name'] == kind]
    testxyz = [(pt['x'],pt['y'],pt['z']) for pt in test if pt['name'] == kind]
    testradii = [pt['radius'] for pt in test if pt['name'] == kind]
    
    if len(refxyz) == 0 or len(testxyz) == 0:
        return 0,len(refxyz)+len(testxyz)
    dists = spatial.distance.cdist(refxyz,testxyz) #rows are ref, columns test
    matches = 0

    while np.any(np.isfinite(dists)):
        #find shortest distance
        r,t = np.unravel_index(np.argmin(dists), dists.shape)  #argmin returns flattened index
        #r is the index of the ref point, t of the test point
        if dists[r,t] < testradii[t]: #a match!
            matches += 1
            dists[r,:] = np.inf #mark points as matched - can only match once
            dists[:,t] = np.inf
        else:
            dists[:,t] = np.inf #there's nothing this can match with
        
        test_len = len(testxyz)
        ref_len = len(refxyz)
        intersection = ref_len + test_len - 2*matches
    return matches, intersection, ref_len, test_len


            
def print_tani(m,u,kind):
    '''print tanimoto'''
    if m+u > 0:
        tani = m/float(m+u)
    else:
        tani = 0;
    
    #print('%16s Tanimoto: %.3f\tMatches: %d\tUnmatches: %d'%(kind,tani,m,u))
    return tani

def print_recall(matches, relevant_points):
    '''prints the recall'''
    if matches+relevant_points > 0:
        recall = matches/relevant_points
    else:
        recall = 0

    return recall

def print_precision(matches, all_positives):
    '''prints precision'''
    if matches+all_positives > 0:
        precision = matches/all_positives
    else:
        precision = 0

    return precision

def splitter(combinedList):
    centerList=[]
    energyList=[]
    for item in combinedList:
        centerList.append(item[0])
        energyList.append(item[1])
    return centerList, energyList

#Finds the center coordinates of each instance of pharmacophore
#Outputs: [[coordinates], energy]
def getCenters(fileName, smartList):
    sdfFileData = rdkit.Chem.rdmolfiles.SDMolSupplier(fileName)
    #Coordlist is a list of coordinates (so a list of multiple [x,y,z] coordinates) of the center
    #of each present instance of the specified pharmacophore type
    coordList=[]
    for mol in sdfFileData:
        minAffinity=0
        for smart in smartList:
            patt = Chem.MolFromSmarts(smart)
            if(mol.HasSubstructMatch(patt)):
                
                #substructIndexList is a list of lists of atom indexes. 
                #Each internal list contains the atom index values of the atoms that make up a
                #pharmacophore
                substructIndexList=mol.GetSubstructMatches(patt)
                #print(substructIndexList)
                #substructIndexList -> [substructIndex1, substructIndex2, ...]
                #substructIndex -> [atomIndex1, atomIndex2, atomIndex3, ...]
                #so: substructIndexList -> [[int1, int2, int3, ...], [int1, int2, int5, ...], ...]
                
                for substructIndex in substructIndexList:
                    molConform=mol.GetConformers()
                    minAffinity=float(mol.GetProp('minimizedAffinity'))
                    
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
                   # print([x,y,z])
                    coordList.append([[x,y,z], minAffinity])
    
    return coordList

#Obtains the indices of each element in each cluster.
def cluster_indices(cluster_assignments):
    n = cluster_assignments.max()
    indices = []
    for cluster_number in range(1, n + 1):
        indices.append(np.where(cluster_assignments == cluster_number)[0].tolist())
    return indices

def calcClusterEnergy(indexList, energyList):
    aveEnergyList=[]
    for cluster in indexList:
        aveEnergy=0.0
        for item in cluster:
            aveEnergy+=energyList[item]
        aveEnergy=aveEnergy/len(cluster)
        aveEnergyList.append(aveEnergy)
    return aveEnergyList

def unSplitter(li1, li2):
    comList=[]
    x=0
    while x<len(li2):
        comList.append([li1[x], li2[x]])
        x+=1
    return comList

#Sorts first by cluster size, then by increasing order of energy.
def sorter(combinedList):
    brick=sorted(combinedList, reverse=True, key=lambda cluster: (len(cluster[0]), -1*cluster[1]))
    return brick


#Converts the raw cluster into something readable, and sorts the results, and 
#Outputs a list [[[coordinatesOfCenterOfCluster], radius,[Average Energy]],[[coordinatesOfCenterOfCluster], radius,[Average Energy]]...]
def sortCenter(coordinates, typeCluster, typeEnergy):
    clustEnergy=calcClusterEnergy(cluster_indices(typeCluster), typeEnergy)
    li=unSplitter(cluster_indices(typeCluster),clustEnergy)
    sortedLi=sorter(li)
    returnList =[]
    for clust in sortedLi:
        clusterIndex=clust[0]
        #print(clusterIndex)
        clusterCoordinates = [coordinates[i] for i in clusterIndex ]
        cent=clusterCenter(clusterCoordinates)
        radi=radius(cent, clusterCoordinates)
        returnList.append([cent, radi, clust[1], len(clusterIndex)])
        #print(len(clusterIndex))
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

def radius(center, coordinates):
    #Distance between them all. Find the farthest distance, that's the radius.
    radi=1.0
    #for coord in coordinates:
       # dist=math.sqrt((coord[0] - center[0])**2 + (coord[1] - center[1])**2 + (coord[2] - center[2])**2)
        #if dist>=radi:
         #   radi=dist
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
        #If number of points is greater than 10 and energy less than -3.0
        if(ite[3]>=numPoints and ite[2]<=energ):
            purgedList.append(ite)
    return purgedList



#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------


"""
Puts it all together.
"""
def main():
    
##

    parser = argparse.ArgumentParser(description='Phinder generates pharmacophore features from receptor structures using fragment docking via GNINA',formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    ### Using a File input of Locations for Pharmacophore Generation
    parser.add_argument('-f', '--inputFile', type=str, help='Can read a csv file containing the locations for your ligands, receptors, and desired output directory in the column format. If expects that the file pairs the correct crystal ligand (.mol2 file) with the corresponding receptor (.pdb file). The format should look as follows (ligand, receptor, outputdirectory). All locations should be relative paths from the ./docker directory')
    
    
    ### Receptor and Ligand Parameters
    parser.add_argument("-r","--receptorFile",type=str, default="receptor.pdb",help="When docking a single Ligand-Receptor Pair, Receptor File. Takes .pdb, and might take other formats.")
    parser.add_argument("-l","--ligandFile",type=str, default="crystal_ligand.mol2",help="When docking a single Ligand-Receptor Pair, Ligand File. Takes .mol2, and might take other formats.")
    
    ### Docking Parameter
    parser.add_argument("-n","--numMode",type=int, default=100, help="Number of dockings to do for each probe. (Default: 100)")
    
    ### Pharmacophore Feature Generation Parameters
    
    ### Directory to Write
    parser.add_argument("-w","--writeDirectory", default='.', type=str,help="When docking a single Ligand-Receptor Pair, Filename to write to. (Default is the Current Directory)")
    
    working_dir = os.getcwd()
    args = parser.parse_args()
    fragmentPrefilter = [-3.49, -3.19, -2.66, -2.33, -3.00, -2.50, -3.00, -1.00, -3.00]
    kinds=["Aromatic", "PositiveIon", "NegativeIon", "HydrogenDonor", "HydrogenAcceptor", "Hydrophobic"]
    
    if args.inputFile == None:
        recept = args.receptorFile
        out_dir = args.writeDirectory

#-------------------------------------------------------------------------------
# RUN GNINA
#-------------------------------------------------------------------------------
        
        for key in fragments:
            f = open("tempFragment.sdf", "w")
            f.write(fragments[key])
            f.close()
            subprocess.run(["gnina", "-r", recept, "-l", "tempFragment.sdf", "-o", out_dir + "/" + key+"Docked.sdf", "--num_modes", str(args.numMode), "--autobox_ligand", args.ligandFile])
            os.remove('tempFragment.sdf')
#-------------------------------------------------------------------------------
# Generate Pharmacophore
#-------------------------------------------------------------------------------
        for kind in kinds:
            if(kind=="Hydrophobic"):
                SMART=HYDROPHOBIC
                nuPo= 10
                maxEne= -2
                criteria = "distance"
                method = "ward"
                clusterThreshold = 5
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
            elif(kind=="Aromatic"):
                SMART= AROMATIC
                nuPo= 10
                maxEne= -2
                criteria = 'distance'
                method = 'ward'
                clusterThreshold = 9
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
            elif(kind=="HydrogenDonor"):
                SMART=HYDROGEN_DONOR
                nuPo= 10
                maxEne= -2
                criteria = 'distance'
                method = 'ward'
                clusterThreshold = 9
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
            elif(kind=="HydrogenAcceptor"):
                SMART=HYDROGEN_ACCEPTOR
                nuPo= 10
                maxEne= -2
                criteria = 'distance'
                method = 'ward'
                clusterThreshold = 3
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
            elif(kind=="NegativeIon"):
                SMART=NEGATIVE_ION
                nuPo= 10
                maxEne= -2
                criteria = 'distance'
                method = 'ward'
                clusterThreshold = 1
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
            elif(kind=="PositiveIon"):
                SMART=POSITIVE_ION
                nuPo=6
                maxEne= -2
                criteria = 'distance'
                method = 'complete'
                clusterThreshold = 9
                THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}

            allSpheres=[]
            thingCenters, thingEnergy=[],[]
            thingFile = []
            filteredCenters, filteredEnergy = [],[]
            fileList = ['benzeneDocked.sdf','isopropylamineDocked.sdf', 'acetamideDocked.sdf', 'isopropanolDocked.sdf', 'imidazoleDocked.sdf', 'acetateDocked.sdf', 'isobutaneDocked.sdf', 'waterDocked.sdf', 'guanidineDocked.sdf']

            for aia in fileList:
                thingCentersTem, thingEnergyTem=splitter(getCenters(os.path.join(out_dir,aia), SMART))
                thingCenters.extend(thingCentersTem)
                thingEnergy.extend(thingEnergyTem)
                for i, k in enumerate(thingEnergyTem):
                    thingFile.append(aia)
                #Cluster info. This isn't really useful without processing it using 

            thingCluster=[]
            if((len(thingCenters)>0)):
                # Defining a prefiltering of Affinities by Feature Type
                things = [thingCenters, thingEnergy, thingFile]
                fragmentPrefilt = fragmentPrefilter
                for i, thing in enumerate(thingCenters):
                    fileIndex = fileList.index(thingFile[i])
                    fragPrefilter = fragmentPrefilt[fileIndex]
                    if thingEnergy[i] < fragPrefilter:
                        filteredCenters.append(thingCenters.pop(i))
                        filteredEnergy.append(thingEnergy.pop(i))

                # Default Criteria: Distance; Default Method: Complete (As defined in the args section above)
                filteredCenters = np.array(filteredCenters)
                filteredEnergy = np.array(filteredEnergy)

                if filteredCenters.size != 0:
                    thingCluster=scip.fclusterdata(filteredCenters, t = clusterThreshold, criterion = criteria, method = method)
                    thingCenterspheres=sortCenter(filteredCenters, thingCluster, filteredEnergy)
                    allSpheres.extend(thingCenterspheres)
                    #If you want to remove all the single molecule clusters
                    allSpheres=purge(allSpheres, THRESHOLD)
                    #print(allSpheres)


                if os.path.exists(out_dir+'/'+'GeneratedPharma.json') == False:
                    #If you want to remove all the single molecule clusters
                    allSpheres=purge(allSpheres, THRESHOLD)
                    pharma=allJsonify(kind, allSpheres)
                    jason = json.dumps(pharma)
                    f = open(out_dir+"/"+"GeneratedPharma.json","w")
                    f.write(jason)
                    f.close()

                else:
                    allspheres = purge(allSpheres, THRESHOLD)
                    pharma = jsonify(kind, allSpheres)
                    # function to add to JSON
                    def write_json(data, filename=out_dir+'GeneratedPharma.json'):
                        with open(filename,'w') as f:
                            json.dump(data, f, indent=4)
                    y = {"receptor": str(open(recept).read()), "recname": recept}
                    with open(out_dir+'GeneratedPharma.json') as json_file:
                        data = json.load(json_file)
                        data.update(y)

                        temp = data['points']

                        temp.extend(pharma)

                    write_json(data)

    else:
        file = np.genfromtxt(args.inputFile, dtype=str, delimiter=',')
        ligands, receptors, outs = file[:, 0], file[:, 1], file[:, 2]

    #-------------------------------------------------------------------------------
    # RUN GNINA
    #-------------------------------------------------------------------------------
        for key in fragments:
            for i, ligand in enumerate(receptors):
                f = open("tempFragment.sdf", "w")
                f.write(fragments[key])
                f.close()
                subprocess.run(["gnina", "-r", recept, "-l", "tempFragment.sdf", "-o", out_dir + "/" + key+"Docked.sdf", "--num_modes", str(args.numMode), "--autobox_ligand", args.ligandFile])
                os.remove("tempFragment.sdf")
    #-------------------------------------------------------------------------------
    # GENERATE PHARMACOPHORE
    #-------------------------------------------------------------------------------
        for receptor in enumerate(receptors):
            os.chdir(outs[i])
            for kind in kinds:
                SMART=[]

                if(kind=="Hydrophobic"):
                    SMART=HYDROPHOBIC
                    nuPo= 6
                    maxEne= -2
                    criteria = "distance"
                    method = "ward"
                    clusterThreshold = 5
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
                elif(kind=="Aromatic"):
                    SMART= AROMATIC
                    nuPo= 6
                    maxEne= -2
                    criteria = 'distance'
                    method = 'ward'
                    clusterThreshold = 6
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
                elif(kind=="HydrogenDonor"):
                    SMART=HYDROGEN_DONOR
                    nuPo= 4
                    maxEne= -2
                    criteria = 'distance'
                    method = 'ward'
                    clusterThreshold = 5
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
                elif(kind=="HydrogenAcceptor"):
                    SMART=HYDROGEN_ACCEPTOR
                    nuPo= 4
                    maxEne= -2
                    criteria = 'distance'
                    method = 'ward'
                    clusterThreshold = 3
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
                elif(kind=="NegativeIon"):
                    SMART=NEGATIVE_ION
                    nuPo= 3
                    maxEne= -3
                    criteria = 'distance'
                    method = 'ward'
                    clusterThreshold = 1
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
                elif(kind=="PositiveIon"):
                    SMART=POSITIVE_ION
                    nuPo=6
                    maxEne= -1
                    criteria = 'distance'
                    method = 'complete'
                    clusterThreshold = 6
                    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}


                allSpheres=[]
                thingCenters, thingEnergy=[],[]
                thingFile = []
                filteredCenters, filteredEnergy = [],[]
                fileList = ['benzeneDocked.sdf','isopropylamineDocked.sdf', 'acetamideDocked.sdf', 'isopropanolDocked.sdf', 'imidazoleDocked.sdf', 'acetateDocked.sdf', 'isobutaneDocked.sdf', 'waterDocked.sdf', 'guanidineDocked.sdf']

                for aia in fileList:
                    thingCentersTem, thingEnergyTem=splitter(getCenters(os.path.join(out_dir, aia), SMART))
                    thingCenters.extend(thingCentersTem)
                    thingEnergy.extend(thingEnergyTem)
                    for i, k in enumerate(thingEnergyTem):
                        thingFile.append(sdfmol)
                    #Cluster info. This isn't really useful without processing it using 

                thingCluster=[]
                if((len(thingCenters)>0)):
                    # Defining a prefiltering of Affinities by Feature Type
                    things = [thingCenters, thingEnergy, thingFile]
                    fragmentPrefilt = args.fragmentPrefilter
                    for i, thing in enumerate(thingCenters):
                        fileIndex = fileList.index(thingFile[i])
                        fragPrefilter = fragmentPrefilt[fileIndex]
                        if thingEnergy[i] < fragPrefilter:
                            filteredCenters.append(thingCenters.pop(i))
                            filteredEnergy.append(thingEnergy.pop(i))

                    # Default Criteria: Distance; Default Method: Complete (As defined in the args section above)
                    filteredCenters = np.array(filteredCenters)
                    filteredEnergy = np.array(filteredEnergy)

                    if filteredCenters.size != 0:
                        thingCluster=scip.fclusterdata(filteredCenters, t = clusterThreshold, criterion = criteria, method = method)
                        thingCenterspheres=sortCenter(filteredCenters, thingCluster, filteredEnergy)
                        allSpheres.extend(thingCenterspheres)
                        #If you want to remove all the single molecule clusters
                        allSpheres=purge(allSpheres, THRESHOLD)
                        
                    if os.path.exists(out_dir+'/'+'GeneratedPharma.json') == False:

                        #If you want to remove all the single molecule clusters
                        allSpheres=purge(allSpheres, THRESHOLD)
                        pharma=allJsonify(kind, allSpheres)
                        jason = json.dumps(pharma)
                        f = open(out_dir+"/"+"GeneratedPharma.json","w")
                        f.write(jason)
                        f.close()

                    else:
                        allspheres = purge(allSpheres, THRESHOLD)
                        pharma = jsonify(kind, allSpheres)
                        # function to add to JSON
                        def write_json(data, filename=out_dir+'GeneratedPharma.json'):
                            with open(filename,'w') as f:
                                json.dump(data, f, indent=4)
                        
                        y = {"receptor": str(open(recept).read()), "recname": recept}

                        with open(out_dir+'GeneratedPharma.json') as json_file:
                            data = json.load(json_file)
                            data.update(y)

                            temp = data['points']

                            temp.extend(pharma)

                        write_json(data)


if __name__ == '__main__':
    main()

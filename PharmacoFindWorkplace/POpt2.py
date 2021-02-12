#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 20:29:01 2020

@author: samcho

NOTE: This is what to use once you have decided on the thresholds
It is almost identical to PharmacoFind, but it will only run once on each receptor.
So, this is what you want to get the final pharmacophore info. 

"""

import os
import subprocess
import csv 
import ast
import json
import numpy as np
import multiprocessing
import time

def get_Tani(energy, numPoints, PharmType, Methods, Criterias, clusterThresholds):
    data=[]
    field_names = ['Benchmark Name', 'Pharmacophore Name', 'Cluster Criteria', 'Cluster Method','Cluster Threshold', 'Energy Threshold', 'Cluster Size Threshold', 'Tanimoto', 'Precision', 'Recall', 'Match', 'Unmatch'] 

    currentDir = os.getcwd()
    parsedcwd = currentDir.split('/')
    benchmark = parsedcwd[-1]
    path = os.path.join(currentDir, "OptimalPharmacophores")
    os.makedirs(path, exist_ok=True)


    ida=subprocess.run(["python", "../PharmacoFindWithTanimoto.py", "-e", str(energy), "-c", Criterias, "-m", Methods, "-j", str(clusterThresholds), "-p", str(numPoints) , "-k", PharmType, '-w', path+'/'+PharmType+"GeneratedPharma.json", '-t', path+'/'+PharmType+"GeneratedPharma.json"], stdout=subprocess.PIPE, text=True)

    with open("../../data/"+PharmType + 'Tani.csv', 'a') as csvfile: 
        writer = csv.DictWriter(csvfile, fieldnames = field_names) 
        writer.writeheader() 
        writer.writerows(data)

if __name__ == '__main__':
    
    p1 = multiprocessing.Process(target = get_Tani(-3, 6, PharmType = 'Aromatic', Methods = 'ward', Criterias = 'distance', clusterThresholds = 6))
    p2 = multiprocessing.Process(target = get_Tani(-2, 4, PharmType = 'HydrogenAcceptor', Methods = 'ward', Criterias = 'distance', clusterThresholds = 3))
    p3 = multiprocessing.Process(target = get_Tani(-2, 4, PharmType = 'HydrogenDonor', Methods = 'ward', Criterias = 'distance', clusterThresholds = 5))
    p4 = multiprocessing.Process(target = get_Tani(-2, 6, PharmType = 'Hydrophobic', Methods = 'ward', Criterias = 'distance', clusterThresholds = 5))
    p5 = multiprocessing.Process(target = get_Tani(-3, 3, PharmType = 'NegativeIon', Methods = 'ward', Criterias = 'distance', clusterThresholds = 1))
    p6 = multiprocessing.Process(target = get_Tani(-1, 6, PharmType = 'PositiveIon', Methods = 'complete', Criterias = 'distance', clusterThresholds = 6))
    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()
    p6.start()

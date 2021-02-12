
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 00:35:41 2019
Program name: PharmacoFindWithTanimoto
What it does: It's basically PharmacoFind, but it also calculates the tanimoto value.
@author: samcho
NOTE: program is the base of GridSearchPharmacoFindWithTanimoto and OptimizedPharmacoFindWithTanimoto, generating a pharmacophore model for one receptor, and giving the tanimoto coefficients.
It'll make a pharmacophore model, then it'll calculate the similarity between the generated model (generatedPharma.json) and the true model (pharma.json)
"""
from __future__ import print_function  # ensures print function compatibility with Python3
import argparse
# safely deal with file paths
import rdkit
import math
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy as scip
import numpy as np
import json, sys
# library for loading the image
from rdkit import Chem, ForceField

import scipy
from scipy import spatial


# library for viewing the image


#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
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
    #pointList.extend(jsonify("Aromatic",aRing))
    #pointList.extend(jsonify("HydrogenAcceptor", hAccept))
   # pointList.extend(jsonify("HydrogenDonor", hDonor))
   # pointList.extend(jsonify("PositiveIon", pIon))
   # pointList.extend(jsonify("NegativeIon", nIon))
   # pointList.extend(jsonify("Hydrophobic", hPhobic))
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

"""
Puts it all together.
"""
def main():
    
    parser = argparse.ArgumentParser(description='Finds pharmacophore',formatter_class=argparse.RawDescriptionHelpFormatter)
    #Dock and include guanidine
    #Also, it looks like imidazole is the one that's trying to be positive.
    parser.add_argument("-r","--readFile",type=list, default=['benzeneDocked.sdf','isopropylamineDocked.sdf', 'acetamideDocked.sdf', 'isopropanolDocked.sdf', 'imidazoleDocked.sdf', 'acetateDocked.sdf', 'isobutaneDocked.sdf', 'waterDocked.sdf', 'guanidineDocked.sdf'], help="Filename to read from. (Default:'benzene.sdf')")
    parser.add_argument("-w","--writeFile", type=str,default="generatedPharma.json" ,help="Filename to write to. (Default: 'generatedPharma.txt')")
    parser.add_argument("-k", "--pharmType", type=str, default="Aromatic", help ="Name of the feature type")
    parser.add_argument("-p", "--numPoints", type=int, default=0, help="Threshold for Minimum number of points")
    parser.add_argument("-e", "--maxEnergy", type=float, default=0.0, help="Threshold for the maximum binding energy (we want it to be low, the more negative the better)")
    parser.add_argument("-f","--referenceFile",type=str, default="pharma.json",help="True Pharmacophore File.")
    parser.add_argument("-t","--testFile",type=str, default="generatedPharma.json",help="Generated Pharmacophore File. Default: generatedPharma.json")
    parser.add_argument("-c", "--clusterCriteria", type=str, default='distance', help="Defines the Cluster Criteria Used. Default: distance")
    parser.add_argument("-m", "--clusterMethod", type=str, default='complete', help="Defines the Cluster Method used. Default: complete")
    parser.add_argument("-j", "--clusterThreshold", type=int, default=1, help= "Threshold to apply when forming flat clusters")
    parser.add_argument("-b", "--fragmentPrefilter", type=float, default=[-3.49, -3.19, -2.66, -2.33, -3.00, -2.50, -3.00, -1.00, -3.00], help='Threshold of Affinity prefiltering. This removes fragments BEFORE clustering. Default: -2.0')
    args=parser.parse_args()
    kind=args.pharmType
    nuPo=args.numPoints
    maxEne=args.maxEnergy
    criteria = args.clusterCriteria
    method = args.clusterMethod
    clusterThreshold = args.clusterThreshold
    THRESHOLD={"type": kind,"numPoint": nuPo, "energy": maxEne}
    SMART=[]
    if(kind=="Hydrophobic"):
        SMART=HYDROPHOBIC
    elif(kind=="Aromatic"):
        SMART=AROMATIC
    elif(kind=="HydrogenDonor"):
        SMART=HYDROGEN_DONOR
    elif(kind=="HydrogenAcceptor"):
        SMART=HYDROGEN_ACCEPTOR
    elif(kind=="NegativeIon"):
        SMART=NEGATIVE_ION
    elif(kind=="PositiveIon"):
        SMART=POSITIVE_ION

    #Mess with distance thresholds, desired energy threshholds, all that good stuff.
    #Then maybe: filter out all with energy more than, let's say -2.
    #Coordinates, Energies
    fileList=args.readFile
    outFile=args.writeFile
    allSpheres=[]
    thingCenters, thingEnergy=[],[]
    thingFile = []
    filteredCenters, filteredEnergy = [],[]

    for aia in fileList:
        thingCentersTem, thingEnergyTem=splitter(getCenters(aia, SMART))
        thingCenters.extend(thingCentersTem)
        thingEnergy.extend(thingEnergyTem)
        for i, k in enumerate(thingEnergyTem):
            thingFile.append(aia)
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
            #print(allSpheres)


    pharma=allJsonify(kind, allSpheres)
    jason = json.dumps(pharma)
    f = open(outFile,"w")
    f.write(jason)
    f.close()
    
if __name__ == '__main__':
    main()

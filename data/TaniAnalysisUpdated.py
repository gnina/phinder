#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python
# coding: utf-8

import pandas as pd ## Table Transformation Library
import numpy as np ## Linear Algebra Library
import os ## Accessing the Operation System (OS)
import matplotlib.pyplot as plt ## Base Plotting Functions
import seaborn as sns ## High-Level Plotting Library


dataDir = os.getcwd()
files = os.listdir(dataDir)
print("Files in Directory: ", files)

csvFiles = []

for file in files:
    if file.endswith(".csv"):
        csvFiles.append(dataDir+"/"+file)
print("CSV Files: ", csvFiles, "\n")

mergedDataList = []

print("Cleaning Output")

for i,file in enumerate(csvFiles):

    data = pd.read_csv(file, sep=',')
    mergedDataList.append(data)

data = mergedDataList[0].append(mergedDataList[1:])
data = pd.DataFrame(data)
data = data[~data['Pharmacophore Name'].str.contains("Pharmacophore Name")]
data.to_csv('CombinedData.csv', index=False)
data = pd.read_csv('CombinedData.csv')

print("Output Cleaned")

print("Creating top 5 settings")
finalTable = []
pharmi = np.unique(data['Pharmacophore Name'])

for i in pharmi:
    top5 = data[data['Pharmacophore Name'] == str(i)].groupby(['Pharmacophore Name','Cluster Affinity Prefilter', 'Cluster Method', 'Cluster Threshold', 'Energy Threshold', 'Cluster Size Threshold', 'Cluster Criteria']).mean().reset_index().sort_values(by='Tanimoto', ascending=False).reset_index().loc[:4, :]
    finalTable.append(top5)

top5 = finalTable[0].append(finalTable[1:])
top5 = pd.DataFrame(top5)
top5 = top5.drop('index', axis=1)
top5.to_csv('top5.csv', index=False)


# In[62]:

print("Creating top settings")
finalTable = []

for i in pharmi:
    top5 = data[data['Pharmacophore Name'] == str(i)].groupby(['Pharmacophore Name', 'Cluster Method','Cluster Affinity Prefilter', 'Cluster Threshold', 'Energy Threshold', 'Cluster Size Threshold', 'Cluster Criteria']).mean('Tanimoto').reset_index().sort_values(by='Tanimoto', ascending=False).reset_index().iloc[:1]
    finalTable.append(top5)

top1 = finalTable[0].append(finalTable[1:])
top1 = pd.DataFrame(top1)
top1 = top1.drop('index', axis=1)
top1.to_csv('top1.csv', index=False)

print("Analysis Complete")

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd


# In[41]:


allTypes = {}
allSubparts = {}
for filePath in glob.glob(r"C:\Users\pakan\OneDrive\Documents\iGEM\*.csv"):
    plateNum = filePath[48:55]
    plate = pd.read_csv(filePath)
    #get part types for type distribution
    partType = plate["  Type "].values.tolist()
    partType = [ele.strip() for ele in partType]
    allTypes[plateNum] = partType 
    #get number of subparts for subpart distribution
    subpartNum = plate["  Subparts "].values.tolist()
    allSubparts[plateNum] = subpartNum
    


# In[47]:


#get all unique part types
types = list()
for plate in allTypes:
    for current in allTypes[plate]:
        newType = True
        for partType in types:
            if current == partType:
                newType = False
        if newType == True:
            types.append(current)
            
            


# In[48]:


#get distribution of part types in 2018 kit
typeCount = {}
for partType in types:
    typeCount[partType] = 0
    for plate in allTypes:
        for current in allTypes[plate]:
            if current == partType:
                typeCount[partType] += 1
print(typeCount)


# In[44]:


#get all unique subpart numbers
subparts = list()
for plate in allSubparts:
    for current in allSubparts[plate]:
        newNum = True
        for num in subparts:
            if current == num:
                newNum= False
        if newNum == True:
            subparts.append(current)


# In[49]:


#get distribution of subpart numbers in 2018 kit
subpartCount = {}
for num in subparts:
    subpartCount[num] = 0
    for plate in allSubparts:
        for current in allSubparts[plate]:
            if current == num:
                subpartCount[num] += 1
print(subpartCount)


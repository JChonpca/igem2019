#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import seaborn as sb


# In[2]:


#clean strain list part IDs
strainListPath = r"C:\Users\pakan\OneDrive\Documents\iGEM\Computation Team\Strain List 7-16.csv"
sL = pd.read_csv(strainListPath)
partID = [ele[1].strip() if type(ele)==list else None for ele in list(sL["PLASMID"].str.split("+"))]
sL["BBa"] = partID


# In[3]:


allTypes = {}
allSubparts = {}
masterPlate=[]
for filePath in glob.glob(r"C:\Users\pakan\OneDrive\Documents\iGEM\*.csv"):
    plateNum = filePath[48:55]
    plate = pd.read_csv(filePath)
    #get part types for type distribution
    partType = plate["  Type "].values.tolist()
    partType = [ele.strip() for ele in partType]
    allTypes[plateNum] = partType 
    plate["  Type "] = partType
    #get number of subparts for subpart distribution
    subpartNum = plate["  Subparts "].values.tolist()
    allSubparts[plateNum] = subpartNum
    plate["  Subparts "] = subpartNum
    #isolate part id, part type,and subpart count
    plate.rename(columns = {"  Part ": "BBa"}, inplace = True)
    parts = plate["BBa"].values.tolist()
    parts = [ele.strip() for ele in parts]
    plate["BBa"] = parts
    #test = pd.merge(sL, plate[["BBa","  Type ","  Subparts "]], on = "BBa")
   # print(plateNum)
    #print(plate[["BBa", "  Type ", "  Subparts "]])
    masterPlate.append(plate)
masterPlate=pd.concat(masterPlate)
masterPlate=masterPlate.drop_duplicates(subset="BBa")


# In[5]:


df=sL.merge(masterPlate[["BBa", "  Type ", "  Subparts "]],on="BBa")


# In[22]:


df.to_csv(r"C:\Users\pakan\Documents\iGEM\strainList.csv")


# In[6]:


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


# In[7]:


#get distribution of part types in 2018 kit
typeCountKit = {}
for partType in types:
    typeCountKit[partType] = 0
    for plate in allTypes:
        for current in allTypes[plate]:
            if current == partType:
                typeCountKit[partType] += 1
print(typeCountKit)
    


# In[8]:


#get distribution of part types from strain list
transformedTypes = df["  Type "]
typeCountProg = {}
for partType in types:
    typeCountProg[partType] = 0
    for current in transformedTypes:
        if current == partType:
            typeCountProg[partType] += 1
print(typeCountProg)


# In[9]:


#get percentages of types from strain list
typePercent = {}
for partType in types:
    fromStrainList = typeCountProg[partType]
    total = len(df["BBa"])
    percent = ((fromStrainList)/total)*100
    typePercent[partType] = percent
    
print(typePercent)
    


# In[10]:


#get percentages of types form 2018 kit
kitTypePer = {}
for partType in types:
    kitCount = typeCountKit[partType]
    total = len(masterPlate["BBa"])
    percent = ((kitCount)/total)*100
    kitTypePer[partType] = percent
print(kitTypePer)


# In[11]:


#round type percentages 
#prepare data frame for plotting
for per in typePercent:
    typePercent[per] = round(typePercent[per],2)
    
for per in kitTypePer:
    kitTypePer[per] = round(kitTypePer[per],2)
    
print(typePercent)
print(kitTypePer)
#from strain list
#typeData = pd.DataFrame(data = typePercent, index = [0])
#typeData = pd.DataFrame(data = kitTypePer, index = [1])
typeData = pd.DataFrame([kitTypePer, typePercent])

    


# In[12]:


print(typeData)


# In[13]:


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


# In[14]:


#get distribution of subpart numbers in 2018 kit
subpartCount = {}
for num in subparts:
    subpartCount[num] = 0
    for plate in allSubparts:
        for current in allSubparts[plate]:
            if current == num:
                subpartCount[num] += 1
print(subpartCount)


# In[15]:


#get distribution of subpart numbers in strain list
transformedSubparts = df["  Subparts "]
subpartProg = {}
for num in subparts:
    subpartProg[num] = 0
    for current in transformedSubparts:
        if current == num:
            subpartProg[num] += 1
print(subpartProg)


# In[16]:


#get ratios of sub part distribution
#get percentages of types
subpartPer = {}
for num in subparts:
    fromStrainList = subpartProg[num]
    total = len(df["  Subparts "])
    percent = ((fromStrainList)/total)*100
    subpartPer[num] = percent
    
print(subpartPer)


# In[17]:


#get percentages from 2018 kit
kitSubpartPer = {}
for num in subparts:
    kitCount = subpartCount[num]
    total = len(masterPlate["  Subparts "])
    percent = ((kitCount)/total)*100
    kitSubpartPer[num] = percent
print(kitSubpartPer)


# In[18]:


#round sub part values and prepare for plotting
for per in subpartPer:
    subpartPer[per] = round(subpartPer[per],2)
    
for per in kitSubpartPer:
    kitSubpartPer[per] = round(kitSubpartPer[per],2)

print(subpartPer)
print(kitSubpartPer)
subpartData = pd.DataFrame([kitSubpartPer, subpartPer])


# In[67]:


#sb.set(style="whitegrid")

# Load the example Titanic dataset
#test = sb.load_dataset("typeData")

# Draw a nested barplot to show survival for class and sex
#plot = sb.catplot(x="Part Types", y="Percent Total", hue="Legend", data= typeData,
                #height=6, kind="bar", palette="muted")
#plot.despine(left=True)
#plot.set_ylabels("Percent Total")


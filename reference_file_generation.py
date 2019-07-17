#!/usr/bin/env python
# coding: utf-8

# In[33]:


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import urllib
import xml.etree.ElementTree as ET
import os
import os.path


# In[2]:


fastaList=list(SeqIO.parse(r"C:\Users\pakan\Documents\iGEM\All_Parts.txt", "fasta"))
fastaDict={ele.id:ele for ele in fastaList}


# In[3]:


with open(r"C:\Users\pakan\Documents\iGEM\pSB1C3.gb",'r') as file_handle:
    record_dict = SeqIO.to_dict(SeqIO.parse(file_handle, 'gb'))
gbkFile = record_dict[list(record_dict.keys())[0]]
#glob.glob(r"C:\Users\pakan\OneDrive\Documents\iGEM\*.csv")


# In[4]:


plates={}
for infile_loc in glob.glob(r"C:\Users\pakan\OneDrive\Documents\iGEM\*.csv"):
    plateNum = infile_loc[48:55]
    plate = pd.read_csv(infile_loc)
    plateParts=plate["  Part "].values.tolist()
    plateParts=[ele.strip() for ele in plateParts]
    plates[plateNum]=plateParts            
    


# In[5]:


allPartsSeqs={}
for plate in plates:
    print(plate)
    for curPart in plates[plate]:
        try:
            sequence=fastaDict[curPart].seq.upper()
            #record = fastaDict[curPart].seq.upper() + gbkFile
        except KeyError:
            print(".",end="")
            igemURL = f'http://parts.igem.org/cgi/xml/part.cgi?part={curPart}'
            tree = ET.parse(urllib.request.urlopen(igemURL))
            root = tree.getroot()
            for seq in root.iter('seq_data'):
                #record = seq.text + gbkFile
                sequence=seq.text
                sequence=sequence.replace("\n", "")
                sequence=Seq(sequence)
        
        allPartsSeqs[curPart]=sequence
        


# In[6]:


sequenceRecords = list()
for x in allPartsSeqs:
    record = allPartsSeqs[x].upper() + gbkFile
    record.id = x;
    sequenceRecords.append(record)


# In[7]:


xmlCheck = {}
for plate in plates:
    print(plate)
    for curPart in plates[plate]:
        print(".",end="")
        igemURL = f'http://parts.igem.org/cgi/xml/part.cgi?part={curPart}'
        tree = ET.parse(urllib.request.urlopen(igemURL))
        root = tree.getroot()
        for seq in root.iter('seq_data'):
            sequence=seq.text
            sequence=sequence.replace("\n", "")
            sequence=Seq(sequence)
        xmlCheck[curPart] = sequence.upper()


# In[8]:


for part in xmlCheck:
    xmlSeq = xmlCheck[part]
    txtSeq = allPartsSeqs[part]
    if xmlSeq != txtSeq:
        #print(part)
        record = xmlCheck[part] + gbkFile
        record.id = part + " (2)"
        sequenceRecords.append(record)
        


# In[9]:


outputPath = r"C:\Users\pakan\Documents\iGEM\referenceFiles"
if not os.path.exists(outputPath):
    os.makedirs(outputPath)


# In[36]:


outputPath = 'C:/Users/pakan/Documents/iGEM/referenceFiles/'
for reference in sequenceRecords:
    fileName = reference.id+".fa"
    completeName = os.path.join(outputPath, fileName)
    file1 = open(completeName, "w")
    SeqIO.write(reference, file1, "fasta")
        


# In[1]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import urllib
import xml.etree.ElementTree as ET
import os
import os.path
from Bio.Alphabet import generic_dna


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
    


# In[23]:


allPartsSeqs={}
for plate in plates:
    print(plate)
    for curPart in plates[plate]:
        igemURL = f'http://parts.igem.org/cgi/xml/part.cgi?part={curPart}'
        tree = ET.parse(urllib.request.urlopen(igemURL))
        root = tree.getroot()
        for seq in root.iter('seq_data'):
            #record = seq.text + gbkFile
            sequence=seq.text
            sequence=sequence.replace("\n", "")
            sequence=Seq(sequence)
        
        allPartsSeqs[curPart]=sequence
        


# In[25]:


from Bio.SeqFeature import SeqFeature, FeatureLocation
sequenceRecords = list()
for x in allPartsSeqs:
    record = allPartsSeqs[x].upper() + gbkFile
    featLoc=FeatureLocation(0,len(allPartsSeqs[x]),1)
    record.features.append(SeqFeature(featLoc, type='Region',qualifiers={"label": x}))
    record.id = x;
    sequenceRecords.append(record)
    


# In[26]:


outputPath = r"C:\Users\pakan\Documents\iGEM\referenceFiles"
if not os.path.exists(outputPath):
    os.makedirs(outputPath)


# In[27]:


from Bio.Alphabet import generic_dna
outputPath = 'C:/Users/pakan/Documents/iGEM/referenceFiles/'
for reference in sequenceRecords:
    reference.seq.alphabet=generic_dna
    fileName = reference.id+".gbk"
    completeName = os.path.join(outputPath, fileName)
    file1 = open(completeName, "w")
    SeqIO.write(reference, file1, "genbank")
    file1.close()
        


#!/usr/bin/env python
# coding: utf-8

# # Importing libraries

# In[3]:


import pandas as pd
import csv
import Bio
import sys
import numpy as np
from Bio.SeqUtils import ProtParamData  # Local
from Bio.SeqUtils import IsoelectricPoint  # Local
from Bio.Seq import Seq
from Bio.Data import IUPACData
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils import GC
import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


# # Import the dataset

# In[21]:


df=pd.read_csv('X.csv')
df


# In[22]:


df['Sequence'] = df['Sequence'].str.replace(r'\s+', '') #REPLCING LINE BREAKS


# In[55]:


seq=df['Sequence'].values.tolist() #creating a list of the sequence
seq


# # Calculating Atomic Forumula

# In[57]:


# Define a function to calculate the atomic formula of a peptide
def calculate_atomic_formula(sequence):
    mol = Chem.MolFromSequence(sequence)
    if mol is not None:
        formula = CalcMolFormula(mol)
        return formula
    else:
        return "Invalid sequence"


# In[56]:


# Create an empty dictionary to store the results
results = {}

# Loop through each sequence and calculate its atomic formula
for i in seq:
    formula = calculate_atomic_formula(i)
    results[i] = formula

# Convert the results to a Pandas DataFrame
import pandas as pd
df1 = pd.DataFrame.from_dict(results, orient="index", columns=["Atomic Formula"])

# Print the results
print(df1)


# In[25]:


# create a dictionary that maps digits to subscript characters
subscripts = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

# loop through each row in the dataframe and convert any string column to subscript form
for col in df.select_dtypes(include='object').columns:
    df1['Atomic Formula'] = df1['Atomic Formula'].str.translate(subscripts)

df1


# In[26]:


AF=df1['Atomic Formula'].values.tolist() #creating a list of the sequence
AF


# In[27]:


df['Atomic Formula']=AF
df


# # Calculating Molecular Weight

# In[28]:


from Bio.SeqUtils import seq1

list_of_MW=[]

for i in seq:
    MW=Bio.SeqUtils.molecular_weight(seq=i, seq_type='protein', double_stranded=False, circular=False, monoisotopic=False)
    if MW==1:
        list_of_MW=[] #empty list of mol weight
    else:
        list_of_MW.append(MW) #adding the MW of the subsequent aa sequence into the list of mol weight
        
print(list_of_MW)


# In[29]:


df['Molecular Weight']=list_of_MW


# In[33]:


df


# In[30]:


seq_analysis=tuple(seq) #converting list into tuple data structure
seq_analysis


# # Calculating GRAVY(Grand Average Value of Hydropathy)

# In[31]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis

GRAVY=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.gravy()
    if y==1:
        GRAVY=[] #empty list of GRAVY
    else:
        GRAVY.append(y) #adding the GRAVY of the subsequent aa sequence into the list of GRAVY
print(GRAVY)


# In[32]:


df['GRAVY']=GRAVY
df


# # Calculating Aromaticity

# In[33]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis

LIST_OF_AROMATICITY=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.aromaticity()
    if y==1:
        LIST_OF_AROMATICITY=[] #empty list of SEQ COUNT
    else:
        LIST_OF_AROMATICITY.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_AROMATICITY)


# In[34]:


df['Aromaticty']=LIST_OF_AROMATICITY


# # Calculating Instability Index

# In[35]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis

LIST_OF_II=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.instability_index()
    if y==1:
        LIST_OF_II=[] 
    else:
        LIST_OF_II.append(y) 
print(LIST_OF_II)


# In[36]:


df['Instability index']=LIST_OF_II


# # Calculating Molecular Extinction Coefficient

# In[37]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis

LIST_OF_MEC=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.molar_extinction_coefficient()
    if y==1:
        LIST_OF_MEC=[] #empty list of SEQ COUNT
    else:
        LIST_OF_MEC.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_MEC)


# In[38]:


MEC_DF = pd.DataFrame(LIST_OF_MEC) 
print(MEC_DF)


# In[39]:


F5=MEC_DF.rename({0:'Reduced Mol extinction coeff',1:'Oxidized Mol extinction coeff'}, axis=1)
F5


# In[40]:


res=[df,F5] #creating the function of set of dataframes to be joined
df=pd.concat(res, axis=1)
df


# # Calculating Isoelectric Point

# In[41]:


LIST_OF_IP=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.isoelectric_point()
    if y==1:
        LIST_OF_IP=[] #empty list of SEQ COUNT
    else:
        LIST_OF_IP.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_IP)


# In[42]:


df['Isoelectric point'] =LIST_OF_IP
df


# # Calculating Peptide Length

# In[43]:


LENGTH_AA=[]

for i in range(0,282):
    y=len(seq_analysis[i])
    
    if y==1:
        LENGTH_AA=[] #empty list of SEQ COUNT
    else:
        LENGTH_AA.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LENGTH_AA)


# In[44]:


df['Peptide length'] =LENGTH_AA
df


# In[45]:


col = df.pop('Peptide length')
df.insert(2, 'Peptide length', col)
df


# In[50]:


df.to_csv('X',index=False)


# # Calculating Amino Acid Frequency

# In[51]:


#Calculation of amino acid count of peptide sequence

from Bio.SeqUtils.ProtParam import ProteinAnalysis

LIST_OF_SEQ_COUNT=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.count_amino_acids()
    if y==1:
        LIST_OF_SEQ_COUNT=[] #empty list of SEQ COUNT
    else:
        LIST_OF_SEQ_COUNT.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_SEQ_COUNT)


# In[52]:


SEQ_COUNT_DF = pd.DataFrame(LIST_OF_SEQ_COUNT) #creating dataframe of GRAVY of the sequences from the list of GRAVY
print(SEQ_COUNT_DF)


# In[48]:


df2=pd.read_csv('X')


# In[54]:


res=[df2,SEQ_COUNT_DF] #creating the function of set of dataframes to be joined
df_seqc=pd.concat(res, axis=1)
df_seqc


# In[55]:


df_seqc.to_csv('AA_Freq.csv',index=False)


# # Calculating Secondary Structure Fractions

# In[56]:


#for Secondary structure fraction of AA


LIST_OF_AA_SS=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.secondary_structure_fraction()
    if y==1:
        LIST_OF_AA_SS=[] #empty list of SEQ COUNT
    else:
        LIST_OF_AA_SS.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_AA_SS)


# In[57]:


#creating dataframe of Secondary Structure Fraction of the peptide sequences 

AA_PER_SS = pd.DataFrame(LIST_OF_AA_SS) 
ss=AA_PER_SS.rename({0:"Helix",1:"Turn",2:"Sheet"},axis=1)
ss


# In[58]:


res2=[df2,ss] #creating the function of set of dataframes to be joined
df_ss=pd.concat(res2, axis=1)
df_ss


# In[59]:


df_ss.to_csv("X_SecondaryStr.csv",index=False)


# # Calculating ProtScale Values for Hydropathy Plot

# In[60]:


LIST_OF_AA_PS=[]

for i in range(0,282):
    analysed_seq = ProteinAnalysis(str(seq_analysis[i]))
    y=analysed_seq.protein_scale(ProtParamData.kd, 9, 0.4)
    if y==1:
        LIST_OF_AA_PS=[] #empty list of SEQ COUNT
    else:
        LIST_OF_AA_PS.append(y) #adding the SEQ COUNT of the subsequent aa sequence into the list of SEQ COUNT
print(LIST_OF_AA_PS)


# In[61]:


AA_PS=pd.DataFrame(LIST_OF_AA_PS)


# In[62]:


res3=[df2,AA_PS] #creating the function of set of dataframes to be joined
df_prot_scale=pd.concat(res3, axis=1)
df_prot_scale


# In[63]:


df_prot_scale.to_csv('V3_Peptide_protscale.csv',index=False)


# In[46]:


df_pfeat=pd.read_csv('V3_Peptide_Pfeature_Physico.csv')
df_pfeat


# In[49]:


res4=[df2,df_pfeat] #creating the function of set of dataframes to be joined
df_pfeat_combined=pd.concat(res4, axis=1)
df_pfeat_combined


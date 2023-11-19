#!/usr/bin/env python
# coding: utf-8

# # Importing Libraries

# In[4]:


import csv
import pandas as pd
import matplotlib.pyplot as plot
import os
import matplotlib.font_manager


# In[5]:


from matplotlib import rcParams
rcParams['axes.titlepad'] = 20 


# # Changing Font

# In[20]:


plot.rcParams['font.family'] = 'DeJavu Serif'
plot.rcParams['font.serif'] = ['Times New Roman']


# # Loading the dataset

# In[7]:


df_pep=pd.read_csv('Peptide_protscale.csv')


# # Removing unnecessary columns

# In[9]:


df_pep2=df_pep.drop([('AAP ID','Sequence')], axis=1)


# # Converting the dataframe into numpy array

# In[22]:


arr=df_pep2.to_numpy()
arr


# # Plotting

# In[24]:


import seaborn as sns
import matplotlib.pyplot as plot

# Set style and context for Seaborn
sns.set_style("whitegrid")
sns.set_context("poster")
sns.set(font='Times New Roman')

# Create the figure and axes
fig, ax = plot.subplots(figsize=(10, 10), dpi=1200)

# Set title and axis labels
ax.set_title("K&D Hydropathy Plot of Peptide", fontsize=20, fontweight='bold', color='blue')
ax.set_ylabel('Hydrophobicity', fontsize=15, color='blue')
ax.set_xlabel('Residue Number', fontsize=15, color='blue')

# Plot the data
x = list(range(1, len(arr[0]) +1))
for i in range(len(arr)):
    sns.lineplot(x=x, y=arr[i], ax=ax)

# Save and show the plot
plot.savefig("Peptide_hydropathy plots_2.png")
plot.show()


# In[11]:


df_pro=pd.read_csv('V3_Protein_protscale.csv')


# In[12]:


df_pro1=df_pro.drop(['AAP ID'], axis=1)


# In[13]:


df_pro2=df_pro1.drop(['Sequence'], axis=1)


# In[14]:


arr2=df_pro2.to_numpy()
arr2


# In[19]:


import seaborn as sns
import matplotlib.pyplot as plot

# Set style and context for Seaborn
sns.set_style("whitegrid")
sns.set_context("poster")
sns.set(font='Times New Roman')

# Create the figure and axes
fig, ax = plot.subplots(figsize=(10, 10), dpi=1200)

# Set title and axis labels
ax.set_title("K&D Hydropathy Plot of Protein", fontsize=20, fontweight='bold', color='blue')
ax.set_ylabel('Hydrophobicity', fontsize=15, color='blue')
ax.set_xlabel('Residue Number', fontsize=15, color='blue')

# Plot the data
x = list(range(1, len(arr2[0]) +1))
for i in range(len(arr2)):
    sns.lineplot(x=x, y=arr2[i], ax=ax)

# Save and show the plot
plot.savefig("Protein_hydropathy plots_2.png")
plot.show()


# In[ ]:





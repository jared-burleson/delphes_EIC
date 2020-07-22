#!/usr/bin/env python
# coding: utf-8

# In[1]:


import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib.ticker import MultipleLocator
import warnings
warnings.filterwarnings('ignore')


# In[80]:


def error_calculate(data_array):
    raw_error = np.sqrt(data_array)
    rel_error = raw_error/data_array
    norm_data_array = data_array/sum(data_array)
    norm_rel_error = rel_error*norm_data_array
    norm_rel_error = np.nan_to_num(norm_rel_error)
    
    return norm_rel_error

def error_bar_hist_plotting(data, bin_count, range_min, range_max, xlabel, ylabel, title, color):
    plt.figure(0)
    (counts, bin_values, patches) = plt.hist(data, bins=bin_count, range=(range_min, range_max), histtype='step', color='white');
    counts = np.append(counts, 0)

    plt.figure(1)
    plt.errorbar(bin_values, counts/sum(counts), yerr=error_calculate(counts),ecolor='black',capsize=2,drawstyle='steps-mid',color=color)
    plt.xlim(range_min, range_max)
    plt.ylim(0,max(max(counts/sum(counts))*1.05,plt.gca().get_ylim()[1]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.minorticks_on()

    return 0

def TrimmedPT(row):
    px_value = row['Jet.TrimmedP4[5].fX']
    py_value = row['Jet.TrimmedP4[5].fY']
    pt_value = np.sqrt(px_value**2 + py_value**2)
    return pt_value

def TrimmedM(row):
    px_value = row['Jet.TrimmedP4[5].fX']
    py_value = row['Jet.TrimmedP4[5].fY']
    pz_value = row['Jet.TrimmedP4[5].fZ']
    e_value  = row['Jet.TrimmedP4[5].fE']
    if((e_value**2 - (px_value**2 + py_value**2 + pz_value**2)) < 0):
        m_value = -1
    else:
        m_value = np.sqrt(e_value**2 - (px_value**2 + py_value**2 + pz_value**2))
    return m_value

def PrunedPT(row):
    px_value = row['Jet.PrunedP4[5].fX']
    py_value = row['Jet.PrunedP4[5].fY']
    pt_value = np.sqrt(px_value**2 + py_value**2)
    return pt_value

def PrunedM(row):
    px_value = row['Jet.PrunedP4[5].fX']
    py_value = row['Jet.PrunedP4[5].fY']
    pz_value = row['Jet.PrunedP4[5].fZ']
    e_value  = row['Jet.PrunedP4[5].fE']
    if((e_value**2 - (px_value**2 + py_value**2 + pz_value**2)) < 0):
        m_value = -1
    else:
        m_value = np.sqrt(e_value**2 - (px_value**2 + py_value**2 + pz_value**2))
    return m_value


# In[93]:


jet_data_CCDIS = file_data_CCDIS.pandas.df(['Jet.PT','Jet.Eta','Jet.Phi','Jet.Mass','Jet.Flavor'], flatten=True)


# In[ ]:





# In[ ]:





# In[42]:


error_bar_hist_plotting(jet_data_NCDIS['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:blue')
error_bar_hist_plotting(jet_data_CCDIS['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:red')
plt.legend(['NC DIS','CC DIS'])
plt.show()

error_bar_hist_plotting(jet_data_NCDIS[jet_data_NCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Charm Jet pT','tab:orange')
error_bar_hist_plotting(jet_data_CCDIS[jet_data_CCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Charm Jet pT','tab:green')
plt.legend(['NC DIS','CC DIS'])
plt.show()

error_bar_hist_plotting(jet_data_NCDIS['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:blue')
error_bar_hist_plotting(jet_data_CCDIS['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:red')
error_bar_hist_plotting(jet_data_NCDIS[jet_data_NCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Charm Jet pT','tab:orange')
error_bar_hist_plotting(jet_data_CCDIS[jet_data_CCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Charm Jet pT','tab:green')
plt.legend(['All Jets (NC)','All Jets (CC)','Charm Jets (NC)','Charm Jets (CC)'])
plt.show()

error_bar_hist_plotting(jet_data_NCDIS[(jet_data_NCDIS['Jet.Flavor'] < 4) | (jet_data_NCDIS['Jet.Flavor'] == 21)]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:blue')
error_bar_hist_plotting(jet_data_CCDIS[(jet_data_CCDIS['Jet.Flavor'] < 4) | (jet_data_CCDIS['Jet.Flavor'] == 21)]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:red')
error_bar_hist_plotting(jet_data_NCDIS[jet_data_NCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:orange')
error_bar_hist_plotting(jet_data_CCDIS[jet_data_CCDIS['Jet.Flavor'] == 4]['Jet.PT'],45,5,50,'Jet pT (GeV)','Frequency','Jet pT','tab:green')
plt.legend(['Light Jets (NC)','Light Jets (CC)','Charm Jets (NC)','Charm Jets (CC)'])
plt.savefig('Graphics/Charm_Light_Jet_pT_NCCC_Comparison.jpg',dpi=300)
plt.show()


# In[ ]:


file_data_CCDIS = uproot.open("D:\College/JetStudy/EICStudy/EIC_CCDIS_SubJet_TrimmedStudy_output-3.root")["Delphes"]


# ### Trimmed Jet Study

# In[78]:


jet_data_CCDIS_trimmed = file_data_CCDIS.pandas.df(['Jet.TrimmedP4[5]'], flatten=True)
jet_data_CCDIS_trimmed.reset_index(inplace=True, drop=True)
jet_data_CCDIS_trimmed["Jet.TrimmedP4[5].pT"] = jet_data_CCDIS_trimmed.apply(TrimmedPT, axis=1)
jet_data_CCDIS_trimmed["Jet.TrimmedP4[5].M"] = jet_data_CCDIS_trimmed.apply(TrimmedM, axis=1)

drop_indexes_CCDIS_trimmed = []
single_drop_indexes_CCDIS_trimmed = []
for i in range(0,len(jet_data_CCDIS_trimmed)):
    if(i % 5 == 0):
        drop_indexes_CCDIS_trimmed.append(i)
    if(jet_data_CCDIS_trimmed['Jet.TrimmedP4[5].fX'][i] == 0.0 and jet_data_CCDIS_trimmed['Jet.TrimmedP4[5].fY'][i] == 0.0):
        drop_indexes_CCDIS_trimmed.append(i)
    if(i % 5 != 0):
        single_drop_indexes_CCDIS_trimmed.append(i)
        
subjet_data_CCDIS_trimmed = jet_data_CCDIS_trimmed.drop(drop_indexes_CCDIS_trimmed)
subjet_data_CCDIS_trimmed.reset_index(inplace=True, drop=True)

single_jet_data_CCDIS_trimmed = jet_data_CCDIS_trimmed.drop(single_drop_indexes_CCDIS_trimmed)
single_jet_data_CCDIS_trimmed.reset_index(inplace=True, drop=True)


# In[ ]:





# ### Pruned Jet Study

# In[86]:


jet_data_CCDIS_Pruned = file_data_CCDIS.pandas.df(['Jet.PrunedP4[5]'], flatten=True)
jet_data_CCDIS_Pruned.reset_index(inplace=True, drop=True)

jet_data_CCDIS_Pruned["Jet.PrunedP4[5].pT"] = jet_data_CCDIS_Pruned.apply(PrunedPT, axis=1)
jet_data_CCDIS_Pruned["Jet.PrunedP4[5].M"] = jet_data_CCDIS_Pruned.apply(PrunedM, axis=1)

jet_data_CCDIS_Pruned

single_drop_indexes_CCDIS_Pruned = []
drop_indexes_CCDIS_Pruned = []
for i in range(0,len(jet_data_CCDIS_Pruned)):
    if(i % 5 == 0):
        drop_indexes_CCDIS_Pruned.append(i)
    if(jet_data_CCDIS_Pruned['Jet.PrunedP4[5].fX'][i] == 0.0 and jet_data_CCDIS_Pruned['Jet.PrunedP4[5].fY'][i] == 0.0):
        drop_indexes_CCDIS_Pruned.append(i)
        single_drop_indexes_CCDIS_Pruned.append(i)
    if(i % 5 != 0):
        single_drop_indexes_CCDIS_Pruned.append(i)
        
subjet_data_CCDIS_Pruned = jet_data_CCDIS_Pruned.drop(drop_indexes_CCDIS_Pruned)
subjet_data_CCDIS_Pruned.reset_index(inplace=True, drop=True)

single_jet_data_CCDIS_Pruned = jet_data_CCDIS_Pruned.drop(single_drop_indexes_CCDIS_Pruned)
single_jet_data_CCDIS_Pruned.reset_index(inplace=True, drop=True)


# In[ ]:





# ### Plotting (Individual and Comparison)

# In[103]:


# Combined Jet pT Plots
error_bar_hist_plotting(jet_data_CCDIS['Jet.PT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'tab:blue')
error_bar_hist_plotting(single_jet_data_CCDIS_trimmed['Jet.TrimmedP4[5].pT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'tab:orange')
error_bar_hist_plotting(single_jet_data_CCDIS_Pruned['Jet.PrunedP4[5].pT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'tab:green')
plt.legend(['Jet pT','TrimmedJet pT','Pruned pT'])
plt.show()

# Combined Jet Mass Plots
error_bar_hist_plotting(jet_data_CCDIS['Jet.Mass'], 51, 0, 10.2, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:blue')
error_bar_hist_plotting(single_jet_data_CCDIS_trimmed['Jet.TrimmedP4[5].M'], 51, 0, 10.2, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:orange')
error_bar_hist_plotting(single_jet_data_CCDIS_Pruned['Jet.PrunedP4[5].M'], 51, 0, 10.2, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:green')
plt.xlim(0,10)
plt.legend(['Jet','TrimmedJet','Pruned'])
plt.show()

# Trimmed Jet Eta Plots

# Trimmed Jet Phi Plots

# Trimmed Jet Mass Plots

# Trimmed SubJet pT Plots

# Trimmed SubJet Eta Plots

# Trimmed SubJet Phi Plots

# Trimmed SubJet Mass Plots
# error_bar_hist_plotting(subjet_data_CCDIS_trimmed['Jet.TrimmedP4[5].M'], 40, 0, 2, 'SubJet M (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Trimmed SubJet M', 'tab:blue')
# plt.show()
# error_bar_hist_plotting(subjet_data_CCDIS_trimmed['Jet.TrimmedP4[5].M'], 25, 0, .5, 'SubJet M (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Trimmed SubJet M', 'tab:blue')
# plt.show()

# Pruned Jet pT Plots

# Pruned Jet Eta Plots

# Pruned Jet Phi Plots

# Pruned Jet Mass Plots

# Pruned SubJet pT Plots

# Pruned SubJet Eta Plots

# Pruned SubJet Phi Plots

# Pruned SubJet Mass Plots
# error_bar_hist_plotting(subjet_data_CCDIS_Pruned['Jet.PrunedP4[5].M'], 40, 0, 2, 'SubJet M (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Pruned SubJet M', 'tab:blue')
# plt.show()
# error_bar_hist_plotting(subjet_data_CCDIS_Pruned['Jet.PrunedP4[5].M'], 25, 0, .5, 'SubJet M (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Pruned SubJet M', 'tab:blue')
# plt.show()


# In[104]:


error_bar_hist_plotting(single_jet_data_CCDIS_trimmed['Jet.TrimmedP4[5].M'], 51, 0, 10.2, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:orange')
plt.show()


# In[ ]:





# This file is a generic file for studying SubJets from our EIC_DIS Simulations using different methods (such as Trimming or Pruning)
# We use the library class "uproot" to import the ROOT file and turn it into a pandas dataframe 
# We can then perform functions/computations, cuts, and plotting on the data

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib.ticker import MultipleLocator
import warnings
warnings.filterwarnings('ignore')

# Functions used for plotting (error calculation and plotting) and generating pT/mass 

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


# Importing ROOT file and making dataframe with Jet variables

file_data = uproot.open("INPUT_FILE.root")["Delphes"]
jet_data = file_data.pandas.df(['Jet.PT','Jet.Eta','Jet.Phi','Jet.Mass','Jet.Flavor'], flatten=True)
jet_data.reset_index(inplace=True, drop=True)

# Trimmed Jet Study

jet_data_trimmed = file_data.pandas.df(['Jet.TrimmedP4[5]'], flatten=True)
jet_data_trimmed.reset_index(inplace=True, drop=True)
jet_data_trimmed["Jet.TrimmedP4[5].pT"] = jet_data_trimmed.apply(TrimmedPT, axis=1)
jet_data_trimmed["Jet.TrimmedP4[5].M"] = jet_data_trimmed.apply(TrimmedM, axis=1)

drop_indexes_trimmed = []
single_drop_indexes_trimmed = []
for i in range(0,len(jet_data_trimmed)):
    if(i % 5 == 0):
        drop_indexes_trimmed.append(i)
    if(jet_data_trimmed['Jet.TrimmedP4[5].fX'][i] == 0.0 and jet_data_trimmed['Jet.TrimmedP4[5].fY'][i] == 0.0):
        drop_indexes_trimmed.append(i)
    if(i % 5 != 0):
        single_drop_indexes_trimmed.append(i)
        
subjet_data_trimmed = jet_data_trimmed.drop(drop_indexes_trimmed)
subjet_data_trimmed.reset_index(inplace=True, drop=True)

single_jet_data_trimmed = jet_data_trimmed.drop(single_drop_indexes_trimmed)
single_jet_data_trimmed.reset_index(inplace=True, drop=True)

# Pruned Jet Study

jet_data_Pruned = file_data.pandas.df(['Jet.PrunedP4[5]'], flatten=True)
jet_data_Pruned.reset_index(inplace=True, drop=True)

jet_data_Pruned["Jet.PrunedP4[5].pT"] = jet_data_Pruned.apply(PrunedPT, axis=1)
jet_data_Pruned["Jet.PrunedP4[5].M"] = jet_data_Pruned.apply(PrunedM, axis=1)

jet_data_Pruned

single_drop_indexes_Pruned = []
drop_indexes_Pruned = []
for i in range(0,len(jet_data_Pruned)):
    if(i % 5 == 0):
        drop_indexes_Pruned.append(i)
    if(jet_data_Pruned['Jet.PrunedP4[5].fX'][i] == 0.0 and jet_data_Pruned['Jet.PrunedP4[5].fY'][i] == 0.0):
        drop_indexes_Pruned.append(i)
        single_drop_indexes_Pruned.append(i)
    if(i % 5 != 0):
        single_drop_indexes_Pruned.append(i)
        
subjet_data_Pruned = jet_data_Pruned.drop(drop_indexes_Pruned)
subjet_data_Pruned.reset_index(inplace=True, drop=True)

single_jet_data_Pruned = jet_data_Pruned.drop(single_drop_indexes_Pruned)
single_jet_data_Pruned.reset_index(inplace=True, drop=True)

# Plotting (Comparison)

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

# More plots can be added to look at individual subjets or such, but the code is omitted because it is very simple to reproduce

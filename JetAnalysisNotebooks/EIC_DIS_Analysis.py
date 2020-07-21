# This file is a generic file for studying Jets from our EIC_DIS Simulations
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

# Importing ROOT file and creating a DataFrame

df_list = []
data = uproot.pandas.iterate("INPUT_FILE.root", "Delphes", ["Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.T", "Jet.Mass", "Jet.DeltaEta", "Jet.DeltaPhi", 'Jet.Flavor', 'Jet.FlavorAlgo', 'Jet.FlavorPhys', 'Jet.BTag', 'Jet.BTagAlgo', 'Jet.BTagPhys', 'Jet.TauTag', 'Jet.TauWeight', 'Jet.Charge', 'Jet.EhadOverEem', 'Jet.NCharged', 'Jet.NNeutrals', 'Jet.Beta', 'Jet.BetaStar', 'Jet.MeanSqDeltaR', 'Jet.PTD', 'Jet.SoftDroppedJet', 'Jet.SoftDroppedSubJet1', 'Jet.SoftDroppedSubJet2', 'Jet.NSubJetsTrimmed', 'Jet.NSubJetsPruned', 'Jet.NSubJetsSoftDropped', 'Jet.ExclYmerge23', 'Jet.ExclYmerge34', 'Jet.ExclYmerge45', 'Jet.ExclYmerge56', 'Jet.Constituents', 'Jet.Particles', 'Jet.Area'], flatten=True)
for dataframe in data:
    df_list.append(dataframe)
jet_data = pd.concat(df_list)
jet_data.reset_index(inplace=True, drop=True)
jet_data.to_csv('jet_data_EIC_DIS.csv')

# Jet Flavor Cutting

# Charm-Jets
jet_data_charm = jet_data[jet_data['Jet.Flavor'] == 4]
jet_data_charm.reset_index(inplace=True, drop=True)

# Light and Gluon Jets
jet_data_light = jet_data[(jet_data['Jet.Flavor'] < 4) | (jet_data['Jet.Flavor'] == 21)]
jet_data_light.reset_index(inplace=True, drop=True)

# Bottom-Jets
jet_data_bottom = jet_data[jet_data['Jet.Flavor'] == 5]
jet_data_bottom.reset_index(inplace=True, drop=True)

# Functions used for plotting (error calculation and plotting) 

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

# Section devoted to plotting

# Analyzing Jet pT
error_bar_hist_plotting(jet_data['Jet.PT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'black')
error_bar_hist_plotting(jet_data_charm['Jet.PT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'tab:blue')
error_bar_hist_plotting(jet_data_light['Jet.PT'], 45, 5, 50, 'Jet pT (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet pT', 'tab:red')
plt.legend(['All Jets','Charm Jets','Light Jets'],title='Jet Flavors')
plt.grid()
plt.savefig('Graphics/jet_pt.jpg',dpi=300)
plt.show()
plt.clf()

# Analyzing Jet eta
error_bar_hist_plotting(jet_data['Jet.Eta'], 50, -5, 5, 'Jet Eta', 'Frequency', 'DelphesPythia8 EIC Event: Jet Eta', 'black')
error_bar_hist_plotting(jet_data_charm['Jet.Eta'], 50, -5, 5, 'Jet Eta', 'Frequency', 'DelphesPythia8 EIC Event: Jet Eta', 'tab:blue')
error_bar_hist_plotting(jet_data_light['Jet.Eta'], 50, -5, 5, 'Jet Eta', 'Frequency', 'DelphesPythia8 EIC Event: Jet Eta', 'tab:red')
plt.legend(['All Jets','Charm Jets','Light Jets'],title='Jet Flavors')
plt.savefig('Graphics/jet_eta.jpg',dpi=300)
plt.show()
plt.clf()

# Analyzing Jet phi
error_bar_hist_plotting(jet_data['Jet.Phi'], 50, -math.pi, math.pi, 'Jet Phi', 'Frequency', 'DelphesPythia8 EIC Event: Jet Phi', 'black')
error_bar_hist_plotting(jet_data_charm['Jet.Phi'], 50, -math.pi, math.pi, 'Jet Phi', 'Frequency', 'DelphesPythia8 EIC Event: Jet Phi', 'tab:blue')
error_bar_hist_plotting(jet_data_light['Jet.Phi'], 50, -math.pi, math.pi, 'Jet Phi', 'Frequency', 'DelphesPythia8 EIC Event: Jet Phi', 'tab:red')
plt.legend(['All Jets','Charm Jets','Light Jets'],title='Jet Flavors')
plt.savefig('Graphics/jet_phi.jpg',dpi=300)
plt.show()
plt.clf()

# Analyzing Jet mass
error_bar_hist_plotting(jet_data['Jet.Mass'], 45, 0, 15, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'black')
error_bar_hist_plotting(jet_data_charm['Jet.Mass'], 45, 0, 15, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:blue')
error_bar_hist_plotting(jet_data_light['Jet.Mass'], 45, 0, 15, 'Jet Mass (GeV)', 'Frequency', 'DelphesPythia8 EIC Event: Jet Mass', 'tab:red')
plt.legend(['All Jets','Charm Jets','Light Jets'],title='Jet Flavors')
plt.grid()
plt.savefig('Graphics/jet_mass.jpg',dpi=300)
plt.show()
plt.clf()

# Analyzing Jet flavor
plt.hist(jet_data['Jet.Flavor'], bins=25, range=(0,25), histtype='step')
plt.xlabel('Jet Flavor')
plt.ylabel('Count')
plt.xlim(0,25)
plt.xticks((0.5,1.5,2.5,3.5,4.5,5.5,21.5),['$\gamma$','d','u','s','c','b','g'])
plt.title('DelphesPythia8 EIC Event: Jet Flavor')
plt.savefig('Graphics/jet_flavor.jpg',dpi=300)
plt.show()

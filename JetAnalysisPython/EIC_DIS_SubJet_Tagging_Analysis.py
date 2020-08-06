# In this file, we analyze ROOT files that are outputted by executing SimpleAnalysis/SubJetModule.cc
# The variables are all produced by the SubJetModule and deals primarily with sIP3D Tagging for Jets and SubJets

# Imports and Functions
import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib.ticker import MultipleLocator
import warnings
warnings.filterwarnings('ignore')

def error_calculate(data_array):
    raw_error = np.sqrt(data_array)
    rel_error = raw_error/data_array
    norm_data_array = data_array/sum(data_array)
    norm_rel_error = rel_error*norm_data_array
    norm_rel_error = np.nan_to_num(norm_rel_error)
    
    return norm_rel_error

def error_bar_hist_plotting(data, bin_count, range_min, range_max, xlabel, ylabel, title, color, error_bars_bool, log_scale_bool):
    plt.figure(0)
    (counts, bin_values, patches) = plt.hist(data, bins=bin_count, range=(range_min, range_max), histtype='step', color='white');
    counts = np.append(counts, 0)
    plt.clf()

    plt.figure(1)
    if(error_bars_bool == True):
        plt.errorbar(bin_values, counts/sum(counts), yerr=error_calculate(counts),ecolor='black',capsize=2,drawstyle='steps-mid',color=color)
    else:
        plt.errorbar(bin_values, counts/sum(counts),drawstyle='steps-mid',color=color)
    
    if(log_scale_bool == True):
        plt.yscale('log')
    
    plt.xlim(range_min, range_max)
    plt.ylim(0,max(max(counts/sum(counts))*1.05,plt.gca().get_ylim()[1]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.minorticks_on()

    return 0
    
#------------------------------------------------------------------------------------------------------------------------------------------------

# Open the ROOT file and place the variables into a DataFrame
file_data = uproot.open("D:\College/JetStudy/EICStudy/EIC_CCDIS_10x275_SubJetTagging_output-merged.root")["tree"]
tagged_data = file_data.pandas.df(['jet.P4','jet.flavor','jet.sIP3D','jet.TaggedsIP3D','subjet1.P4','subjet1.sIP3D','subjet1.TaggedsIP3D','subjet2.P4','subjet2.sIP3D','subjet2.TaggedsIP3D','subjet3.P4','subjet3.sIP3D','subjet3.TaggedsIP3D','subjet4.P4','subjet4.sIP3D','subjet4.TaggedsIP3D'], flatten=False)

#Arrays and Counter-Variables to use
charm_jet_count = []
charm_tagged_count = []
tagged_charm_jet_sIP3D_values = []
tagged_light_jet_sIP3D_values = []
untagged_charm_jet_sIP3D_values = []
untagged_light_jet_sIP3D_values = []
tagged_subjet_charm_sIP3D_values = []
tagged_subjet_light_sIP3D_values = []
untagged_subjet_charm_sIP3D_values = []
untagged_subjet_light_sIP3D_values = []
subjet_mass_tagged_charm = []
subjet_mass_total = []
tagged_subjet_charm_counter = 0
tagged_subjet_light_counter = 0
untagged_subjet_charm_counter = 0
untagged_subjet_light_counter = 0

#Finding how many jets per event (row) are charm jets and how many are tagged by the sIP3D algorithm
for i in range(len(tagged_data['jet.flavor'])):
    charm_jet_counter = 0
    charm_tagged_counter = 0
    
    for j in range(len(tagged_data['jet.flavor'][i])):
        if(tagged_data['jet.flavor'][i][j] == 4):
            charm_jet_counter+= 1
        if(tagged_data['jet.TaggedsIP3D'][i][j] == 1):
            charm_tagged_counter+= 1
                   
    charm_jet_count.append(charm_jet_counter)
    charm_tagged_count.append(charm_tagged_counter)
                   
tagged_data['n_charm_jet'] = charm_jet_count
tagged_data['n_charm_jet_TaggedsIP3D'] = charm_tagged_count

#Creating subset DataFrames that contain specific types of the data
truth_charm_jet_events = tagged_data[tagged_data['n_charm_jet'] > 0]
truth_charm_jet_events.reset_index(inplace=True, drop=True)

tagged_charm_jet_events = tagged_data[tagged_data['n_charm_jet_TaggedsIP3D'] > 0]
tagged_charm_jet_events.reset_index(inplace=True, drop=True)

true_tagged_charm_jet_events = tagged_data[(tagged_data['n_charm_jet'] > 0) & (tagged_data['n_charm_jet_TaggedsIP3D'] > 0)]
true_tagged_charm_jet_events.reset_index(inplace=True, drop=True)

#Simple printing out of the statistics from the events
print('Total Events:                 '+str(len(tagged_data['jet.P4'])))
print('True Charm Events:            '+str(len(truth_charm_jet_events['jet.P4'])))
print('Tagged Charm Events:          '+str(len(tagged_charm_jet_events['jet.P4'])))
print('True and Tagged Charm Events: '+str(len(true_tagged_charm_jet_events['jet.P4'])))

#Now we want to study the sIP3D values for certain jets, whether tagged or not and based on the truth-flavor
for i in range(len(tagged_data['jet.flavor'])):
    for j in range(len(tagged_data['jet.flavor'][i])):
        #Depending on the type of jet, put into an array the sIP3D values of the Jet
        if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['jet.TaggedsIP3D'][i][j] == 1):
            for k in range(len(tagged_data['jet.sIP3D'][i][j])):
                tagged_charm_jet_sIP3D_values.append(tagged_data['jet.sIP3D'][i][j][k])
        if(tagged_data['jet.flavor'][i][j] != 4 and tagged_data['jet.TaggedsIP3D'][i][j] == 1):
            for k in range(len(tagged_data['jet.sIP3D'][i][j])):
                tagged_light_jet_sIP3D_values.append(tagged_data['jet.sIP3D'][i][j][k])
        if(tagged_data['jet.flavor'][i][j] == 4 and tagged_data['jet.TaggedsIP3D'][i][j] == 0):
            for k in range(len(tagged_data['jet.sIP3D'][i][j])):
                untagged_charm_jet_sIP3D_values.append(tagged_data['jet.sIP3D'][i][j][k])
        if(tagged_data['jet.flavor'][i][j] != 4 and tagged_data['jet.TaggedsIP3D'][i][j] == 0):
            for k in range(len(tagged_data['jet.sIP3D'][i][j])):
                untagged_light_jet_sIP3D_values.append(tagged_data['jet.sIP3D'][i][j][k])
        
        #Depending on the type of jet, put into an array the sIP3D values of the SubJet (searching for if SubJets are tagged by sIP3D)
        #True Charm Jets, Looking for SubJets that pass sIP3D Tagging
        if(len(tagged_data['subjet1.TaggedsIP3D'][i][j]) != 0):
            subjet_mass_total.append(tagged_data['subjet1.P4'][i][j][3])
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet1.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_charm_counter+=1
                subjet_mass_tagged_charm.append(tagged_data['subjet1.P4'][i][j][3])
                for k in range(len(tagged_data['subjet1.sIP3D'][i][j])):
                    tagged_subjet_charm_sIP3D_values.append(tagged_data['subjet1.sIP3D'][i][j][k])
        if(len(tagged_data['subjet2.TaggedsIP3D'][i][j]) != 0):
            subjet_mass_total.append(tagged_data['subjet2.P4'][i][j][3])
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet2.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_charm_counter+=1
                subjet_mass_tagged_charm.append(tagged_data['subjet2.P4'][i][j][3])
                for k in range(len(tagged_data['subjet2.sIP3D'][i][j])):
                    tagged_subjet_charm_sIP3D_values.append(tagged_data['subjet2.sIP3D'][i][j][k])
        if(len(tagged_data['subjet3.TaggedsIP3D'][i][j]) != 0):
            subjet_mass_total.append(tagged_data['subjet3.P4'][i][j][3])
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet3.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_charm_counter+=1
                subjet_mass_tagged_charm.append(tagged_data['subjet3.P4'][i][j][3])
                for k in range(len(tagged_data['subjet3.sIP3D'][i][j])):
                    tagged_subjet_charm_sIP3D_values.append(tagged_data['subjet3.sIP3D'][i][j][k])
        if(len(tagged_data['subjet4.TaggedsIP3D'][i][j]) != 0):
            subjet_mass_total.append(tagged_data['subjet4.P4'][i][j][3])
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet4.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_charm_counter+=1
                subjet_mass_tagged_charm.append(tagged_data['subjet4.P4'][i][j][3])
                for k in range(len(tagged_data['subjet4.sIP3D'][i][j])):
                    tagged_subjet_charm_sIP3D_values.append(tagged_data['subjet4.sIP3D'][i][j][k])
                    
        #Not True Charm Jets, Looking for SubJets that pass sIP3D Tagging despite being not from a Charm Jet
        if(len(tagged_data['subjet1.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet1.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet1.sIP3D'][i][j])):
                    tagged_subjet_light_sIP3D_values.append(tagged_data['subjet1.sIP3D'][i][j][k])
        if(len(tagged_data['subjet2.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet2.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet2.sIP3D'][i][j])):
                    tagged_subjet_light_sIP3D_values.append(tagged_data['subjet2.sIP3D'][i][j][k])
        if(len(tagged_data['subjet3.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet3.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet3.sIP3D'][i][j])):
                    tagged_subjet_light_sIP3D_values.append(tagged_data['subjet3.sIP3D'][i][j][k])
        if(len(tagged_data['subjet4.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet4.TaggedsIP3D'][i][j][0] == 1):
                tagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet4.sIP3D'][i][j])):
                    tagged_subjet_light_sIP3D_values.append(tagged_data['subjet4.sIP3D'][i][j][k])

        #True Charm Jets, Looking for SubJets that fail sIP3D Tagging
        if(len(tagged_data['subjet1.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet1.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_charm_counter+=1
                for k in range(len(tagged_data['subjet1.sIP3D'][i][j])):
                    untagged_subjet_charm_sIP3D_values.append(tagged_data['subjet1.sIP3D'][i][j][k])
        if(len(tagged_data['subjet2.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet2.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_charm_counter+=1
                for k in range(len(tagged_data['subjet2.sIP3D'][i][j])):
                    untagged_subjet_charm_sIP3D_values.append(tagged_data['subjet2.sIP3D'][i][j][k])
        if(len(tagged_data['subjet3.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet3.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_charm_counter+=1
                for k in range(len(tagged_data['subjet3.sIP3D'][i][j])):
                    untagged_subjet_charm_sIP3D_values.append(tagged_data['subjet3.sIP3D'][i][j][k])
        if(len(tagged_data['subjet4.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] == 4  and tagged_data['subjet4.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_charm_counter+=1
                for k in range(len(tagged_data['subjet4.sIP3D'][i][j])):
                    untagged_subjet_charm_sIP3D_values.append(tagged_data['subjet4.sIP3D'][i][j][k])
                    
        #Not True Charm Jets, Looking for SubJets that fail sIP3D Tagging
        if(len(tagged_data['subjet1.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet1.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet1.sIP3D'][i][j])):
                    untagged_subjet_light_sIP3D_values.append(tagged_data['subjet1.sIP3D'][i][j][k])
        if(len(tagged_data['subjet2.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet2.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet2.sIP3D'][i][j])):
                    untagged_subjet_light_sIP3D_values.append(tagged_data['subjet2.sIP3D'][i][j][k])
        if(len(tagged_data['subjet3.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet3.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet3.sIP3D'][i][j])):
                    untagged_subjet_light_sIP3D_values.append(tagged_data['subjet3.sIP3D'][i][j][k])
        if(len(tagged_data['subjet4.TaggedsIP3D'][i][j]) != 0):
            if(tagged_data['jet.flavor'][i][j] != 4  and tagged_data['subjet4.TaggedsIP3D'][i][j][0] == 0):
                untagged_subjet_light_counter+=1
                for k in range(len(tagged_data['subjet4.sIP3D'][i][j])):
                    untagged_subjet_light_sIP3D_values.append(tagged_data['subjet4.sIP3D'][i][j][k])
                    
#Print out the statistics of SubJets in tagging and flavor information
print('------------------------------------------------------------------------')
print('Total # of Tagged SubJets in Truth Charm Jets  :',tagged_subjet_charm_counter)
print('Total # of Tagged SubJets in Truth Light Jets  :',tagged_subjet_light_counter)
print('Total # of Untagged SubJets in Truth Charm Jets:',untagged_subjet_charm_counter)
print('Total # of Untagged SubJets in Truth Light Jets:',untagged_subjet_light_counter)

print('Purity in Jets:'   ,((len(true_tagged_charm_jet_events['jet.P4']))/((len(tagged_charm_jet_events['jet.P4'])))))
print('Purity in SubJets:',tagged_subjet_charm_counter/(tagged_subjet_charm_counter+tagged_subjet_light_counter))

#Plotting Section - Plots on ANY of the arrays can be made, this is just a simple comparison of two of the arrays
error_bar_hist_plotting(tagged_charm_jet_sIP3D_values, 60, -20, 20, 'sIP3D', 'Frequency', 'DelphesPythia8 EIC Event: Jet sIP3D Values', 'tab:blue', False)
error_bar_hist_plotting(tagged_light_jet_sIP3D_values, 60, -20, 20, 'sIP3D', 'Frequency', 'DelphesPythia8 EIC Event: Jet sIP3D Values', 'tab:orange', False)
plt.legend(['Tagged+Charm','Tagged+Light'])
plt.vlines(3,0,1,linestyle='dashed',color='black') #3.0 is the sIP3D cut-off value, sIP3D above 3.0 contribute to whether a jet/subjet is tagged as true
plt.grid()
plt.show()

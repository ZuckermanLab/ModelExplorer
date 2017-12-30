#!/bin/python3
# Min Finder Script
# By August George 2017
# This script loads the evolver_rates data and finds the local minimuma below a threshold
import numpy as np
from numpy import genfromtxt
from scipy.signal import argrelmin
import pprint
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
from shutil import copyfile
import subprocess
import sys
#needs analyze model and evolver_rates in directory

#manually enter alpha, seed, n values
min_factor = 0.25 #(min value * min_factor = cutoff point for min finder )
#imported_alpha = 0.5
#imported_seed = 3456789
#imported_n = 1e5
#imported_analysis = "vary_dmu_S"

used with automation
imported_alpha = sys.argv[1]
imported_seed = sys.argv[2]
imported_n = sys.argv[3]
imported_analysis = sys.argv[4]

def graph_analysis_data(sub_imported_alpha,sub_imported_seed,sub_imported_n, sub_model_mcn, sub_analysis_type):
    analysis_alpha = sub_imported_alpha
    analysis_seed = sub_imported_seed
    analysis_nsteps = sub_imported_n
    analysis_mcn = sub_model_mcn
    analysis_type = sub_analysis_type

    if analysis_type == "vary_dmu_N":
        analysis_file_to_graph = "analysis-vary_dmu_N-dmu_init__-6__to__dmu_fin__0"
        varying_dmu = "Sodium"
    elif analysis_type == "vary_dmu_S":
        analysis_file_to_graph = "analysis-vary_dmu_S-dmu_init__-2__to__dmu_fin__4"
        varying_dmu = "Substrate"
    elif analysis_type == "vary_dmu_W":
        analysis_file_to_graph = "analysis-vary_dmu_W-dmu_init__-2__to__dmu_fin__4"
        varying_dmu = "Toxin"

    raw_analysis_data = np.genfromtxt(analysis_file_to_graph)
    analysis_x_dmu_n = raw_analysis_data[:,0]
    analysis_n = raw_analysis_data[:,1]
    analysis_s = raw_analysis_data[:,2]
    analysis_w = raw_analysis_data[:,3]
    analysis_y_sn = analysis_s / analysis_n
    analysis_y_wn = analysis_w / analysis_n
    analysis_y_ws = analysis_w / analysis_s


    fig = plt.figure()
    #ax = fig.add_subplot(111)
    big_title = "Proof-maker:  Sodium, Substrate, and Toxin Flux"
    title = "Reference Model: n = %s \nMC Parameters: Alpha = %s, Seed = %s, N-steps = %s" % (analysis_mcn, analysis_alpha, analysis_seed, analysis_nsteps)
    ylabel_text = "Flux (Out->In) [$mol$ $s^{-1}$ $m^{-2}$)]"
    xlabel_text = "$d\mu_{%s}$ (change in chemical potential of %s) [$k T$]" % (varying_dmu, varying_dmu)
    analysis_graph_name = "analysis_graph_a%s_s%s_mcn%s_modeln%s_vary_dmu_%s.png" % (analysis_alpha, analysis_seed, analysis_nsteps, analysis_mcn,varying_dmu)
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    fig.suptitle(big_title,fontweight = 'bold', fontsize = 12)
    plt.title(title, fontsize = 10)
    plt.grid(True)
    plt.ylabel(ylabel_text, fontweight = 'bold')
    plt.xlabel(xlabel_text, fontweight = 'bold')
    #ax.plot( analysis_x_dmu_n, analysis_y_sn, "r--", analysis_x_dmu_n, analysis_y_wn, "b--", analysis_x_dmu_n, analysis_y_ws, "g--" )
    plt.plot( analysis_x_dmu_n, analysis_y_sn, "r--", label="Substrate:Sodium")
    plt.plot( analysis_x_dmu_n, analysis_y_wn, "b", label="Toxin:Sodium")
    plt.plot( analysis_x_dmu_n, analysis_y_ws, "g--", label="Toxin:Substrate")
    #plt.legend(["Substrate:Sodium", "Toxin:Sodium", "Toxin:Substrate"])
    plt.legend()
    axes = plt.gca()
    axes.set_ylim([-2,2])
    fig.savefig(analysis_graph_name,dpi=300, bbox_inches='tight', )
    plt.clf()
    plt.close()

    sub_model_meta_data = "a%ss%smcn%s     %s     %s     %s     %s     %s\n" %(analysis_alpha, analysis_seed, analysis_mcn, analysis_type, analysis_x_dmu_n[0],analysis_y_sn[0],analysis_y_wn[0],analysis_y_ws[0])
    return sub_model_meta_data



#used with automation
#imported_alpha = sys.argv[1]
#imported_seed = sys.argv[2]
#imported_n = sys.argv[3]
#imported_analysis = sys.argv[4]
#compare_width = int(float(imported_n)*0.01)
#print("compare width = %s \n") % compare_width

#directory = "C:\Users\georgeau\Box Sync\August\ZuckermanLab\Proof-maker\RUNS\proof_maker_alpha0_s123456"
file_name = "evolver_rates.dat"
slice_initial = 500 # how many steps at the beginning to cut (eliminates noise and improves average calculation)

raw_mc_data = np.genfromtxt(file_name) #load data from file
copy_raw_data = raw_mc_data
print "\nraw mc data:\n"
pprint.pprint(raw_mc_data)



#plt.show()

### Finds some (not all) local minima from MC Energy Landscape
mc_energy_data = raw_mc_data[:,1] #gets only the second column (mc energies)
print "\nmc energy data:\n"
pprint.pprint(mc_energy_data)

#mc_energy_data[0:slice_initial] = 0 # data up to slice point is set to 0. This removes initial noise and helps with avereage calculation

mc_energy_mean  = np.mean(mc_energy_data)
print "\n mean = \n"
pprint.pprint(mc_energy_mean)

mc_energy_min = np.min(mc_energy_data)
print "\n min = \n"
pprint.pprint(mc_energy_min)

distance = np.sqrt((mc_energy_min-mc_energy_mean)**2)
print "\n distance = \n"
pprint.pprint(distance)


#min_threshold = (-0.05*distance) + mc_energy_mean
min_threshold = min_factor*mc_energy_min
print "\n min threshold = \n"
pprint.pprint(min_threshold)

below_threshold = mc_energy_data #copy data into new array
below_threshold[below_threshold > 0 ] = 0 #filter out positive values
below_threshold[below_threshold > min_threshold] = 0 #set values > threshold to zero
print "\n below threshold matrix\n"
pprint.pprint(below_threshold)

min_index = argrelmin(below_threshold) #finds min values
print "\nmin indexes:\n"
pprint.pprint(min_index)

min_value = mc_energy_data[min_index]
print "\nmin values:\n"
pprint.pprint(min_value)

markers_on = min_index[0][:].tolist()

#markers_on = min_index.tolist()

copy_raw_data =np.genfromtxt(file_name) #load data from file

x = copy_raw_data[:,0]
y = copy_raw_data[:,1]
ymax = np.amax(y)
ymin = np.amin(y)
#pprint.pprint(ymax)
#pprint.pprint(ymin)

fig = plt.figure()
ax = fig.add_subplot(111)
title = "Energy Fuction: -sflow*|sflow/wflow|^alpha\nMC Parameters: Alpha = %s, Seed = %s, N-steps = %s" % (imported_alpha, imported_seed, imported_n)
#ax.plot(x,y, label = 'MC Run')
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
fig.suptitle("Proof-maker: Monte Carlo Energy Landscape",fontweight = 'bold', fontsize = 12)
plt.title(title, fontsize = 10)
plt.grid(True)
plt.ylabel('MC Energy', fontweight = 'bold')
plt.xlabel('MC Step', fontweight = 'bold')
ax.plot(x,y,  linewidth=0.5, label = 'Approx. Local Min', marker='*', markerfacecolor='red', markeredgecolor='black', markevery=markers_on)
ax.legend()
mc_fig_title = "mc_energy_graph_a%s_s%s_mcn%s.png" % (imported_alpha, imported_seed, imported_n)
fig.savefig(mc_fig_title,dpi=300, bbox_inches='tight', )
plt.clf()
plt.close()
#plt.show()

meta_data_file = open("meta_data.dat", "a")
meta_header = "model     analysis     dmu N     S/N     W/N     W/S\n"
meta_data_file.write(meta_header)

#models_directory = "models_from_run/"
for array in min_index:
    for value in array:
        analysis_directory = "analysis_n%s" % value
        old_barriers_file = "models_from_run/barriers-%s" % value
        old_energies_file = "models_from_run/energies-%s" % value
        new_barriers_file = "analysis_n%s/barriers-%s" % (value, value)
        new_energies_file = "analysis_n%s/energies-%s" % (value, value)
        old_analysis_file = "analyze-model.prl"
        new_analysis_file = "analysis_n%s/%s" % (value, old_analysis_file)
        print(analysis_directory,old_barriers_file,new_barriers_file)
        if not os.path.exists(analysis_directory):
            os.makedirs(analysis_directory)
            copyfile(old_barriers_file, new_barriers_file)
            copyfile(old_energies_file, new_energies_file)
            copyfile(old_analysis_file, new_analysis_file)
            temp_file = "analysis_config.txt"
            config_text = "efile_init energies-%s \nbfile_init barriers-%s \nanalysis %s" % (value, value,imported_analysis)
            with open(os.path.join(analysis_directory, temp_file), 'wb') as config:
                config.write(config_text)
            pipe = subprocess.Popen(["perl", "analyze-model.prl"],cwd=analysis_directory)
            pipe.wait()
            current_directory = os.getcwd()
            print current_directory
            os.chdir(analysis_directory)
            model_meta_data = graph_analysis_data(imported_alpha,imported_seed,imported_n,value, imported_analysis)
            os.chdir(current_directory)
            meta_data_file.write(model_meta_data)

meta_data_file.close()

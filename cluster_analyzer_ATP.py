###  This script clusters models from a ModelExplorer simulation, and runs the
###  analysis perl script.
###  August George - 1/7/2018

import numpy as np
import sys
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.signal import argrelextrema, savgol_filter
from scipy import spatial
import networkx as nx
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import os
from shutil import copyfile
import subprocess
import csv
import datetime
import time
from decimal import Decimal
np.set_printoptions(threshold=sys.maxsize)


###  Imports a comma delimited float datafile and returns a numpy array
### 'invalid raise' skips line if an error is raised (i.e. unequal columns)
def import_data(datafile, proof=0):
    #print("importing data...\n")
    data = np.genfromtxt(datafile, delimiter=',', dtype=float, invalid_raise = False, skip_header=1)
    if proof ==1:
        assert np.shape(data)[1] == 45, "Data column size mismatch. Should be 45 (flows + n, emc, s,w,n flows) for proofreading."
    return data


### Slices array into three seperate arrays (1st col, 2nd col, remaining col).
### Returns col1 as a row, col2 as a row, remaning columns.
### This is useful for processing the imported cluster_data.dat file [mc_n, mc_e,...flows...]
def slice_array(array, run_length,offset,sample_rate, agg = False):
    """
    Slices array into 6 sections: MC step number, MC energy, state flows,
    N flow, S flow, and W flow.
    """
    #print("slicing array...\n")
    array_T = np.transpose(array)  # transpose to run col into row

    n = array_T[:1, :]  # first col of array
    emc = array_T[1:2, :]  # second col of array
    state_flows = array[:, 2:-3]  # remaining columns of array
    #flows = state_flows[:, :len(state_flows[0])/2]  # symmetry?
    N_flow = array[:, -3:-2]
    S_flow = array[:, -2:-1]
    W_flow = array[:, -1:]
    n = n.flatten()  # remove nesting
    emc = emc.flatten()

    # for aggregate 
    #run_length = 1e3
    #sample_rate = 1e1
    samples_per_run = run_length / sample_rate +1  # extra point sampled
    offset = 1/offset

    if agg:
        print("labeling for aggregation")
        for i in range(len(n)):
            n[i] = n[i] + (offset*np.floor(i/samples_per_run))  # add .0x to label data to run x
    
    return n, emc, state_flows, N_flow, S_flow, W_flow


###  Creates a heirarchical cluster using euclidian distance metric and single
###  nearest neighbor linkage. Creates a dendrogram based on cluster.
###  'truncate mode' is either 'lastp' merged clusters or p = 'level's
def create_cluster(matrix, minima):
    #print("creating a cluster...\n")
    dist = pdist(matrix, metric='euclidean')
    if dist.size == 0:
        print("distance matrix is empty. No clustering possible.")
        print("returning minima:%s" % (minima))
        return minima
    #square_dist_matrix = squareform(dist)
    Z = linkage(dist, 'complete')  # complete clustering
    fig = plt.figure(figsize=(24,8))
    with plt.rc_context({'lines.linewidth': 1}):
        dn = dendrogram(
            Z,
            truncate_mode = 'level',
            p = 6,
            color_threshold = 0.65,
            distance_sort = 'ascending',
            #no_labels = 'True',
            labels = minima,
            leaf_rotation = 90,
            leaf_font_size = 16,
            #orientation = 'right',
            #labels=np.asarray(minima).astype(int),
            show_contracted = 'True'
        )
    #dn = dendrogram(
        #Z,
        #distance_sort = 'ascending',
        #labels=minima,
    #)
    plt.title('Clustering hierarchy (truncated)\nUsing only filtered models\n', fontsize=32)
    plt.ylabel('Cluster distance [Euclidean]', fontsize=26)
    plt.xlabel('Model (cluster size)', fontsize=26)
    #plt.figtext(.80,.85,"Found %s Models" % len(minima), fontsize=30, ha='center')
    plt.tick_params(axis='both', which='major', labelsize=22)
    plt.margins(x=0.01, y=0.01)
    #plt.tick_params(axis='x', which='major', labelsize=20)
    #plt.xticks(fontsize=12)
    # plt.annotate("cluster A", xy = (0.10,0.30), xycoords = 'axes fraction', fontsize=22)
    # plt.annotate("cluster B", xy = (0.35,0.37), xycoords = 'axes fraction', fontsize=22)
    # plt.annotate("cluster C", xy = (0.60,0.23), xycoords = 'axes fraction', fontsize=22)
    # plt.annotate("cluster D", xy = (0.85,0.46), xycoords = 'axes fraction', fontsize=22)
    plt.savefig("model_heirarchy_cluster_complete.png", format='png', bbox_inches='tight')
    #plt.show()
    leaves = dn['ivl']  # get leaves of cluster
    #print(leaves)
    cleaned_leaves = [x for x in leaves if str(x).replace(
        ".", "", 1).isdigit()]  # remove truncated '(x)'
    #print(cleaned_leaves)
    #exit() # 
    return(cleaned_leaves,Z)


### Plots the energy trajectory data w/ optional minima markers
### Use marker indices for markers (not model numbers themselves)
### Note: can problably take out the optional markers since it is always passed markers
def graph_mc_data(x,y, markers={}):
    #print("graphing data...\n")
    fig = plt.figure(figsize=(50,30))
    ax = fig.add_subplot(111)
    #thresh = max(markers)
    if markers:  # for matrix use: markers.any()
        x2 = []
        y2 = []
        x2[:] = x[markers[:]]
        y2[:] = y[markers[:]]
        x2 = np.transpose(x2)
        y2 = np.array(y2)
        ax.plot(x2,y2, ls="", marker="o", label="Models below threshold", color='black')
        plt.figtext(.75, .80, "%s (of %s) Models below E_mc = %.1E" % (
            len(markers), len(x), (max(y2)+0.1*max(y2))), fontsize=30, ha='center')
    plt.title('Monte Carlo Energy Trajectory', fontsize=30)
    plt.ylabel('MC Energy', fontsize=26)
    plt.xlabel('MC Step', fontsize=26)
    ax.plot(x,y,  linewidth=2, color='orange')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    plt.tick_params(axis='both', which='major', labelsize=20)
    #for xy in zip(x2, y2):
    #    ax.annotate('(%s)' % xy[0], xy=xy, textcoords='data')
    plt.legend(fontsize=30)
    plt.margins(x=0.01, y=0.01)
    plt.ylim([-1e-4, 1e-4])
    plt.savefig("energy_trajectory_with_minima.png", format='png', bbox_inches='tight')
    plt.clf()


### Finds local min values of an array. Array is list (1D) of values.
### 'n' determines filter window length and number of neighbors to check
### for min-finding. 'thresh' filters out mins above certain value.
def find_min(array, n, thresh, filter=False):
    #print("finding local mins...\n")
    #print("raw number of models: %s" %np.size(array))
    min_index = []
    if filter:
        min_index = np.where(array < thresh)[0].tolist()
    else:
        mins = argrelextrema(array, np.less, order=n)
        mins = mins[0]

        for value in mins: # do in a more pythonic way
            if array[value] < thresh:
                min_index.append(value)

        if not min_index:
            min_i = np.argmin(array)
            min_v = array[min_i]
            print("no mins found...using argmin: i: %s a[i]: %s" % (min_i , min_v))
            min_index.append(np.argmin(array))

    return min_index


### Returns an array that has the flow values for each inputed index (minima)
def get_flows(array, indices):
    #print("getting flows...\n")
    a = array
    b = a[indices]
    return b


### Runs analysis script (ModelAnalyzer) for each local minima model in list.
def analyze_models(models_list):
    #print("analyzing models...\n\n")
    for value in models_list:
        analysis_directory = "analysis_n%s" %value
        old_barriers_file = "models_from_run/barriers-%s" % value
        old_energies_file = "models_from_run/energies-%s" % value
        new_barriers_file = "analysis_n%s/barriers-%s" % (value, value)
        new_energies_file = "analysis_n%s/energies-%s" % (value, value)
        old_analysis_file = "ModelAnalyzer.prl"
        new_analysis_file = "analysis_n%s/%s" % (value, old_analysis_file)
        old_cycle_check_file = "cycle_checker.py"
        new_cycle_check_file = "analysis_n%s/%s" % (value, old_cycle_check_file)
        old_analysis_config = "analysis_config.txt"
        new_analysis_config = "analysis_n%s/%s" % (value, old_analysis_config)
        old_pathways_file = "flux_grapher.py"
        new_pathways_file = "analysis_n%s/%s" % (value, old_pathways_file)
        config_text = "efile_init energies-%s \nbfile_init barriers-%s \ntimestamp %s\nanalysis vary_dmu_N\nanalysis_steps 1e1\n" % (value, value, str(datetime.datetime.now()))
        if not os.path.exists(analysis_directory):
            os.makedirs(analysis_directory)
            copyfile(old_barriers_file, new_barriers_file)
            copyfile(old_energies_file, new_energies_file)
            copyfile(old_analysis_file, new_analysis_file)
            copyfile(old_cycle_check_file, new_cycle_check_file)
            copyfile(old_analysis_config, new_analysis_config)
            copyfile(old_pathways_file, new_pathways_file)
            top_directory = os.getcwd()
            os.chdir(analysis_directory)
            with open(old_analysis_config, 'a') as file:
                file.write(config_text)
            pipe = subprocess.Popen(["perl", "ModelAnalyzer.prl"])
            pipe.wait()
            #print ("waiting...\n")
            os.chdir(top_directory)


### Takes an array of net flows (flux) and normalizes them (divides by max).
### Any flows below threshold are set to zero. Note that this assumes symmetry
### where flow A->B = - flow B->A
def normalize_and_threshold_flows(array):  # need to test
    #print("processing (normalizing and thresholding) flows...\n")
    a = array
    thresh = 0.01 # set anything below this to 0
    new_list = []
    for row in a:
        if max(row)>0:  # avoid dividing by zero (there shouldn't be a negative max due to symmetry)
            new_list.append([item*1./max(row) for item in row])  # divides each element by max in row
        else:
            new_list.append(row)
    b = np.array(new_list)
    b[b<thresh] = 0  # removes negative numbers (opposite direction) and small numbers
    #for (num,item) in enumerate(b):
        #print(num+1,item)
    return b

def normalize_flows(a): 
    '''
    Normalizes an array of flows
    '''
    b = np.linalg.norm(a, axis=1, keepdims=True)
    zero_rates = np.where(b == 0)[0]
    if zero_rates.size > 0:
        print("zeros in normalization matrix: %s" % zero_rates)
    a2 = np.delete(a, zero_rates, axis=0)
    b2 = np.linalg.norm(a2, axis=1, keepdims=True)
    zero_rates2 = np.where(b2 == 0)[0]
    a_norm = a2/b2   
    return a_norm




def indices_to_models(index_list, reference_array):
    """ Convert indices (minima) to models for analysis.
    Useful for when index doesn't correspond to MC n.
    (i.e. mc_n[0:2] = 0, 500, 1000)
    Return a list of model numbers (from ModelExplorer run)
    """
    #print("index to model...")
 
    models_list = []
    models_list = reference_array[index_list]
  
    models_list = [x for x in models_list]  # convert to int for use as index ?
    #for (num,item) in enumerate(models_list):
        #print(num+1,item)
    return models_list


def find_proofreading(s_flows, w_flows, n_flows, minima_indices, ddg=1.0, n=1, thresh = 1e-10):
    proofreading_models = []
    #print(len(minima_indices))
    
    for i in minima_indices:  # list of indices of models below threshold
        # print(i)
        # print(i*1e1)
        #print('1 checking for w > 1e-18, w < 1e+18')
        if np.abs(w_flows[i]) >= 1e-18 and np.abs(w_flows[i]) <= 1e+18:
            # print('!w > 1e-18, w < 1e+18!')
            #print(' 2 checking for abs(s_flows[i] / w_flows[i]) > (n*np.exp(ddg))')
            # print(n_flows[i])
            # print(s_flows[i])
            # print(w_flows[i])
            # print(abs(s_flows[i] / w_flows[i]))
            # print(n*np.exp(ddg))
            
            if ( (abs(s_flows[i] / w_flows[i]) > (n*np.exp(ddg)))):
                # print(' !abs(s_flows[i] / w_flows[i]) > (n*np.exp(ddg))!')
                #print('     3 checking s_flows[i] > thresh')
                
                if (s_flows[i] > thresh):
                    # print('     !s_flows[i] > thresh!')
                    #print('         4 checking n_flows[i] > 1e-12')
                    
                    if (n_flows[i] > 1e-12):
                        # print('         !n_flows[i] > 1e-12!')
                        #print('             5 checking abs(s_flows[i] / n_flows[i]) >= 0.1')
                        
                        if(abs(s_flows[i] / n_flows[i]) >= 0.1):
                            # print('             !abs(s_flows[i] / n_flows[i]) >= 0.1!')
                            #print('FOUND PROOFREADING MODEL')
                            # print(n_flows[i])
                            proofreading_models.append(i)
                            
    
    return proofreading_models


### Package all data processing (before clustering) into one subroutine
def process_data(datafile, graph=0, proof=0, threshold = 0, proof_thresh = 1e-10, run_length =1e3, offset = 1, sample_rate=1e1):
    print("processing data...\n")

    data = import_data(datafile, proof)  # cluster_data.dat is [mc_n, mc_e,...flows...]
    # mc_n, mc_e, flow_data, n_flow, s_flow, w_flow = slice_array(data[:int(run_length/d_print)])
    mc_n, mc_e, flow_data, n_flow, s_flow, w_flow = slice_array(data[:], run_length=run_length,offset=offset,sample_rate=sample_rate)
  
    n_models = len(mc_n)
    #thresh = np.max(mc_e)+1  # all models
    #thresh = -1e-7
    thresh = threshold
    minima_indices = find_min(mc_e, 1, thresh, filter=True)  # check n nearest neighbors below threshold
    n_thresh_models = len(minima_indices)
    thresh_ratio = 1.0 * n_thresh_models/n_models
    n_proof_models = 0
    proof_ratio = 0.0
    ddg_sw = 1.0
    proof_n = 10
    proof_t = proof_thresh

    if proof == 1:
        proof_thresh2 = 1e-05
        min_idx2 = find_proofreading(
            s_flow, w_flow, n_flow, minima_indices, ddg=ddg_sw, n=proof_n, thresh=proof_thresh2)

        minima_indices = find_proofreading(
            s_flow, w_flow, n_flow, minima_indices, ddg=ddg_sw, n =proof_n, thresh = proof_t )
        n_proof_models = len(minima_indices)
        proof_ratio = 1.0*n_proof_models/n_models


        proof_thresh3 = 1e-15
        min_idx3 = find_proofreading(
            s_flow, w_flow, n_flow, minima_indices, ddg=ddg_sw, n=proof_n, thresh=proof_thresh3)     
    minima_models_list = indices_to_models(minima_indices, mc_n)
  
    minima_models_list2 = indices_to_models(min_idx2, mc_n)
    minima_models_list3 = indices_to_models(min_idx3, mc_n)
    minima_flows = get_flows(flow_data, minima_indices)
    #old_processed_flows = normalize_and_threshold_flows(minima_flows)
    processed_flows = normalize_flows(minima_flows)
   
    #processed_flows = processed_flows[:int(2e6/500)]  # threshold on/off


    # def vector_autocorrelate(t_array):
    #     n_vectors = np.shape(t_array)[0]
    #     # correlate each component indipendently
    #     acorr = np.array([np.correlate(t_array[:, i], t_array[:, i], 'full')
    #                       for i in range(40)])[:, n_vectors-1:]
    #     # sum the correlations for each component
    #     acorr = np.sum(acorr, axis=0)
    #     # divide by the number of values actually measured and return
    #     acorr /= (n_vectors - np.arange(n_vectors))
    #     return acorr

    # # testing ACF for only good models at fixed min temp
    # # make sure looking at all models...
    # if np.shape(data)[0] == np.shape(processed_flows)[0]:
    #     print("looking at all samples\n")
    #     temp_period = 2000  # how many mc steps tempering cycle is
    #     initial_step_at_temp = 500  # first mc step at desired temp
    #     print_index = 500  # data printed every 500 mc steps
    #     good_models = range(int(500/500),np.shape(processed_flows)[0],int(2000/500))
    #     processed_flows2 = np.zeros_like(processed_flows)
    #     negative_emc_models = []
    #     for i, j in enumerate(good_models):
    #         processed_flows2[int(i)] = processed_flows[int(j)][:]
    #         negative_emc_models.append(mc_e[int(j)])
    #     processed_flows3 = np.delete(processed_flows2, np.s_[len(good_models):], axis=0)

    #     #print(len(good_models))
    #    # print(np.shape(processed_flows3)[0])
    #     #print(negative_emc_models)
    #     #print("selected models below mc e = 0: %s" % sum(n > 0 for n in negative_emc_models))
    #     correlation = acf(processed_flows3)
    #     # correlation = vector_autocorrelate(processed_flows)
    # else:
    #     correlation = acf(processed_flows)
    #     # correlation = vector_autocorrelate(processed_flows)
    

    # correlation = vector_autocorrelate(processed_flows)

    # plt.close()
    # ax1 = plt.subplot(111)
    # ax1.plot(mc_n,mc_e, linewidth=1)
    # ax1.set_ylim(-1e-3, 1e-3)
    # #ax1.set_yscale('log')
    # ax1.set_xlim(0, 3e6)
    # ax1.set_title("ModelExplorer trajectory in model space (truncated)")
    # ax1.set_ylabel("Monte Carlo energy")
    # ax1.set_xlabel("Monte Carlo iteration number [n]")
    # ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    # ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    # plt.tight_layout()
    # plt.show()
    #exit()

    #correlation = acf(processed_flows)
    #print("ACF[0] = %s" % correlation[0])
    if graph == 1:  # graph mc energy
        graph_mc_data(mc_n, mc_e, minima_indices)
        # fig = plt.figure(figsize=(25, 10))
        # ax = fig.add_subplot(211)
        # x = np.arange(correlation.shape[0])
        # #x = 500*x
        # ax.plot(x,correlation)
        # plt.minorticks_on()
        # #half_index = 500*(np.abs(correlation - 0.5)).argmin()
        # #zero_index = 500*(np.abs(correlation - 0.0)).argmin()
        # #ax.plot(half_index, correlation[int(half_index/500)], marker='o', markersize=7, color="red", label = half_index)
        # #ax.plot(zero_index, correlation[int(zero_index/500)], marker='o',
        #          #markersize=7, color="green", label=zero_index)
        # plt.legend()
        # #ax.hlines(y=0.5, xmin=0, xmax=len(correlation)-1, linestyles='dotted')
        # #ax.hlines(y=0, xmin=0, xmax=len(correlation)-1, linestyles='dotted')
        # plt.title(
        #     "Autocorrelation of MC Run\n(MCn=5e6, seed=456789, demax=0.6, b_low=125, b_hi=1450, b_inc=100, b_dec=325)", size=14)
        # plt.ylabel("Autocorrelation", size=12)
        # #plt.xlabel("Lag time", size=12)
        # plt.xlim(0,10)
        # plt.ylim(-1,1)
        # #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        # ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        # ax.tick_params(axis='both', which='major', labelsize=10)
        # ax.tick_params(axis='both', which='minor', labelsize=8)

        # plt.savefig("mc_auto_correlation.png", format='png', bbox_inches='tight')
        fig = plt.figure()
        plt.plot(minima_models_list2, np.arange(len(minima_models_list2)),
                label="abs(s/w flow) > %s*exp(%s)) & s flow > %s" % (proof_n, ddg_sw, proof_thresh2))
        plt.plot(minima_models_list, np.arange(len(minima_models_list)),
                label="abs(s/w flow) > %s*exp(%s)) & s flow > %s" % (proof_n, ddg_sw, proof_thresh))

        plt.plot(minima_models_list3, np.arange(len(minima_models_list3)),
                label="abs(s/w flow) > %s*exp(%s)) & s flow > %s" % (proof_n, ddg_sw, proof_thresh3))
        plt.ylabel("Number of 'good' models")
        plt.xlabel("MC step number")
        plt.title("Number of 'good' models during an MC run")
        plt.legend(loc='best')
        #plt.show()
        plt.savefig("good_models_mc_run.png",
                    format='png', bbox_inches='tight')

    if proof == 0:
        report = "Total Models: %s\nModels below MC energy threshold (%s) : %s, Ratio: %s\n" % (
            n_models, thresh, n_thresh_models, thresh_ratio)
    else:
        report = "Total Models: %s\nModels below MC energy threshold (%s) : %s, Ratio: %s\nModels above proofreading threshold %s*(e^%s): %s, Ratio: %s" % (
            n_models, thresh, n_thresh_models, thresh_ratio, proof_n, ddg_sw, n_proof_models, proof_ratio)
    return processed_flows, minima_models_list, report


def model_correlation(norm_flows, model_list, thresh, max_d=3e6):
    dist_mat = squareform(pdist(norm_flows, metric='euclidean'))
    #print(dist_mat.shape)
    n = int(dist_mat.shape[0])
    avg = []
    std_err = []
    count = []

    if not model_list:
        print("no models found. avg = 2.5e6 +/- 0")
        avg = [max_d] * len(thresh)
        std_err = [0.0] * len(thresh)
        count = [0] * len(thresh)
    else:
        for t in thresh:
            d = []
            for i in np.arange(n):
                # find first element over threshold
                idx = np.argmax(dist_mat[i, i:] > t)
                if idx == 0 or np.isnan(idx):  # if there is no model > thresh
                    d.append(max_d-model_list[i])  # distance = distance to last model
                else:
                    d.append(model_list[int(idx + i)]-model_list[i])

            # if not d:  # if empty (no models)
            #     d.append(0.0)
            #     count.append(0)
            # else:
            #     count.append(int(len(d)))
            count.append(int(len(d)))
            dist = np.asarray(d)
            avg.append(float(np.average(dist)))
            std_err.append(float(np.std(dist)/np.sqrt(dist.shape[0])))

    return avg, std_err, count

### Compares each consecutive model in dataset1 with models in dataset2
### to find the minimum distance
def compare_two_runs(datafile1, datafile2):
    #print("comparing %s to %s...\n" % (datafile1, datafile2))
    processed_flows1, minima_models_list1, report1 = process_data(datafile1, proof=1)
    processed_flows2, minima_models_list2, report2 = process_data(datafile2, proof=1)
    # cdist creates a distance matrix between all elements in two matrices
    # np.min on axis 1 finds the min for for each col (each model in run1)
    mindist = np.min(spatial.distance.cdist(processed_flows1, processed_flows2), axis=1)
    return mindist, minima_models_list1, report1, report2


### Plot 2 histograms: Run1 models distance from Run2 models and vice versa
### Should clean up later (don't repeat yourself!)
def plot_two_histograms(run1_data, run2_data):
    #print("graphing histograms...\n")
    fig = plt.figure(figsize=(30, 10))
    plt.suptitle(
        "Comparing Models of Two Runs\n(demax = 0.6, temp scale = 1, mc n = 5e6) vs (demax = 1, temps scale = 1, mc n = 3e6)")
    

    ax1 = fig.add_subplot('121')  # top left histogram (run1)
    ax1.tick_params(axis='x', which='major', labelsize=10)
    ax1.title.set_text('Run1 vs Run2')
    ax1.set_ylabel('Frequency')
    ax1.set_xlabel('Minimum Distance (run1 model to run2 model)')
    bin_size1 = int(round(len(run1_data)**(0.5)))  # quick estimate of bin size is sqrt(n)
    bin_step1 = (max(run1_data) - min(run1_data))/bin_size1
    if bin_step1 == 0:
        bin_step1 = 1
    ax1.hist(run1_data, bins=bin_size1)
    ax1.set_xticks(np.arange(min(run1_data), max(run1_data)+bin_step1, bin_step1))  # use arrange instead of range to handle floats
    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax1.tick_params(axis='y', which='minor', bottom=False)
    plt.figtext(.40,.85,"Run 1 has %s Models" % len(run1_data))

    ax2 = fig.add_subplot('122')  # top right histogram (run2)
    ax2.tick_params(axis='x', which='major', labelsize=10)
    ax2.title.set_text('Run2 vs Run1')
    ax2.set_ylabel('Frequency')
    ax2.set_xlabel('Minimum Distance (run2 model to run1 model)')
    bin_size2 = int(round(len(run2_data)**(0.5)))  # quick estimate of bin size is sqrt(n)

    bin_step2 = (max(run2_data) - min(run2_data))/bin_size2
    if bin_step2 == 0:
        bin_step2 = 1
    ax2.hist(run2_data, bins=bin_size2)
    ax2.set_xticks(np.arange(min(run2_data), max(run2_data)+bin_step2, bin_step2))  # use arrange instead of range to handle floats
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax2.tick_params(axis='y', which='minor', bottom=False)
    plt.figtext(.825,.85,"Run 2 has %s Models" % len(run2_data))
    #plt.show()
    plt.savefig("run_comparison_histogram.png", format='png', bbox_inches='tight')


### finds models that are above (or below) a thresholdself.
### need to fix to be more pythonic (remove for loop)
def find_outlier_models(mindist, model_list, thresh, below=False):
    outlier_models = []
    if below:
        for i, x in enumerate(mindist):
            if x < thresh:
                #print("i:%s, x:%s, model:%s" % (i, x, model_list[i]))
                outlier_models.append(model_list[i])
    else:
        for i, x in enumerate(mindist):
            if x > thresh:
            #print("i:%s, x:%s, model:%s" % (i, x, model_list[i]))
                outlier_models.append(model_list[i])
    #print(outlier_models)
    return list(outlier_models)


def acf(data):
    '''
    Calculates the normalized Autocorrelation function (ACF)
    Parameters:
        data: ndarray
            n x m matrix (n of vectors) of data
    Returns:
        out: ndarray
            ACF of data
    Notes: CF_k = <(x(t)-x.mean).(x(t+k)-x.mean)> / <(x(t)-x.mean).(x(t)-x.mean)>
    where k is the lag time, and < > denotes the average over t. 
    '''

    N = np.shape(data)[0]
    mean = np.mean(data)
    acf = np.zeros(N)

    # calculate den = <(x(i)-x.mean).(x(i)-x.mean)>
    den = 0.0
    for t in range(N):  # t = 0 to N-1
        den = den + np.dot((data[t]-mean), (data[t]-mean))
    den = den/N

    # calculate num = <(x(t)-x.mean).(x(t+k)-x.mean)>
    for k in range(N):  # k = 0 to N-1 (lag time)
        num = 0.0
        for t in range(N-k):  # for t = 0 to N-k-1 (window)
            num = num + np.dot((data[t]-mean), (data[t+k]-mean)) 
        num = num/(N-k)
        acf[k] = num/den   

    return acf


def avg_cluster_change_time(clusters, models_list, max = 2.5e6):
    for i in np.arange(clusters):
        pass
        # move until cluster changes
        # calculate time difference
        # repeat until end of cluster list




### Main program follows simple pipeline architecture:
### 1) import raw data, 2) seperate data into mc_n, mc_e, and flows,
### 3) find local minima from mc_n, 4) get net flows of minima,
### 5) normalize net flows of minima, 6) graph mc_n and mc_e w/ minima
### 7) cluster models based on normalized net flows of energy minima models
### 8) run analysis script for each energy minima model
def cluster_and_analyze_one_run(datafile, analyze=None, graph=1, proof=0, threshold=1e+100, run_length=1e6, n_runs = 1, sample_rate=1e1):
        print("clustering and analyzing one run: %s\n" % datafile)
        start_time = time.time()  # testing for runtime data
        print("processing data...\n")
        offset = 1*10**(int(np.log10(n_runs)+1))
        print("offest = %s" %offset)
        processed_flows, minima_models_list, report = process_data(
            datafile, graph=graph, proof=proof, threshold=threshold, run_length=run_length, offset=offset, sample_rate=sample_rate)
        n_steps = run_length
        agg_data_analysis(minima_models_list, offset = offset, n_steps=n_steps, n_runs=n_runs)    
        print("clustering data...\n")
        cluster_models = create_cluster(processed_flows, minima_models_list)
        #print(cluster_models)
        # array M_i = cluster number for sample at i
        clustering = fcluster(cluster_models[1], t=0.65, criterion='distance')
        #print(minima_models_list)
        #print(clustering)
        #print(np.max(clustering))
        #print(np.max(np.asarray(minima_models_list)))
        #print(np.bincount(clustering))
        #print(clustering.shape)
        plt.close()
        ax1 = plt.subplot(111)
        ax1.step(np.asarray(minima_models_list), clustering, where='post', linewidth=1.5)
        ax1.set_ylim(np.min(clustering), 8)
        ax1.set_xlim(0,3e6)
        ax1.set_title("ModelExplorer cluster trajectory ")
        ax1.set_ylabel("Cluster index [i]")
        ax1.set_xlabel("Monte Carlo iteration number [n]")
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        plt.tight_layout()
        #plt.show()
        
        if analyze == 'Run1':
            print("analyzing models...\n")
            analyze_models(cluster_models[0])
        runtime = time.time()-start_time  # testing for runtime data
        print(report)
        #print(cluster_models[0])
        if analyze == 'Run1':
            print("Models analyzed: %s" % (len(cluster_models[0])))
        print("Runtime: %s (s)" % (runtime))
        print(time.strftime("%Y-%m-%d %H:%M"))

        return (minima_models_list, processed_flows)


def compare_and_analyze_two_runs(datafile1, datafile2, thresh, analyze=None):
    start_time = time.time()  # testing for runtime data
    print("Comparing two runs: %s and %s" % (datafile1, datafile2))
    print("Min distance threshold = %s" % thresh)
    print("Analyze (None if blank): %s \n" % analyze)
    r1_vs_r2, r1_labels, report1, report2 = compare_two_runs(datafile1, datafile2)
    r2_vs_r1, r2_labels, _, _ = compare_two_runs(datafile2, datafile1)
    plot_two_histograms(r1_vs_r2, r2_vs_r1)
    exit()

    outliers1 = find_outlier_models(r1_vs_r2, r1_labels, thresh, below=True)
    outliers2 = find_outlier_models(r2_vs_r1, r2_labels, thresh, below=True)
    if analyze == 'Run1':
        analyze_models(outliers1)
    elif analyze == 'Run2':
        analyze_models(outliers2)
    runtime = time.time()-start_time
    print("Run1 had %s Models w/ %s outliers %s: %s" % (len(r1_labels), len(outliers1), thresh, outliers1))
    print(report1)
    print("Run2 had %s Models w/ %s outliers %s: %s" % (len(r2_labels), len(outliers2), thresh, outliers2))
    print(report2)
    print("Total Runtime: %s (s)" % runtime)
    print(time.strftime("%Y-%m-%d %H:%M"))

def compare_n_runs(run_list, thresh, flow_thresh = None):
    start_time = time.time()
    flow_thresh = [1e-10]
    fig = plt.figure(figsize=(30, 10))
    for f in flow_thresh:
        d = []
        for run in run_list:
            processed_flows, minima_models_list, report = process_data(
                run[0], graph=0, proof=1, threshold=1e100, proof_thresh=f)
            avg, std_err, count = model_correlation(processed_flows, minima_models_list, thresh)
            d.append((avg,std_err,run[1],count))
        print("demax optimization")
        for t in range(len(thresh)):
            print("threshold: %s" % thresh[t])
            x = []
            y = []
            y_err = []
            for i in range(len(run_list)):
                print("  demax: %s. avg: %s +/- %s. %s models" %
                    (d[i][2], d[i][0][t], 2*d[i][1][t], d[i][3][t]))
                if d[i][0][t] != 0:
                    x.append(d[i][2])
                    y.append(d[i][0][t])
                    y_err.append(2*d[i][1][t])
            # fmt='none',
            plt.errorbar(x, y, y_err,  label="dist > %s & s flow > %s" %
                        (thresh[t],f), linewidth=2.0)
    plt.title(
        "Demax and S flow threshold optimization (temp. scale = 1)\nusing models where abs(s/w flow) > 10*exp(1)) & s flow > x")
    plt.ylabel("Avg. MC steps to model above dist. thresh.")
    plt.xlabel("Demax")
    plt.legend()
    plt.yscale('log')
    plt.grid(b=True, which='both', axis='both')
    plt.tight_layout()
    #plt.show()
    runtime = time.time()-start_time
    print("Total Runtime: %s (s)" % runtime)
    print(time.strftime("%Y-%m-%d %H:%M"))
    

def aggregate_analysis():

   #D_list = ["cluster_data_d02.dat", "cluster_data_d08.dat"]

    run_size = 1e4
    step_size = 1e1

    d1 = import_data("cluster_data_run0.dat", proof=1)  # demax=1 temp_scale=1
    d2 = import_data("cluster_data_run1.dat",
                     proof=1)  # demax=0.5 temp_scale=0.5
    d3 = import_data("cluster_data_run2.dat", proof=1)  # demax=0.2 temp_scale=1
    d4 = import_data("cluster_data_run3.dat",
                     proof=1)  # demax=1 temp_scale=2

    # d1 = import_data("cluster_data_d1_ts1_s123456.dat", proof=1)
    # d2 = import_data("cluster_data_d1_ts1_s234567.dat", proof=1)
    # d3 = import_data("cluster_data_d1_ts1_s345678.dat", proof=1)
    # d4 = import_data("cluster_data_d1_ts1_s456789.dat", proof=1)


    D1 = d1[:int(run_size/step_size)]
    D2 = d2[:int(run_size/step_size)]
    D3 = d3[:int(run_size/step_size)]
    D4 = d4[:int(run_size/step_size)]



    step_size = 1e1
    D1[:, 0] = np.arange(0, D1.shape[0]*step_size, step_size) + 0.1
    D2[:, 0] = np.arange(0, D2.shape[0]*step_size, step_size) + 0.2
    D3[:, 0] = np.arange(0, D3.shape[0]*step_size, step_size) + 0.3
    D4[:, 0] = np.arange(0, D4.shape[0]*step_size, step_size) + 0.4

    D = np.concatenate(
        (D1, D2, D3, D4), axis=0)
    # D = np.concatenate(
    #     (D1[:int(1.2e6/step_size)], D2[:int(1.2e6/step_size)]), axis=0)

    
    #print(D[:, 0])

    # cluster_data.dat is [mc_n, mc_e,...flows...]
    mc_n, mc_e, flow_data, n_flow, s_flow, w_flow = slice_array(D)
    n_models = len(mc_n)
    #thresh = np.max(mc_e)+1  # all models
    #thresh = -1e-7
    thresh = 1e100
    # check n nearest neighbors below threshold
    minima_indices = find_min(mc_e, 1, thresh, filter=True)
    n_thresh_models = len(minima_indices)
    thresh_ratio = 1.0 * n_thresh_models/n_models
    n_proof_models = 0
    proof_ratio = 0.0
    ddg_sw = 1.0
    proof_n = 10
    proof_t = 1e-10

    minima_indices = find_proofreading(s_flow, w_flow, n_flow, minima_indices, ddg=ddg_sw, n=proof_n, thresh=proof_t)
    #print(minima_indices)
    n_proof_models = len(minima_indices)
    proof_ratio = 1.0*n_proof_models/n_models
    minima_models_list = indices_to_models(minima_indices, mc_n)
    print(minima_models_list)
    


    #exit()
    #print(minima_models_list[0:])
    
    minima_flows = get_flows(flow_data, minima_indices)
    #old_processed_flows = normalize_and_threshold_flows(minima_flows)
    processed_flows = normalize_flows(minima_flows)

    report = "Total Models: %s\nModels below MC energy threshold (%s) : %s, Ratio: %s\nModels above proofreading threshold %s*(e^%s) & s flow > %s: %s, Ratio: %s" % (
        n_models, thresh, n_thresh_models, thresh_ratio, proof_n, ddg_sw, proof_t, n_proof_models, proof_ratio)
    print(report)

    cluster_models = create_cluster(processed_flows, minima_models_list)
    clustering = fcluster(cluster_models[1], t=0.65, criterion='distance')  # array M_i = cluster number for sample at i
    print("%s clusters using euclidean distance threshold = 0.65" % np.max(clustering))
    print(np.bincount(clustering))
    print(clustering.shape)


    # fix later (don't repeat yourself...)
    d1_idx = np.where(np.logical_and(np.greater_equal(
        np.asarray(minima_indices), 0), np.less_equal(np.asarray(minima_indices), D1.shape[0])))[0].astype(int)
    d1_x = np.asarray(minima_models_list[d1_idx[0]:d1_idx[-1]+1])-0.1
   
    d2_idx = np.where(np.logical_and(np.greater_equal(
        np.asarray(minima_indices), D1.shape[0]), np.less_equal(np.asarray(minima_indices), D1.shape[0]+D2.shape[0])))[0].astype(int)
    d2_x = np.asarray(minima_models_list[d2_idx[0]:d2_idx[-1]+1])-0.2

    d3_idx = np.where(np.logical_and(np.greater_equal(
        np.asarray(minima_indices), D1.shape[0]+D2.shape[0]), np.less_equal(np.asarray(minima_indices), D1.shape[0]+D2.shape[0]+D3.shape[0])))[0].astype(int)
    d3_x = np.asarray(minima_models_list[d3_idx[0]:d3_idx[-1]+1])-0.3

    d4_idx = np.where(np.logical_and(np.greater_equal(
        np.asarray(minima_indices), D1.shape[0]+D2.shape[0]+D3.shape[0]), 
        np.less_equal(np.asarray(minima_indices), D1.shape[0]+D2.shape[0]+D3.shape[0]+D4.shape[0])))[0].astype(int)
    d4_x = np.asarray(minima_models_list[d4_idx[0]:d4_idx[-1]+1])-0.4

    d1_y = clustering[0:d1_x.shape[0]]
    d2_y = clustering[d1_x.shape[0]:(d1_x.shape[0]+d2_x.shape[0])]
    d3_y = clustering[(d1_x.shape[0]+d2_x.shape[0]):(d1_x.shape[0]+d2_x.shape[0]+d3_x.shape[0])]
    d4_y = clustering[(d1_x.shape[0]+d2_x.shape[0]+d3_x.shape[0]):(d1_x.shape[0]+d2_x.shape[0]+d3_x.shape[0]+d4_x.shape[0])]

    x_s = [d1_x,d2_x,d3_x,d4_x]
    y_s = [d1_y, d2_y, d3_y, d4_y]

    #x_s = [d1_x, d2_x]
    #y_s = [d1_y, d2_y]

    #print(run2_cluster_traj)
    plt.close()

    ### KEEP
    plt.suptitle("ModelExplorer trajectory in model space", fontsize=36)

    for i in range(0,4):
        
        x = mc_n[int(i*run_size/step_size):int((i+1)*run_size/step_size)]
        y  = mc_e[int(i*run_size/step_size):int((i+1)*run_size/step_size)]
        plot_n = 221+i
        ax1 = plt.subplot(plot_n)
        ax1.plot(x, y, linewidth=0.5)
        ax1.set_ylim(-1e-3, 1e-3)
        ax1.set_xlim(0,run_size)
        ax1.set_title("Run %s" % (i+1), fontsize=24)
        if i == 0:
            ax1.set_ylabel(
                "Monte Carlo 'energy' [arb. flux units]", fontsize=20)
        if i == 2:
            ax1.set_ylabel("Monte Carlo 'energy' [arb. flux units]", fontsize = 20)
            ax1.set_xlabel("Monte Carlo iteration number [n]", fontsize=20)
        if i == 3:
            ax1.set_xlabel("Monte Carlo iteration number [n]", fontsize=20)
        #plt.yticks(np.arange(0, 13, 1))
        plt.xticks(np.arange(0, run_size*1.1, 1e5))
        ax1.tick_params(axis='both', which='major', labelsize=14)
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    #plt.tight_layout()
    #plt.show()

    plt.close()

    ### KEEP
    plt.suptitle("ModelExplorer trajectory in cluster space", fontsize=32)

    for i in range(0, 4):
        print(np.bincount(y_s[i]))
        x = np.append(x_s[i], run_size)
        y = np.append(y_s[i], y_s[i][-1])
        plot_n = 221+i
        ax1 = plt.subplot(plot_n)
        #ax1.step(x, y, where='post', linewidth=1.5)
        ax1.plot(x,y)
        ax1.set_ylim(0, max(clustering))
        ax1.set_xlim(0, run_size)
        ax1.set_title("Run %s" % (i+1), fontsize=24)
        if i == 0:
            ax1.set_ylabel("Cluster index [i]", fontsize=20)
        if i == 2:
            ax1.set_ylabel("Cluster index [i]", fontsize=20)
            ax1.set_xlabel("Monte Carlo iteration number [n]", fontsize=20)
        if i == 3:
            ax1.set_xlabel("Monte Carlo iteration number [n]", fontsize=20)

        ax1.tick_params(axis='both', which='major', labelsize=14)
        plt.yticks(np.arange(0, max(clustering)+2, 1))
        plt.xticks(np.arange(0, 1.1e6, 1e5))
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    #plt.tight_layout()
    #plt.show()



    # plt.suptitle("ModelExplorer cluster distributions", fontsize=20)
    # for i in range(0, 4):
    #     #print(np.bincount(y_s[i]))
    #     plot_n = 331+i
    #     ax1 = plt.subplot(plot_n)
    #     ax1.bar(np.arange(np.bincount(y_s[i]).shape[0]), np.bincount(y_s[i]))
    #     #ax1.set_ylim(0, np.max(clustering))
    #     ax1.set_xlim(0, 16)
    #     ax1.set_title("Run %s" % (i+1))
    #     ax1.set_ylabel("Counts for cluster i [n]")
    #     ax1.set_xlabel("Cluster index [i]")
    #     plt.xticks(np.arange(0, 16, 1.0))
    #     #ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    # plt.tight_layout()
    # plt.show()


def agg_data_analysis(minima_models_list, n_runs = 6*50, offset = 1000, n_steps = 1e3):
    print("aggregate analysis: %s runs, offset = %s, n steps: %s\n" %(n_runs, offset,n_steps))

    def avg(list_of_times):
        return sum(list_of_times)/len(list_of_times)

    # models in run
    steps_to_min = []
    avg_per_run = []
    min_per_run = []
    run_n_0 = 0
    # histogram of runs
    run_hist = np.zeros(n_runs)
    for model in minima_models_list:
        model_n = int(offset*np.floor(model))/(offset)
        run_n = int(np.round(model*offset)) - int(offset*np.floor(model))  # get run number
        #print(run_n)
        run_hist[run_n] = run_hist[run_n] + 1  # add count to histogram
        if run_n == run_n_0:
            steps_to_min.append(model_n)
        else:
            try:
                avg_per_run.append(sum(steps_to_min)/len(steps_to_min))
                min_per_run.append(min(steps_to_min))
            except:
                pass         
            run_n_0 = run_n
            steps_to_min.clear()
            steps_to_min.append(model_n)
    try:
        avg_per_run.append(sum(steps_to_min)/len(steps_to_min)) # store last one
        min_per_run.append(min(steps_to_min))
        print("min time to find filtered model (using best run): %s" % min(min_per_run))

        print("avg MIN mc n to find filtered models (using only good runs): %s" %avg(min_per_run))

        min_per_run.extend([n_steps for i in range(int(n_runs - len(min_per_run)))])
        print("avg MIN mc n to find filtered models (using ALL runs): %s" %avg(min_per_run))

        print("avg mc n to find filtered models (using only good runs): %s" %avg(avg_per_run))

        avg_per_run.extend([n_steps for i in range(int(n_runs - len(avg_per_run)))])
        print("avg mc n to find filtered models (using ALL runs): %s" %avg(avg_per_run))
    except:
        print("no models found")
  
    #print(run_hist)
    print("highly fit models in %s/%s runs" % (np.count_nonzero(run_hist),n_runs))
    print(run_hist)
    
    if np.count_nonzero(run_hist) != 0:
        fig = plt.figure(figsize=(30, 10))
        ax1 = fig.add_subplot('111') 
        ax1.set_title("Distribution of highly fit models\n%s/%s runs have highly fit models" % (np.count_nonzero(run_hist),n_runs))
        ax1.set_ylabel('Count')
        ax1.set_xlabel('Run number')
        ax1.bar(range(n_runs),run_hist, align='edge')
        plt.savefig("filtered_models_distribution.png", format='png', bbox_inches='tight')

    


def analyze_experiments():
    # analyze/cluster aggregate data from several experiments

    print("analyzing experiments...aggregate")
    # initialize data
    print("experiment 1")
    (min_list1, flows1) = cluster_and_analyze_one_run("cluster_data_agg1.dat", analyze='None',
                              graph=0, proof=1, threshold=1e+100, run_length=5e3, n_runs = 6*50, sample_rate=1e1)
    exp1= [str(x) + '.1' for x in min_list1]

    print("experiment 2")
    (min_list2, flows2) = cluster_and_analyze_one_run("cluster_data_agg2.dat", analyze='None',
                              graph=0, proof=1, threshold=1e+100, run_length=1e3, n_runs = 6*50, sample_rate=1e1)
    exp2= [str(x) + '.2' for x in min_list2]

    print("experiment 3")
    (min_list3, flows3) = cluster_and_analyze_one_run("cluster_data_agg3.dat", analyze='None',
                              graph=0, proof=1, threshold=1e+100, run_length=1e4, n_runs = 6*50, sample_rate=1e1)
    exp3= [str(x) + '.3' for x in min_list3]

    print("experiment 4")
    (min_list4, flows4) = cluster_and_analyze_one_run("cluster_data_agg4.dat", analyze='None',
                              graph=0, proof=1, threshold=1e+100, run_length=5e2, n_runs = 6*50, sample_rate=1e1)
    exp4= [str(x) + '.4' for x in min_list4]

    flows_agg = np.vstack([flows1,flows2, flows3,flows4])
    exp_label_agg = exp1 + exp2 + exp3 + exp4

    print("aggregate data...")
    cluster_models = create_cluster(flows_agg, exp_label_agg)
    clustering = fcluster(cluster_models[1], t=0.65, criterion='distance')
    #exit()
    

def plot_coupling(H_flow, ATP_flow, n_runs, run_length, sample_rate):
    ### prints the mc_e for each run
    fig = plt.figure(figsize=(9, 3))
    ax1 = fig.add_subplot('111') 
    ax1.set_title("Coupling ratios (n = %s)" %(n_runs))
    ax1.set_ylabel('Coupling ratio (H+/ATP)')
    ax1.set_xlabel('MC step number')
    x = np.arange(0,run_length + sample_rate,sample_rate)
   
    # for each run 
    for i in range(n_runs):
        start_idx = i*int(run_length/sample_rate + 1)
        end_idx = start_idx + int(run_length/sample_rate + 1)
        y = (H_flow[start_idx:end_idx])/(ATP_flow[start_idx:end_idx])
        ax1.plot(x,y, alpha=0.7)
    #plt.yscale('symlog')
    plt.savefig("coupling_ratio.png", format='png', bbox_inches='tight')


   

def plot_mc_c(mc_e, n_runs, run_length, sample_rate):
    ### prints the mc_e for each run
    fig = plt.figure(figsize=(9, 3))
    ax1 = fig.add_subplot('111') 
    ax1.set_title("MC Energy Trajectories (n = %s)" %(n_runs))
    ax1.set_ylabel('MC energy')
    ax1.set_xlabel('MC step number')
    x = np.arange(0,run_length + sample_rate,sample_rate)
   
    # for each run 
    for i in range(n_runs):
        start_idx = i*int(run_length/sample_rate + 1)
        end_idx = start_idx + int(run_length/sample_rate + 1)
        y = mc_e[start_idx:end_idx]
        ax1.plot(x,y, alpha=0.7)
    #plt.yscale('symlog')
    plt.savefig("mc_trajectories.png", format='png', bbox_inches='tight')


def filter_atp_models(minima_indices, n_flow, s_flow):
    filtered_models = []
    for i in minima_indices:  # list of indices of models below threshold
        if (s_flow[i] < 0 and n_flow[i] > 0):
            if ( (abs(s_flow[i] / n_flow[i]) >= 2.9)):
                if (abs(s_flow[i]) > 1e-12 and abs(n_flow[i]) > 1e-12):
                    filtered_models.append(i)
    return filtered_models

def scale_flows(a):
    '''
    Scale a numpy array of flows
    '''
    return np.true_divide(a,np.max(a))


def compare_model_clusters(a, b, raw_data, flow_data, run_length, offset, sample_rate):
    print(a)
    print(b)
    
    def find_flows(index, flow_data, run_length, offset, sample_rate):
        x = np.modf(index)
        mc_n = int(x[1])
        run = int(x[0]*offset+1)
        i = run*int(run_length/sample_rate +1) + int(mc_n/sample_rate)
        return flow_data[i]

    def get_labels(datafile):
        data = np.genfromtxt(datafile, delimiter=',', dtype=None, names = True, invalid_raise = False)
        labels = list(data.dtype.names)
        labels = labels[2:-3]
        return labels

    def graph_model(model, a_dict,s_list, s2_list):
        # input dict: key = 'transition edge name' and value = scaled (by max) net flow
        # input lists: 1) strings to remove for visual clarity (i.e. H1o)
        #              2) strings to insert hyphen for visual clarity (i.e. H1b)
     

        # create a directed graph object
        G = nx.DiGraph()
        # divide into states A and B based on name
        # make state labels easier to read
        # build up graph from each state pair and transition edge
        for key in a_dict:
            #print(key,a_dict[key])  # output for testing
            # 'A_TO_B' --> A and B
            s = key 
            w = ("%1.2f" % a_dict[key])
            a,b = s.split('_TO_',1)
            unbinding = False
            #print(a,b)
            # simplify string (remove i/o/f states)
            for text in s_list:
                a = a.replace(text,"")
                b = b.replace(text,"")
            for text in s2_list:
                if text in a and text not in b:
                    unbinding = True
                a = a.replace(text,"-%s" % text)
                b = b.replace(text,"-%s" % text)
            #print(a,b)

            # check for existing nodes/edges
            # if they don't exist, add them
            if a not in G:
                G.add_node(a)
            if b not in G:
                G.add_node(b)
            if not G.has_edge(a,b):
                if unbinding: # if b in a (i.e. state 'a' w/o a sub-state), then the transition is unbinding
                    G.add_edge(a,b,weight=w, color='r')
                else:
                    G.add_edge(a,b,weight=w, color='b')
        # graph viz. formatting
        plt.close()
        fig = plt.figure(figsize=(25, 10))
        plt.title("Scaled reaction network (%s) of a 3H+ ATPase\nNet flows < 0.25 are omitted" %(model))
        pos = nx.spring_layout(G,k=10/np.sqrt(G.order()))  # more distance
        colors = nx.get_edge_attributes(G,'color').values()
        edges = G.edges()
        widths = [str(int((float(G[u][v]['weight'])*3)+1)) for u,v in edges]
        nx.draw_networkx(G,pos, with_labels=True, edge_color =colors, width=widths)
        labels = nx.get_edge_attributes(G,'weight')
        nx.draw_networkx_edge_labels(G,pos, edge_labels=labels)
        plt.show()
   


    #### Pipeline
    thresh = 0.250
    a_flows = find_flows(a, flow_data, run_length, offset, sample_rate)
    scaled_a_flows = scale_flows(a_flows)
    labels = get_labels(raw_data)
    (a_flows_filt, a_index) = (scaled_a_flows[scaled_a_flows>thresh],np.nonzero(scaled_a_flows > thresh)[0])
 
    a_dict = {}
    for i in range(np.size(a_index)):
        a_dict[labels[a_index[i]]]= a_flows_filt[i]

    b_flows = find_flows(b, flow_data, run_length, offset, sample_rate)
    scaled_b_flows = scale_flows(b_flows)
    labels = get_labels(raw_data)
    (b_flows_filt, b_index) = (scaled_b_flows[scaled_b_flows>thresh],np.nonzero(scaled_b_flows > thresh)[0])

    b_dict = {}
    for i in range(np.size(b_index)):
        b_dict[labels[b_index[i]]]= b_flows_filt[i]
   
    # similar 
    a_keys = list(a_dict.keys())
    b_keys = list(b_dict.keys())
    similar_list = [string for string in a_keys if string in b_keys]
    # different - B
    diff_list_b = [string for string in b_keys if string not in a_keys]
    # different - A
    diff_list_a = [string for string in a_keys if string not in b_keys]
    
    ### output
    print("\nComparing significant scaled net flows (> 0.25) for model %s and %s"  %(a,b))
    print("Shared transitions:")
    print("Transition, net flow A, net flow B")
    for i in similar_list: 
        print("%s, %.2g, %.2g" % (i, a_dict[i], b_dict[i]))
        #print("%-30s, %4.2g, %4.2g" %(i, a_dict[i], b_dict[i]))
    print("\nUnshared transitions:")
    print("Transition, net flow A")
    for i in diff_list_a: 
        print("%s, %.2g" % (i, a_dict[i]))
    print("\nUnshared transitions:")
    print("Transition, net flow B")
    for i in diff_list_b: 
        print("%s, %.2g" % (i, b_dict[i]))

    s_list = ['H1o','H2o','H3o','ATPf']
    s2_list = ['H1b','H2b','H3b','ATPb', 'ADPPIb']
    graph_model(a,a_dict,s_list, s2_list)
    graph_model(b,b_dict,s_list, s2_list)

 # for each experiment
    # get good models (need label x.ab where a is the experiment number and b is the run number)
    # append to good models list
    # cluster these models


 # for each experiment
    # get good models (need label x.ab where a is the experiment number and b is the run number)
    # append to good models list
    # cluster these models


def cluster_atp():
   
      # config
    datafile = "cluster_data_agg.dat"
    run_length=5e3
    n_runs = 10
    sample_rate=1e1
    offset = 1*10**(int(np.log10(n_runs)+1))
    agg = True

    ### preprocessing
    ## input data from file
    data = import_data(datafile)  # cluster_data.dat is [mc_n, mc_e,...flows...]
    ## divide data
    mc_n, mc_e, flow_data, atp_flow, h_flow, _ = slice_array(data[:], run_length=run_length,offset=offset,sample_rate=sample_rate, agg=agg)
    #print(s_flow*abs(s_flow)/(abs(n_flow)+1e-15))  #$emc = $sflow*abs($sflow)/(abs($wflow)+1e-15); #eps for 'fudge factor'
    ## filter data
    filtered_models = filter_atp_models(list(range(0, np.shape(mc_e)[0])), atp_flow, h_flow)  # all models
    ## get list of flows for the filtered models
    minima_flows = get_flows(flow_data, filtered_models)
    minima_models_list = indices_to_models(filtered_models, mc_n)
    ## normalize flows
    try: 
        processed_flows = scale_flows(minima_flows)
    except:
        print ("Error: could not models")
        exit(1)
    ## cluster data
    try:
        cluster_models = create_cluster(processed_flows, minima_models_list)
        clustering = fcluster(cluster_models[1], t=0.65, criterion='distance')
        compare_model_clusters(cluster_models[0][0],cluster_models[0][-1], datafile, flow_data, run_length=run_length,offset=offset,sample_rate=sample_rate)
    except:
        print ("Error: could not cluster models")
        exit(1)
    plot_mc_c(mc_e, n_runs, run_length, sample_rate)
    plot_coupling(h_flow, atp_flow, n_runs, run_length, sample_rate)



def main():
    cluster_atp()
    # print("analyze later...")
    # exit()
    # analyze_experiments()
    # exit()
    # (minima_models_list, processed_flows) = cluster_and_analyze_one_run("cluster_data_1.dat", analyze='None',
    #                           graph=0, proof=1, threshold=1e+100, run_length=5e3, n_runs = 1, sample_rate=1e1)  # "Run1" to analyze

    
    #compare_and_analyze_two_runs("cluster_data.dat","cluster_data_ts1_d1.dat", .01)

    # runs = [("cluster_data_d0_0_1.dat", 0.01),
    #         ("cluster_data_d0_1.dat", 0.1),
    #         ("cluster_data_d02.dat", 0.2),
    #         ("cluster_data_d0_4v2.dat", 0.4),
    #         ("cluster_data_d0_5.dat", 0.5),
    #         ("cluster_data_d0_6v2.dat", 0.6),
    #         ("cluster_data_d08.dat", 0.8),
    #         ("cluster_data_d1.dat",1),
    #         ("cluster_data_d1_25.dat",1.25)]
            
    # #runs = [("cluster_data_d0_1.dat", 0.1)]
    # thresh = [ 0.1, 0.25, 0.5, 0.75, 1.0]
    # #thresh = [0.5]
    # compare_n_runs(runs, thresh)

    #aggregate_analysis()
   

if __name__ == "__main__":
    main()

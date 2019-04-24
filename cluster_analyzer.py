###  This script clusters models from a ModelExplorer simulation, and runs the
###  analysis perl script.
###  August George - 1/7/2018

import numpy as np
import sys
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
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
    data = np.genfromtxt(datafile, delimiter=',', dtype=float, invalid_raise = False)
    if proof ==1:
        assert np.shape(data)[1] == 45, "Data column size mismatch. Should be 45 (flows + n, emc, s,w,n flows) for proofreading."
    return data


### Slices array into three seperate arrays (1st col, 2nd col, remaining col).
### Returns col1 as a row, col2 as a row, remaning columns.
### This is useful for processing the imported cluster_data.dat file [mc_n, mc_e,...flows...]
def slice_array(array):
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
    fig = plt.figure(figsize=(20, 10))
    dn = dendrogram(
        Z,
        truncate_mode = 'level',
        p = 7,
        color_threshold = 0.25,
        distance_sort = 'ascending',
        labels=minima,
    )
    #dn = dendrogram(
        #Z,
        #distance_sort = 'ascending',
        #labels=minima,
    #)
    plt.title('Model Hierarchy [truncated]', fontsize=34)
    plt.ylabel('Cluster Distance', fontsize=30)
    plt.xlabel('Model Number [MC step]', fontsize=30)
    plt.figtext(.80,.85,"Found %s Models" % len(minima), fontsize=30, ha='center')
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.margins(x=0.01, y=0.01)
    plt.savefig("model_heirarchy_cluster_complete.png", format='png', bbox_inches='tight')
    leaves = dn['ivl']  # get leaves of cluster
    cleaned_leaves = [x for x in leaves if str(x).isdigit()]  # remove truncated '(x)'
    return(cleaned_leaves)


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
    plt.title('Monte Carlo Energy Trajectory', fontsize=34)
    plt.ylabel('MC Energy', fontsize=30)
    plt.xlabel('MC Step', fontsize=30)
    ax.plot(x,y,  linewidth=2, color='orange')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    plt.tick_params(axis='both', which='major', labelsize=20)
    #for xy in zip(x2, y2):
    #    ax.annotate('(%s)' % xy[0], xy=xy, textcoords='data')
    plt.legend(fontsize=30)
    plt.margins(x=0.01, y=0.01)
    plt.ylim([-.75e-4, .75e-4])
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
    a_norm = np.zeros_like(a)
    for i in range(np.shape(a)[0]):  # for each row vector
        x_max = np.max(a[i])
        x_min = np.min(a[i])
        x_ptp = x_max - x_min
        for j in range(a[i].size):  # for each element in row vector
            a_norm[i][j] = 2*(a[i][j]-x_min)/(x_ptp) -1
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
    models_list = [int(x) for x in models_list]  # convert to int for use as index
    #for (num,item) in enumerate(models_list):
        #print(num+1,item)
    return models_list


def calculate_s_and_w_flows(flow_data):
    pass


def find_proofreading(s_flows, w_flows, minima_indices, ddg=1.0, n=1):
    proofreading_models = []
    for i in minima_indices:  # list of indices of models below threshold
        if w_flows[i] != 0:
            if (s_flows[i] / w_flows[i]) > (n*np.exp(ddg)):
                proofreading_models.append(i)
    return proofreading_models


### Package all data processing (before clustering) into one subroutine
def process_data(datafile, graph=0, proof=0, threshold = 0):
    #print("processing data...\n")
    data = import_data(datafile, proof)  # cluster_data.dat is [mc_n, mc_e,...flows...]
    mc_n, mc_e, flow_data, n_flow, s_flow, w_flow = slice_array(data)
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
    n = 10
    if proof == 1:
        minima_indices = find_proofreading(
            s_flow, w_flow, minima_indices, ddg=ddg_sw, n =10 )
        n_proof_models = len(minima_indices)
        proof_ratio = 1.0*n_proof_models/n_models
    minima_models_list = indices_to_models(minima_indices, mc_n)
    minima_flows = get_flows(flow_data, minima_indices)
    old_processed_flows = normalize_and_threshold_flows(minima_flows)
    processed_flows = normalize_flows(minima_flows)

    correlation = acf(processed_flows)
    #print(correlation)

    if graph == 1:  # graph mc energy
        graph_mc_data(mc_n, mc_e, minima_indices)
        plt.figure(figsize=(5,5))
        plt.plot(correlation)
        #half_index = (np.abs(correlation - 0.5)).argmin()
        #plt.plot(half_index, correlation[half_index], marker='o', markersize=3, color="red", label = half_index)
        #plt.legend()
        plt.hlines(y=0.5, xmin=0, xmax=len(correlation)-1, linestyles='dotted')
        plt.title("Autocorrelation of MC Run")
        plt.ylabel("Autocorrelation")
        plt.xlabel("Lag time [500 MC steps]")
        #plt.xlim(0,25)
        plt.savefig("mc_auto_correlation.png", format='png', bbox_inches='tight')
    if proof == 0:
        report = "Total Models: %s\nModels below MC energy threshold (%s) : %s, Ratio: %s\n" % (
            n_models, thresh, n_thresh_models, thresh_ratio)
    else:
        report = "Total Models: %s\nModels below MC energy threshold (%s) : %s, Ratio: %s\nModels above proofreading threshold %s*(e^%s): %s, Ratio: %s" % (
            n_models, thresh, n_thresh_models, thresh_ratio, n, ddg_sw, n_proof_models, proof_ratio)
    return processed_flows, minima_models_list, report


### Compares each consecutive model in dataset1 with models in dataset2
### to find the minimum distance
def compare_two_runs(datafile1, datafile2):
    #print("comparing %s to %s...\n" % (datafile1, datafile2))
    processed_flows1, minima_models_list1, report1 = process_data(datafile1)
    processed_flows2, minima_models_list2, report2 = process_data(datafile2)
    # cdist creates a distance matrix between all elements in two matrices
    # np.min on axis 1 finds the min for for each col (each model in run1)
    mindist = np.min(spatial.distance.cdist(processed_flows1, processed_flows2), axis=1)
    return mindist, minima_models_list1, report1, report2


### Plot 2 histograms: Run1 models distance from Run2 models and vice versa
### Should clean up later (don't repeat yourself!)
def plot_two_histograms(run1_data, run2_data):
    #print("graphing histograms...\n")
    fig = plt.figure(figsize=(30, 10))
    plt.suptitle("Comparing Models of Two Runs (under the same conditions)")

    ax1 = fig.add_subplot('121')  # top left histogram (run1)
    ax1.title.set_text('Run1 vs Run2')
    ax1.set_ylabel('Frequency')
    ax1.set_xlabel('Minimum Distance (run1 model to run2 model)')
    bin_size1 = int(round(len(run1_data)**(0.5)))  # quick estimate of bin size is sqrt(n)
    bin_step1 = (max(run1_data) - min(run1_data))/bin_size1
    if bin_step1 == 0:
        bin_step1 = 1
    ax1.hist(run1_data, bins=bin_size1)
    ax1.set_xticks(np.arange(min(run1_data), max(run1_data)+bin_step1, bin_step1))  # use arrange instead of range to handle floats
    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
    plt.figtext(.40,.85,"Run 1 has %s Models" % len(run1_data))

    ax2 = fig.add_subplot('122')  # top right histogram (run2)
    ax2.title.set_text('Run2 vs Run1')
    ax2.set_ylabel('Frequency')
    ax2.set_xlabel('Minimum Distance (run2 model to run1 model)')
    bin_size2 = int(round(len(run2_data)**(0.5)))  # quick estimate of bin size is sqrt(n)

    bin_step2 = (max(run2_data) - min(run2_data))/bin_size2
    if bin_step2 == 0:
        bin_step2 = 1
    ax2.hist(run2_data, bins=bin_size2)
    ax2.set_xticks(np.arange(min(run2_data), max(run2_data)+bin_step2, bin_step2))  # use arrange instead of range to handle floats
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
    plt.figtext(.825,.85,"Run 2 has %s Models" % len(run2_data))
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


# def autocorrelate(a):
#     n = np.shape(a)[0]  # number of 'rows' = number of time steps
#     auto_corr = np.zeros(n)
#     print(n)
#     for t in range(n-1):
#         #print(t)
#         sum = 0.0
#         num = 0.0
#         dem = 0.0
#         #factor = (1.0)/(n-t)
#         factor = 1.0
#         for i in range(n-1-t):
#             num += (factor) * (np.dot((a[i] - np.mean(a)),(a[i+t] - np.mean(a))))
#             dem += (factor) * (np.dot((a[i] - np.mean(a)),(a[i] - np.mean(a))))
#             #sum += factor * (num/dem)
#         sum = num/(1.0*dem)
#         print(sum)
#         auto_corr[t] = sum
#     return auto_corr

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
    #N = data.size
    mean = np.mean(data)
    acf = np.zeros(N)

    # calculate den = <(x(i)-x.mean).(x(i)-x.mean)>
    den = 0.0
    for t in range(N):  # t = 0 to N-1
        den = den + np.dot((data[t]-mean), (data[t]-mean))
    den = den/N

    # calculate num = <(x(t)-x.mean).(x(t+k)-x.mean)>
    for k in range(N):  # k = 0 to N-1
        num = 0.0
        for t in range(N-k):
            num = num + np.dot((data[t]-mean), (data[t+k]-mean))
        num = num/(N-k)
        acf[k] = num/den
    return acf


### Main program follows simple pipeline architecture:
### 1) import raw data, 2) seperate data into mc_n, mc_e, and flows,
### 3) find local minima from mc_n, 4) get net flows of minima,
### 5) normalize net flows of minima, 6) graph mc_n and mc_e w/ minima
### 7) cluster models based on normalized net flows of energy minima models
### 8) run analysis script for each energy minima model
def cluster_and_analyze_one_run(datafile, analyze=None, graph=1, proof=0, threshold=1e+100):
        print("clustering and analyzing one run: %s\n" % datafile)
        start_time = time.time()  # testing for runtime data
        processed_flows, minima_models_list, report = process_data(
            datafile, graph=graph, proof=proof, threshold=threshold)
        cluster_models = create_cluster(processed_flows, minima_models_list)
        if analyze == 'Run1':
            analyze_models(cluster_models)
        runtime = time.time()-start_time  # testing for runtime data
        print(report)
        print(cluster_models)
        if analyze == 'Run1':
            print("Models analyzed: %s" % (len(cluster_models)))
        print("Runtime: %s (s)" % (runtime))
        print(time.strftime("%Y-%m-%d %H:%M"))


def compare_and_analyze_two_runs(datafile1, datafile2, thresh, analyze=None):
    start_time = time.time()  # testing for runtime data
    print("Comparing two runs: %s and %s" % (datafile1, datafile2))
    print("Min distance threshold = %s" % thresh)
    print("Analyze (None if blank): %s \n" % analyze)
    r1_vs_r2, r1_labels, report1, report2 = compare_two_runs(datafile1, datafile2)
    r2_vs_r1, r2_labels, _, _ = compare_two_runs(datafile2, datafile1)
    plot_two_histograms(r1_vs_r2, r2_vs_r1)
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


def main():
    cluster_and_analyze_one_run("cluster_data.dat", analyze=None,
                                graph=1, proof=0, threshold=1e+100)  # "Run1" to analyze
    #compare_and_analyze_two_runs("cluster_data.dat","cluster_data2.dat", .01)
    pass

if __name__ == "__main__":
    main()

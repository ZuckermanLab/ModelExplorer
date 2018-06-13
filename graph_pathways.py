# Transporter Flux Pathway Grapher - August George - 2018
# Takes in a .CSV of state A, state B, and flux A->B (with first row reserved for labels) and makes a directed graph
# of the state space pathways (using relative flux)

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from numpy import genfromtxt
from decimal import Decimal
import pprint
import datetime
import warnings # there is a discrepnecy between python and numpy when comparing strings
warnings.simplefilter(action='ignore', category=FutureWarning) #to ignore the warnings about the discrepency
# Nodes are hardcoded with name/position.
def graph_pathways(sub_weighted_edge_list, sub_image_file_name, sub_analysis_labels = {}):


    weighted_edge_list = sub_weighted_edge_list
    analysis_labels = sub_analysis_labels
    date = datetime.date.today()
    title = "%s [A.G.-%s]\n%s" % (analysis_labels[0],date,analysis_labels[1])
    image_file_name = sub_image_file_name
    G=nx.DiGraph()


    #manually add and position nodes. Not ideal.
    #checks first node pair for 'OF-No-So-Wo' and/or 'OF-No-So'. Not ideal since that states could appear in a different order
    if 'OF-No-So-Wo' in weighted_edge_list[0]: #for proofreading.
        G.add_node("OF-No-So-Wo",pos=(100,150))
        G.add_node("OF-No-Sb-Wo",pos=(50,75))
        G.add_node("OF-Nb-So-Wo",pos=(100,100))
        G.add_node("OF-No-So-Wb",pos=(150,125))
        G.add_node("OF-Nb-Sb-Wo",pos=(75,25))
        G.add_node("OF-Nb-So-Wb",pos=(125,50))
        G.add_node("IF-No-So-Wo",pos=(300,150))
        G.add_node("IF-No-Sb-Wo",pos=(250,75))
        G.add_node("IF-Nb-So-Wo",pos=(300,100))
        G.add_node("IF-No-So-Wb",pos=(350,125))
        G.add_node("IF-Nb-Sb-Wo",pos=(275,25))
        G.add_node("IF-Nb-So-Wb",pos=(325,50))
    elif 'OF-No-So' in weighted_edge_list[0]: #for non-proofreading.
        G.add_node("OF-No-So",pos=(100,150))
        G.add_node("OF-No-Sb",pos=(50,75))
        G.add_node("OF-Nb-So",pos=(150,125))
        G.add_node("OF-Nb-Sb",pos=(100,25))
        G.add_node("IF-No-So",pos=(300,150))
        G.add_node("IF-No-Sb",pos=(250,75))
        G.add_node("IF-Nb-So",pos=(350,125))
        G.add_node("IF-Nb-Sb",pos=(300,25))

    G.add_weighted_edges_from(weighted_edge_list)

    pos=nx.get_node_attributes(G,'pos')
    e_labels = nx.get_edge_attributes(G,'weight')

    nx.draw_networkx_nodes(G,pos, node_size=600)
    nx.draw_networkx_edges(G,pos,  alpha=0.50)
    nx.draw_networkx_edge_labels(G,pos, label_pos=0.5, edge_labels=e_labels)
    nx.draw_networkx_labels(G,pos)
    plt.axis('off')

    #scales figure so it fits
    l,r = plt.xlim()
    b,t = plt.ylim()
    plt.xlim(l-30,r+15)
    plt.ylim(b-5,t+5)
    plt.suptitle("Transporter (Relative) Flux Pathway",)
    plt.title(title, fontsize=6)
    #plt.show()
    plt.savefig("%s.png" % image_file_name, dpi = 150)
    plt.close()



# Data is imported from a .csv with format: (state A, state B, Flux AB)
def import_data(sub_data_file):
    data_file = sub_data_file
    raw_data = np.genfromtxt(data_file, delimiter=',', dtype=None)
    data_labels = []
    data_no_labels = []
    data_labels.append(raw_data[0][0])
    data_labels.append(raw_data[0][1])
    for i in range(len(raw_data)-1):
        data_no_labels.append(raw_data[i+1])
    pprint.pprint(raw_data)
    return data_no_labels, data_labels

# We want data in (state a, state b, flow ab) where the flow is positive (since the arrows correspond to a->b).
# Since every transition has a reverse, go through each tranisition and filter out negative flux.
def sort_data(sub_data):
    data = sub_data
    flows = []
    sorted_data = []

    # Find max flux and calculate minimum threshold for a significant flux (1%)
    for i in range(len(data)):
        #flows.insert(i, data[i][2])
        flows.append(data[i][2])
    max = np.nanmax(flows)
    threshold = 0.01 * max

    # Find transitions with flux > 0 and > threshold. Store in new array.
    for i in range(len(data)):
        if (data[i][2] > 0 and data[i][2] > threshold):
            sorted_data.append(data[i])
    pprint.pprint(sorted_data)

    # Calulate (and store) relative flux
    for i in range(len(sorted_data)):
        sorted_data[i][2] = sorted_data[i][2]/max
        sorted_data[i][2] = '{:.2f}'.format(sorted_data[i][2])
    pprint.pprint(sorted_data)

    return sorted_data

#since analysis.prl prints out 10 flux_n.csv files
for n in range(0,110,10):

    flux_name = "flux_n%s" %n
    analysis_file = "%s.csv" %flux_name

    (data, data_labels) = import_data(analysis_file)
    edge_list = sort_data(data)
    graph_pathways(edge_list, flux_name, data_labels)

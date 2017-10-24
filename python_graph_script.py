# Proof-maker Graph Cycle Script
# Written by August George 2017
# Warning: This software is underdevelopment. Use at your own risk.
#
# This script uses .csv files (state adjacency matrix and relative energies matrix)
# generated from proof-maker to generate graph(s) of the connected states which
# are then checked for consistency using the relative enegies matrix.
#
# Inputs: adjacency matrix .CSV and relative energies matrix .csv
# Output: Error Flag, Dictionary of subgraph energy landscapes (also a dictionary), list of connected states in subgraph
# Note: .CSV files currently have row/column labels and null character at the end of each rows
# Note: Unused subroutine to draw and save a graph image (.png) as well as generate test data
# In progress: fine tuning graph layout for images

from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import json

# creates adjacency graph from adjacency matrix
def adjacency_matrix_to_graph(ADJACENCY_MATRIX):
	adj = ADJACENCY_MATRIX
	rows, cols = np.where(adj == 1) # finds row & col elements that are connected
	edges = zip(rows.tolist(), cols.tolist()) #creates a list of tuples from row & col lists
	gr = nx.Graph() # creates a networkx graph var
	gr.add_edges_from(edges) # adds edges to graph (i.e. creates the graph)
	return gr

# creates list of connected components (i.e. subgraphs, retaining graph info)
def graph_to_subgraph(ADJACENCY_GRAPH):
	gr = ADJACENCY_GRAPH
	subgraph = list(nx.connected_component_subgraphs(gr, copy=True))
	return subgraph

# checks consistency of graph cycle (given a graph and energies matrix)
def check_and_prune_cycle(SUB_ADJACENCY_GRAPH,ENERGIES_MATRIX):
	gr = SUB_ADJACENCY_GRAPH
	en = ENERGIES_MATRIX
	en_sum = 0
	consistent = "x"
	try: #check if there is a cycle at all
		cycle = list(nx.find_cycle(gr)) #lists cycle by edges (i.e. [(a,b),(b,c),(c,a)])
	except nx.exception.NetworkXNoCycle: #catch error if there is no cycle found
		cycle = False #no cycles => cannot check for consistency, exit
		consistent = "No cycles to check for consistency"
		return gr, consistent
	while cycle: #while there is a cycle in the cycle list
		for j in xrange(len(cycle)): #sum energy around loop
			en_sum = en_sum + en[cycle[j][0]][cycle[j][1]] #sum energy for given edge, note: using energies[x][y] => edge (x->y)
		if (en_sum !=0): #if sum of energy is not zero, cycle is inconsistent
			en_sum = 0
			cycle = False #exit
			consistent = ""
			return gr, consistent
		else:
			gr.remove_edge(cycle[0][0],cycle[0][1]) #removes first edge
			en_sum = 0
			try:
				cycle = list(nx.find_cycle(gr))
			except nx.exception.NetworkXNoCycle:
				cycle = False
				consistent = "cycle is consistent"
				return gr, consistent #cycle is consistent, exit

# Load .csv data files and convert them to somthing useful
def import_and_format_data(ADJACENCY_FILE,ENERGIES_FILE):
	raw_adjacency_data = np.genfromtxt(ADJACENCY_FILE, delimiter=',', dtype=None) #guesses data type (is slower)
	raw_energies_data = np.genfromtxt(ENERGIES_FILE, delimiter=',', dtype=None)
	label_data = raw_adjacency_data[1:, 0:1].tolist() #labels from only first column, each row after first row
	n = len(label_data) #dimension of nxn matrix = number of labels
	init_adjacency_matrix = raw_adjacency_data[1:n+1,1:n+1].astype(np.int) #skip labels, ends before null char, make list
	init_energies_matrix = raw_energies_data[1:n+1,1:n+1].astype(np.float) # matrix[x][y] is value at row x, col y
	init_labels = dict(zip(xrange(n), sum(label_data, [])[:])) #"flattens" list of lists into a single list
	return (init_adjacency_matrix, init_energies_matrix, init_labels)

# Draws and saves (.png) graph image. Can add text node labels. Default is numeric node labels.
def draw_and_save_graph(GRAPH, FILENAME, TEXT_LABELS={}):
	gr = GRAPH
	filename = FILENAME #subroutine needs file name to save image to
	text_labels = TEXT_LABELS #list of labels (from initialization subroutine)
	mapping = {} #mapping is a dictionary of nodes and label names. Used to draw text labels
	for i in gr.nodes(): #to avoid errors, only map the nodes used in this graph (subgraphs use less nodes)
		mapping[i] = text_labels[i]
	plt.clf() #clear plot
	if not TEXT_LABELS:
		nx.draw( #draw graph with numeric node labels
			gr,
			node_size=1,
			with_labels=True,
			pos = nx.spring_layout(gr,k=0.15,iterations=20)) #position nodes using spring layout (can be fine-tuned)
		plt.savefig("%s.png" %filename, dpi=250 )
	else:
		nx.draw( #draw graph with text node labels
			gr,
			node_size=1,
			with_labels=True,
			labels = mapping, #use mapping dictionary {[node1:text1],[node2:text2]...} for labels
			pos = nx.spring_layout(gr,k=0.15,iterations=20))  #position nodes using spring layout (can be fine-tuned)
		plt.savefig("%s.png" %filename, dpi=250 )

#Goes through tree graph and creates an "energy landscape" of a subgraph
#Specifically, this creates a dictionary of states and relative energies to a reference state (first node)
#where the key is the reference state (first node). Also, updates master dictionary  and master list of states (for use in multiple subgraphs)
#note: for a tree graph there is only one path between two nodes
def graph_to_data(TREE_GRAPH, ENERGIES, LABELS, MASTER_DICTIONARY, MASTER_SUBGRAPH_LIST):
	gr = TREE_GRAPH
	en = ENERGIES
	lbls = LABELS
	master_dictionary = MASTER_DICTIONARY # dictionary with keys being first node of a subgraph and values being dictionary of states and relative energies
	master_subgraph_list = MASTER_SUBGRAPH_LIST # list of lists containng tied state names for a subgraph
	data = {} #dictionary containing energy landscape
	labeled_data = {} #replaces node numbers with state names (using labels list)
	label_list = [] #list of state labels
	node_list = list(gr.nodes()) # creates list of nodes in subgraph
	for node in node_list: # goes through each node
		label_list.append(lbls[node]) #creates list of labels in subgraph
		if node == node_list[0]: #first node in list
			data[node_list[0]] = 0.0 #set reference node to 0 energy compared to itself (by definition)
			labeled_data[lbls[node_list[0]]] = 0.0
		else:
			path = nx.shortest_path(gr, source = node_list[0], target = node) #find path between reference node and target node
			relative_energy = 0.0 #initlilize relative energy running total
			for i in xrange(len(path)-1): #for each edge in path
				relative_energy = relative_energy + en[path[i]][path[i+1]] # add energy along path to previous energy. note using en[x][y] => edge(x->y)
			data[path[len(path)-1]] = relative_energy #add to dictionary. key = last node in path, value = relative energies added along path
			labeled_data[lbls[path[len(path)-1]]] = relative_energy
	master_dictionary[lbls[node_list[0]]] = labeled_data # add to master dictionary. key = first node, value = dictionary of nodes and relative energies (energy landscape)
	master_subgraph_list.append(label_list) # adds to master list of "lists of connected states in subgraph"

#takes master dictionary, list, and error flag and prints it to a file to be
#read by perl
def export_to_file(DATA,DATA2,FLAG, FILENAME):
	data = DATA
	data2 = DATA2
	flag = FLAG
	filename = FILENAME

	file = open(filename,"w")
	file.write("ERROR FLAG = %s\n\n" %flag)
	file.write("%s\n\n" %data)
	file.write("%s\n\n" %data2)
	file.close()

	with open('data.json', 'w') as outfile:
		json.dump(data, outfile)
#generate data for testing. under development
#def test_module():

#intialializes variables/flags and controls program flow
#Flow: import/format data -> create adjacency graph -> create subgraphs -> check each subgraph for consistency
#->prune subgraph to tree -> create energy landscape data -> export to file
def main():
	adjacency_file = "adjacent_matrix.csv"
	energies_file = "energies_matrix.csv"
	error = False #error flag to send back to perl
	master_dictionary = {} # store all subgraph data into one dictionary that then gets exported into perl. The subgraph "energy landscape" to be used by perl
	master_subgraph_list = [] # store a list of states of each connected subgraph, to be used by perl
	(adjacency_matrix, energy_matrix, state_labels) = import_and_format_data(adjacency_file, energies_file) #import and format data
	adjacency_graph = adjacency_matrix_to_graph(adjacency_matrix) #create adjacency graph
	subgraphs = graph_to_subgraph(adjacency_graph) #create subgraphs

	for k in xrange(len(subgraphs)): #loop through each subgraph
		print "checking subgraph %s" % k
		(pruned_graph, consistent) = check_and_prune_cycle(subgraphs[k],energy_matrix) #set constistency flag for each subgraph
		if not consistent: #check_and_prune_cycle function will set conistent ("" = false, string type used to allow text for "no cycles")
			error = True #set error flag
			print "ERROR! Cycles are not self-consistent!"
		else:
			print consistent #cycle is consistent or there are no cycles
		graph_to_data(pruned_graph,energy_matrix,state_labels, master_dictionary, master_subgraph_list) #adds to master dictionary and list
	export_to_file(master_dictionary,master_subgraph_list,error,"python_to_perl") #exports master dictionary and list

if __name__ == "__main__": #best practice to use main
    main()

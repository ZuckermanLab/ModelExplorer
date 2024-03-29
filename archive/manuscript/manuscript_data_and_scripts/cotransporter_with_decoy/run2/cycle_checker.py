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
# Possible Improvements: fine tuning graph layout for images
#
# IN PROGRESS: Adding transition/barriers  

#todo: check  subgraphs, pruned graphs, and output data

from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import json
import pprint

# creates adjacency graph from adjacency matrix
def adjacency_matrix_to_graph(ADJACENCY_MATRIX):
	adj = ADJACENCY_MATRIX

	print("adj matrix")
	pprint.pprint(adj)

	rows, cols = np.where(adj == 1) # finds row & col elements that are connected

	print("rows")
	pprint.pprint(rows)
	print("cols")
	pprint.pprint(cols)

	edges = zip(rows.tolist(), cols.tolist()) #creates a list of tuples from row & col lists

	print("edges")
	pprint.pprint(edges)

	gr = nx.Graph() # creates a networkx graph var

	print ("gr")
	pprint.pprint(gr)

	gr.add_edges_from(edges) # adds edges to graph (i.e. creates the graph)

	print("gr add edges")
	pprint.pprint(gr)


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
		cycle = False
		consistent = "No cycles to check for consistency"
		return gr, consistent #no cycles => cannot check for consistency, exit
	while cycle: #while there is a cycle in the cycle list
		for j in xrange(len(cycle)): #sum energy around loop
			en_sum = en_sum + en[cycle[j][0]][cycle[j][1]] #sum energy for given edge, note: using energies[x][y] => edge (x->y)
		if (en_sum !=0): #if sum of energy is not zero, cycle is inconsistent
			en_sum = 0
			cycle = False
			consistent = ""
			return gr, consistent #cycle is inconsistent, exit
		else:
			gr.remove_edge(cycle[0][0],cycle[0][1]) #removes first edge
			en_sum = 0
			try:
				cycle = list(nx.find_cycle(gr)) #lists cycle by edges (i.e. [(a,b),(b,c),(c,a)]) (there is one less edge now)
			except nx.exception.NetworkXNoCycle: #catch error if there is no cycle found
				cycle = False
				consistent = "cycle is consistent"
				return gr, consistent #cycle is consistent, exit

# Load .csv data files and convert them to somthing useful (IN PROGRESS, WILL MAKE MORE GENERAL)
def import_and_format_data(ADJACENCY_FILE,ENERGIES_FILE, TRANSITION_ADJACENCY_FILE, TRANSITION_ENERGY_FILE):


	raw_adjacency_data = np.genfromtxt(ADJACENCY_FILE, delimiter=',', dtype=None) #guesses data type (is slower)
	raw_energies_data = np.genfromtxt(ENERGIES_FILE, delimiter=',', dtype=None)
	raw_transition_adjacency_data = np.genfromtxt(TRANSITION_ADJACENCY_FILE, delimiter=',', dtype=None)
	raw_transition_energies_data = np.genfromtxt(TRANSITION_ENERGY_FILE, delimiter=',', dtype=None)

	print("raw_transition_adjacency_data")
	np.set_printoptions(threshold='nan')
	pprint.pprint(raw_transition_adjacency_data)

	label_data = raw_adjacency_data[1:, 0:1].tolist() #labels from only first column, each row after first row
	n = len(label_data) #dimension of nxn matrix = number of labels

	transition_label_data = raw_transition_adjacency_data[1:, 0:1].tolist()
	transition_n = len(transition_label_data)

	print("transition_label_data")
	np.set_printoptions(threshold='nan')
	pprint.pprint(transition_n)
	pprint.pprint(transition_label_data)

	init_adjacency_matrix = raw_adjacency_data[1:n+1,1:n+1].astype(np.int) #skip labels, ends before null char, make list
	init_energies_matrix = raw_energies_data[1:n+1,1:n+1].astype(np.float) # matrix[x][y] is value at row x, col y
	init_transition_adjacency_matrix = raw_transition_adjacency_data[1:transition_n+1,1:transition_n+1].astype(np.int) #skip labels, ends before null char, make list
	init_transitions_energies_matrix = raw_transition_energies_data[1:transition_n+1,1:transition_n+1].astype(np.float) # matrix[x][y] is value at row x, col y

	print("init_transition_adjacency_matrix")
	np.set_printoptions(threshold='nan')
	pprint.pprint(init_transition_adjacency_matrix)

	init_labels = dict(zip(xrange(n), sum(label_data, [])[:])) #"flattens" list of lists into a single list
	init_transition_labels = dict(zip(xrange(transition_n), sum(transition_label_data, [])[:]))

	print("init transition labels")
	np.set_printoptions(threshold='nan')
	pprint.pprint(init_transition_labels)

	return (init_adjacency_matrix, init_energies_matrix, init_labels, init_transition_adjacency_matrix, init_transitions_energies_matrix, init_transition_labels)

# Draws and saves (.png) graph image. Can add text node labels. Default is numeric node labels.

def draw_and_save_graph(GRAPH, FILENAME, TEXT_LABELS={}):
	gr = GRAPH
	filename = FILENAME #subroutine needs file name to save image to
	text_labels = TEXT_LABELS #list of labels (from initialization subroutine)
	mapping = {} #mapping is a dictionary of nodes and label names. Used to draw text labels
	if TEXT_LABELS:
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
	gr = TREE_GRAPH #must use tree type data structure
	en = ENERGIES
	lbls = LABELS
	master_dictionary = MASTER_DICTIONARY # dictionary with keys being first node of a subgraph and values being dictionary of states and relative energies
	master_subgraph_list = MASTER_SUBGRAPH_LIST # list of lists containng tied state names for a subgraph
	data = {} #dictionary containing energy landscape
	label_list = [] #list of state labels
	labeled_data = {} #replaces node numbers with state names (using labels list)
	node_list = list(gr.nodes()) # creates list of nodes in subgraph

	for node in node_list: # goes through each node
		label_list.append(lbls[node]) #creates list of labels in subgraph
		if node == node_list[0]: #first node in list
			data[node_list[0]] = 0.0 #set reference node to 0 energy compared to itself (by definition)
			labeled_data[lbls[node_list[0]]] = 0.0
		else:
			path = nx.shortest_path(gr, source = node_list[0], target = node) #find path between reference node and target node (list of nodes)
			print("path %s" % path)
			relative_energy = 0.0 #initlilize relative energy running total
			for i in xrange(len(path)-1): #for each edge in path
				relative_energy = relative_energy + en[path[i]][path[i+1]] # add energy along path to previous energy. note using en[x][y] => edge(x->y)
				print("relative_energy %s" % relative_energy)
			data[path[len(path)-1]] = relative_energy #add to dictionary. key = last node in path, value = relative energies added along path
			labeled_data[lbls[path[len(path)-1]]] = relative_energy
	master_dictionary[lbls[node_list[0]]] = labeled_data # add to master dictionary. key = first node, value = dictionary of nodes and relative energies (energy landscape)
	master_subgraph_list.append(label_list) # adds to master list of "lists of connected states in subgraph"

#takes master dictionary, list, and error flag and prints it to a file to be read by perl
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

	with open('%s.json' %filename, 'w') as outfile:
		json.dump(data, outfile)

#generate data for testing. under development
#def test_module():

#intialializes variables/flags and controls program flow
#Flow: import/format data -> create adjacency graph -> create subgraphs -> check each subgraph for consistency
#->prune subgraph to tree -> create energy landscape data -> export to file
def main():
	adjacency_file = "adjacent_matrix.csv"
	energies_file = "energies_matrix.csv"
	transition_adjacency_file = "transition_state_adjacent_matrix.csv"
	transition_energy_file = "transition_state_energy_matrix.csv"

	error = False #error flag to send back to perl
	master_dictionary = {} # store all subgraph data into one dictionary that then gets exported into perl. The subgraph "energy landscape" to be used by perl
	master_subgraph_list = [] # store a list of states of each connected subgraph, to be used by perl
	master_transition_dictionary = {} # store all subgraph data into one dictionary that then gets exported into perl. The subgraph "energy landscape" to be used by perl
	master_transition_subgraph_list = [] # store a list of states of each connected subgraph, to be used by perl

	print("importing and formating data...")
	(adjacency_matrix, energy_matrix, state_labels, transition_adjacency_matrix, transition_energy_matrix, transition_labels) = import_and_format_data(adjacency_file, energies_file, transition_adjacency_file, transition_energy_file) #import and format data

	np.set_printoptions(threshold='nan')
	print("transition_adjacency_matrix main")
	pprint.pprint(transition_adjacency_matrix)

	print("state adjacency matrix to graph...")
	adjacency_graph = adjacency_matrix_to_graph(adjacency_matrix) #create adjacency graph

	print("transition adjacency matrix to graph...")
	transition_adjacency_graph = adjacency_matrix_to_graph(transition_adjacency_matrix)

	print("transition_adjacency_graph main")
	pprint.pprint(transition_adjacency_graph)

	print("transition energy matrix")
	pprint.pprint(transition_energy_matrix)

	draw_and_save_graph(transition_adjacency_graph, "transition_adjacency_graph")
	draw_and_save_graph(adjacency_graph, "adjacency_graph")

	print("graphs to subgraphs...")
	subgraphs = graph_to_subgraph(adjacency_graph) #create subgraphs
	transition_subgraphs = graph_to_subgraph(transition_adjacency_graph)

	print("check and prune subgraphs...")
	print("checking state subgraphs...")
	for k in xrange(len(subgraphs)): #loop through each subgraph
		print "checking subgraph %s" % k
		(pruned_graph, consistent) = check_and_prune_cycle(subgraphs[k],energy_matrix) #set constistency flag for each subgraph
		if not consistent: #check_and_prune_cycle function will set conistent ("" = false, string type used to allow text for "no cycles")
			error = True #set error flag
			print "ERROR! Cycles are not self-consistent!"
		else:
			print consistent #cycle is consistent or there are no cycles
		graph_to_data(pruned_graph,energy_matrix,state_labels, master_dictionary, master_subgraph_list) #adds to master dictionary and list
	export_to_file(master_dictionary,master_subgraph_list,error,"python_to_perl_state") #exports master dictionary and list

	error = False
	print("checking transition state subgraphs...")
	for k in xrange(len(transition_subgraphs)): #loop through each subgraph
		print "checking subgraph %s" % k
		(transition_pruned_graph, consistent) = check_and_prune_cycle(transition_subgraphs[k],transition_energy_matrix) #set constistency flag for each subgraph
		if not consistent: #check_and_prune_cycle function will set conistent ("" = false, string type used to allow text for "no cycles")
			error = True #set error flag
			print "ERROR! Cycles are not self-consistent!"
		else:
			print consistent #cycle is consistent or there are no cycles
		graph_to_data(transition_pruned_graph,transition_energy_matrix,transition_labels, master_transition_dictionary, master_transition_subgraph_list) #adds to master dictionary and list
	export_to_file(master_transition_dictionary,master_transition_subgraph_list,error,"python_to_perl_transition") #exports master dictionary and list


if __name__ == "__main__": #best practice to use main
    main()

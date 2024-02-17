This folder contains data and scripts for four Monte Carlo simulations of a transporter with a decoy substrate, searching for a high substrate/decoy flux ratio.

Each simulation run folder contains:

Scripts:
*ModelExplorer.prl: Main perl script which uses Markov chain Monte Carlo to generate kinetic models (set of rate constants and associated state/barrier energies)
*ModelAnalyzer.prl: Secondary perl script used to analyze existing models (i.e. state/barrier energies) to provide flow data over a range of simulation conditions
*cycle_checker.py: Python script used within both ModelExplorer and ModelAnalyzer to check for cycle consistency (i.e. energy detail balance wthin a chemical cycle)
*flux_grapher.py: Python script used w/ ModelAnalyzer to automate the flux pathway diagrams
*python_solver.py: Python script used w/ ModelExplorer to aid in solving for rate matrix in pathological cases
*graph_mc_traj.py: Python script used to graph MC trajectory in model space using cluster_data

Files:
* config, and analysis_config files set the parameters used by ModelExplorer and ModelAnalyzer respectively. 
* cluster_data contains the MC step, MC energy, and netflows for each model. Used for cluster analysis based on Euclidean distance of the net flows of different models. 
** cluster_data format: MC n, MC energy, ... net flows between states ..., ion (N) flow, substrate (S) flow, decoy (S) flow
* models_from_run contains files for both state and barrier energies for each model. 

The aggregate folder contains:

Scripts:
*cluster_analyzer.py: Python script that imports cluster_data from four runs, generates/plots a hierarchical cluster, and plots the MC and cluster trajectories.

Files:
* cluster_data files of four seperate MC simulations

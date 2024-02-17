This folder contains data and scripts for two Monte Carlo simulations of a cotransporter: a symporter and an antiporter

Each simulation run folder contains:

Scripts:
*ModelExplorer.prl: Main perl script which uses Markov chain Monte Carlo to generate kinetic models (set of rate constants and associated state/barrier energies)
*ModelAnalyzer.prl: Secondary perl script used to analyze existing models (i.e. state/barrier energies) to provide flow data over a range of simulation conditions
*cycle_checker.py: Python script used within both ModelExplorer and ModelAnalyzer to check for cycle consistency (i.e. energy detail balance wthin a chemical cycle)
*graph_mc_traj.py: Python script used to graph MC trajectory in model space using cluster_data

Files:
* config, and analysis_config files set the parameters used by ModelExplorer and ModelAnalyzer respectively. 
* cluster_data contains the MC step, MC energy, and netflows for each model. Used for cluster analysis based on Euclidean distance of the net flows of different models. 
** cluster_data format: MC n, MC energy, ... net flows between states ...
* models_from_run contains files for both state and barrier energies for each model. 
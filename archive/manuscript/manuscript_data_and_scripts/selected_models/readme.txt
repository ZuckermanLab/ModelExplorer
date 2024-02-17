This folder contains the flux data for 7 cotransporter models and an accompanying graphing script.

The flux data for a given model is given for a range a ion chemical potential differences (dmu). 
Each flux data file is formated: dmu_ion value, ion (N) flux, substrate (S) flux, decoy substrate (W) flux 

Scripts:
*graph_relative_flow.py: Python script that graphs the flux (flow) of each model

Data:
*antiporter_flux: flux data for an antiporter (no W flux)
*symporter_flux: flux data for a symporter (no W flux)
*modelA_flux: flux data for representative model of cluster A (run 1, MC n 3500)
*modelB_flux: flux data for representative model of cluster B (run 1, MC n 29000)
*modelB_flux_noleak: flux data for representative model of cluster B (run 1, MC n 29000) with the ion leak removed
*modelC: flux data for representative model of cluster C (run 3, MC n 3000)
*modelD: flux data for representative model of cluster D (run 3, MC n 829000)

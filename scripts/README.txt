Transporter State Space Modeling README (outdated - need to update!)

NOTE: The software is frequently updated to include better pipelining, automation, configuration, and analysis capabilities, as well as any bugs that we find.
The readme might not be 100% up-to-date on the latest changes to the GitHub repo. Automation is still a WIP. Use at your own risk!  

README.txt written by August George and Dr. Zuckerman.
Proof-maker.prl written by Dr. Zuckerman and August George.
Analyze-model.prl written by Dr. Zuckerman and by August George.

Readme TODO: METHODOLOGY,more assumptions?, update with bug fix and improvements, update everything!.

TABLE OF CONTENTS:
1. INTRODUCTION
2. REQUIREMENTS
3. INPUT/OUTPUT
4. ASSUMPTIONS (in progress)
5. USAGE
6. METHODOLGY/ALGORITHM/notes (in progress)
7. THEORY
8. BUGS/To-Do
9. CHANGELOG (major updates)
10. ACKNOWLEDGEMENTS


1. INTRODUCTION:

This program package includes: proof-maker.prl, analyze-model.prl, python_graph_script.py, README.txt, and a LICENSE.
For legal/copyright concerns, please read LICENSE.

Proof-maker.prl uses monte carlo techniques to find a range of models which satisfy the given state space system constraints.
Analyze-model.prl uses "successful" models from proof-maker.prl to calculate the flows with varying chemical potential.


This software was developed to investigate the behavior of cell membrane transporter proteins. Specifically, the programs
creates an abstract state space representation of a transporter and then finds potential cycle pathways which optimize sugar
flow in and toxin flow out. These models can then be analyzed to give the flows of sugar, sodium, and toxin flow for
a varying chemical potential. The overall goal is the use this data to better classify transporter cycle regimes. These regimes
can then be used to aid experimental work and the theoretical understanding of transport cycles.


2. REQUIREMENTS:

All programs use Perl 5 and were tested/run on *nix machines.
*Proof-maker.prl uses Perl PDL to solve system of linear equations
*Analyze-model.prl uses Python/Numpy to "solve" system of linear equations (left over from proof-maker base code)


3. INPUT/OUTPUT:

Now using config file! Will update readme soon...

Proof-maker:
*Inputs: None. Parameters (i.e. number of monte carlo steps) adjusted by user in program file.
*Outputs:
**Console Output: Outputs starting parameters and state space configurations during monte carlo loop.
**rate_matrix.tmp: temporary file which contains running log of rates (k) in matrix form
**prob_ss.dat: temporary file which stores steady state values in a matrix.
**evolver_flows.dat: temporary file which stores running log of configuration states for trial pertubation
**evovler_rates.dat: fiel which stores monte carlo step number and associated energy function calculation
**models_from_run: folder which contains energies-x and barriers-x
***energies-x: stores state energy values at model (monte carlo) number 'x'.
***barriers-x: stores barrier energy values at model (monte carlo) number 'x'.

Proof-maker:
*Inputs: None.
*Outputs:
**Console Output: outputs initial parameters and state space configurations during monte carlo loop.
**Rate_matrix.tmp: temporary file which contains running log of rates (k) in matrix form
**Prob_ss.dat: temporary file which stores steady state values in a matrix.
**Evolver_flows.dat: temporary file which stores running log of configuration states for trial pertubation
**Evovler_rates.dat: data file which stores monte carlo step number and associated energy function calculation
**Mc_out.dat: data log file which stores state space and energy information for each monte carlo run
**Models_from_run: folder which contains energies-x and barriers-x
***Energies-x: stores state energy values at model (monte carlo) number 'x'.
***Barriers-x: stores barrier energy values at model (monte carlo) number 'x'.

Analyze-model:
Note: analyze-maker uses the same code base as proof-maker so it is has similar outputs, not all of
which are useful...

*Inputs:
**energies-x: stores state energy values at model (monte carlo) number 'x'. Outputed by proof-maker.
**barriers-x: stores barrier energy values at model (monte carlo) number 'x'. Outputed by proof-maker.

*Outputs:
**Analysis-vary_dmu_N-dmu_init__-6__to__dmu_fin__0: file that stores Sodium, Sugar, and Toxin flows
**at a given range of dMu (ussualy dMu_Na).
**Console Output: outputs initial parameters and state space configurations during monte carlo loop.
**Rate_matrix.tmp: temporary file which contains running log of rates (k) in matrix form
**Prob_ss.dat: temporary file which stores steady state values in a matrix.
**Evolver_flows.dat: temporary file which stores running log of configuration states for trial pertubation
**Evovler_rates.dat: data file which stores monte carlo step number and associated energy function calculation
**Mc_out.dat: data log file which stores state space and energy information for each monte carlo run
**Models_from_run: folder which contains energies-x and barriers-x
***Energies-x: file that stores state energy values at model (monte carlo) number 'x'.
***Barriers-x: file stores barrier energy values at model (monte carlo) number 'x'.


DEFAULT PARAMATERS:

### Simulation type
$proof = 1; # 1 if proofreading (additional states and transitions present compared to simple transport)
$na_first = 1; # 1 if 'sodium' (N) is forced to bind first

### MONTE CARLO PARAMETERS ###
$nsteps = 1e4; # number of Monte Carlo steps
$dprint = 1e0; # interval between prints
$n_beta = 2e2; # number of MC steps after which beta (inverse temperature) changes
$seed = 456789; # seed for random number generator

$demax = 1.0; # maximum change in energy of state or transition state
$alpha = 2.0; # exponent for toxin flow in MC 'energy' function which sets fitness of model
$beta_init = 1e1; # initial value of beta = inverse temperature
$fbeta = 1e3; # max factor by which beta can change during the MC run, which includes 'tempering'
$tol = 0.3; # fractional energy change deemed significant in tempering
$pbeta_stay = 0.2; # probability to stay at same beta during tempering
#$beta_early = 1e1;  # alternate between two values
#$beta_late = 1e3;
### PHYSICAL PARAMETERS (everything in kT units, effectively)
$dmu_N = -4; # chemical potential change (mu_i - mu_o) of driving "ion" N (e.g., sodium)
$dmu_S = 2; # chemical potential change of "substrate" S (e.g., sugar)
$dmu_W = 2; # chemical potential change of wrong "substrate" W (e.g., toxin)
$fmu_N = 0.5; # fraction of dmu attributed to outside - i.e., log of concnetration increase relative to equil
$fmu_S = 0.5; # fraction of dmu attributed to outside - i.e., log of concnetration increase relative to equil
$fmu_W = 0.5; # fraction of dmu attributed to outside - i.e., log of concnetration increase relative to equil
$ebarslow = 99; # minimum energy value of barrier height (above max of two states) for slow transitions
$kzero = 1e-3; # rate constant prefactor: full first-order rate constant = k_ij = kzero * exp( - barrier_energy )
$ebump = 1.0;  # amount by which initial transition-state (barrier) energy exceeds max of pair


4. ASSUMPTIONS/QUESTIONS(in progress):

# Sb states are tied to Wb states to maintain the dg_sw difference between them.

*WE ARE TIEING BARRIERS FOR Sb and Wb (OF-Nb->IF-Nb = OF-Wb->IF-Wb ) AND (No->Nb = Wo->Wb)
By drawing the free energy levels, you will see that this is equivalent to the usual assumption (e.g., Hopfield's in the original kinetic proofreading) 
that the substrates differ only in their off-rates. Also, tying the barriers together makes the substrates smoothly equivalent as dg_SW -> 0.  

*Should OF-IF conformational transistions be equivalent (when N, S, and W are held constant)?
*Same barrier for W and S (un)binding. Should barrier difference = 0 for s & w?
*Na must bind first (based on experimental insight)
*Steady State Flow: constant source/drain of chemical potenial (i.e. sodium)
*Should toxin be strictly 2KT higher energy than sugar? (dg_SW=2?)
*Is W vs S binding on/off rate being changed? Does is it depend on relative values?



5. USAGE:

Note: Both programs were tested and used on *nix machines.
Note: Negative Monte Carlo energies correspond to a "fit" model according to our energy function
Note: Positive flow is defined as outside->inside

Programs can be run from terminal using "$ perl -w program_name.prl" (or similar)
It is advised to save console output using "$ perl -w program_name.prl > datalog" (or similar)
The srand(x) function is used for repeatability. A simulation for the same x value
should give the same results.
Proof-reading and Na binding can be turned on/off by setting their associated variables to 0.

Proof-maker:
In general, set the alpha, nsteps, and random seed value and then run the program.

Some useful variables/parameters:
*$nsteps: sets how many steps in a simulation
*$emc: defines mc energy function. Default is -$sflow *abs($sflow/$wflow)^$alpha.
*$alpha: meta parameter which controls the emc defined above.
*$seed: random number seed used in srand() function
*$proof: turns proofreading on/off (default is 1 = ON)
*$na_first: requires Na to bind first (default is 1 = ON)
*$dMu_N: chemical potential change of driving force (Mu_i - Mu_o) Na (default is -4).
*$dMu_S: chemical potential change of sugar (default is 2)
*$dMu_W: chemical potential change of toxin (default is 2)
*$dg_SW: binding energy difference between S and W (default is 2) <--Asumption!

*$analysis: variable in analyze-model which determines which dMu to vary (Na or toxin)
*$efile_init: variable in analyze-model which stores state energy filename
*$bfile_init: variable in analyze-model which stores barrier energy filename
*$dmu_init: variable in analyze-model which stores initital dMu value (default -6 for Na)
*$dmu_fin: variable in analyze-model which stores final dMu value (default 0 for Na)

Analyze-model:
In general, find "successful" state and barrier energy files (located in models_from_run folder).
Put these file names into script along with ?? alpha, nsteps, and random seed value from
proof-maker run ?? then run the program.

Simple Analysis:
Plot evolver_rates.dat (gnuplot is one simple method) and look for negative valleys. Pick the
local min (approx). These negative values correspond "successful" models given the monte carlo
energy function. Note the mc n number for these local minima and find the associated energies-x
and barriers-x files.
Use the files to run the analyze-model script which will output
Analysis-vary_dmu_N-dmu_init__-6__to__dmu_fin__0. This file contains flows of N, S, and W.
Plot these flows vs dMu. To investigate stoichiometry, plot the ratio of the flows (i.e. N/W vs dMu)


6. METHODOLOGY (in progress):

High Level algorithm:
Subroutines:

Note: $nchange is used to see if two states are equivalent (for binding of the other other species)
### (i.e. for OF: NoSiWo->NoSbWo->NoSoWo differs by one "equivalent" unit Ni/No, No/Ni from NiSiWo->NiSbWo->NiSoWo)

7. THEORY:

Transporters:
*These proteins consist of: uniporters, symporters, and antiporters.
*Uniporters transport one molecule across the cell membrane.
*Symporters transport two (or more) molecules across the cell membrane, in the same direction.
*Antiporters transport two (or more) molecules across the cell mambrane, in opposite direction.
*These transporter proteins are powered by the chemical potential differences (outside/inside cell) of each molecule.
*By utilizing a large chemical potential difference of one molecule (ex Na), transporters are able to move another
*molecule (ex Sugar) against its smaller (relative) chemical potential difference.

Cycles and State Space Representation:
*In general, a simple complete symporter cycle consists of molecules binding to the protein, the protein changing
*conformation, the molecules unbinding, and the protein changing conformation again.
*In this way, a particular cycle can be described by a molecule's binding state (bound/unbound) and the transport
*protein's conformation state (facing in/out).
*Finally, a matrix of every possible configuration can be constructed using these two states (conformation and binding)
*creating the transporter "state space".

Complex cycle dynamics:
*Transporters in nature do not neccesarily follow the simple cycle described above.
*Slippage might occur, in which the cycle will not produce a 1:1 ratio of molecules
*(i.e. a sugar is transported for each sodium).
*This slippage may be considered a leak, and is less energy efficient.
*These cycle ineffiencies may have a functional use, providing an extra cycle step to filter out (proofread)
*unwanted toxins.

For a more in-depth discussion on cell membrane transporter proteins please see:
http://www.physicallensonthecell.org/ - Online textbook written by Dr. Zuckerman

For a more in-depth discussion on Transporter Cycles and Energies:
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill

For a more in-depth discussion on statistcal mechanics methods applied to cell biology:
"Statistical Physics of Biomolecules: An Introduction" by DM Zuckerman


8. BUGS/To-Do:

W vs S binding on or off rate being changed? Does it depend on relative values?

Does Monte Carlo change the slowies?

Fix equivalent transitions (proof-reading section) subroutine to be cleaner (some how replace "n_bases - 2" )

Correct Initialized Model?

Better was to fix randomness issue?

test that all tied members are consistent with any/all constraints ( encode and run through all constraints ...) 

small bug: cannot use strict
small bug: transition states are printed into evolver_rates.dat

*BUG FIXED* Incorrect equivalent energy state groupings -> python script written to  configure/validate energy state groupings
*BUG FIXED* Incorrect equivalent barrier energy groupings -> proof-maker algorithm added to configure groupings. python script written to validate energy state groupings
*BUG FIXED* Runs were not repeatable -> sorted @barkeys so each run has the same ordered list for trial moves.

Proof-maker Unit Tests to Add:
(1) Sum of steady state probabilities is one
(2) Ratio of rates for a state pair give Boltz fac of energy difference
(3) Energy differences among tied states match constraints (state "energy landscape" is set in python script, then initialized in perl)
(4) Energies of fixed states remain at constrained values (state "energy landscape" is set in python script, then initialized in perl)
(5) Reject/Restore step properly restores energies (done in perl)
(6) Tied states graphs are self-consistent (done in python, but can be better integrated into perl code, which currently loads boolean "consistent" variable from file) 

Proof-maker Improvements:
*Should OF-IF conformational transistions be equivalent (when N, S, and W are held constant)?
*"use strict" errors out (currently commented out)
*Distance subroutine is a bottleneck (only needs to be called once, before MC loop starts)
*Optimize PDL matrix solving funtion
*Rewrite much of the code into subroutines that can be aclled in main
*Optimize code to run on cluster
*Port (or rewrite) programs to python/numpy
*Consolidate all programs into one program with streamlined input/output files
*Expand abstraction to include more complex systems

Analyze model:
*still calls python/numpy to solve matrix (bottle neck)

General Tweaks:
*Clean code/comments
*Fine tune tempering
*Other energy functions?
*other ways to generate more models per run?

9. CHANGE LOG:

Update 5: (2017-04-01): Added improved pipeline, automation, analysis, and graphing capabilities. Added functionality to tie together states for dg_SW ( might still have bugs). Added unit tests 
for energy restoration, cycle consitency, energy consistency and list consistency. 

Update 4 (2017-12-6): "Transition  energies not restored after rejection " bug has been fixed using new python script. Similar issue as "State energies not restored after rejection" bug.
Implemented new proof-maker algorithm to correctly group equivalent transitions and reused python script to validate.
"Randomness" bug has been fixed. Runs with same config did not give identical output. Issue was that perturbed barriers were being picked from an unsorted list @barkeys so each run was different.
Sorted @barkeys before it is used.

Update 3 (2017-11-1): "State energies not restored after rejection" bug has been fixed using new python script.
Issue was the equivalent states were not all completely connected.
Python creates connected state graphs based on adjacency matrix, checks cycle for consistency using relative energies,
and creates new connected states dictionary (of dictionary) containing state and connected states with relative energies (i.e. energy landscape).
Config file functionality is added.

Update 2 (2017-08-13): Perl PDL used in place of Python/Numpy to elimate system IO call slow down in proof-maker. Uploaded as proof-maker2.prl (will keep original for historical purposes).
Uploaded analyze-model.prl. Added more to README.

Update 1 (2017-07-13): README.txt created by August George. Added table of
contents and changelog. Minor changes to comments in proof-maker.prl.
Github repo setup.

Update 0 (2016-17): Proof-maker.prl program created by Dr. Dan Zuckerman


10. ACKNOWLEDGEMENTS:

Thanks to Dr. Zuckerman for his work writing the proof-maker and analyze-model programs, and also for his
helpful advice and mentoring.

Oregon Health and Sciences University (OHSU)

http://www.physicallensonthecell.org/ - Useful online resource written by Dr. Zuckerman
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill

stackexchange

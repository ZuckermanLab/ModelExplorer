Transporter State Space Modeling README

README.txt written by August George with help from Dr. Zuckerman.
Proof-maker.prl written by Dr. Zuckerman and updated/modified by August George.
Analyze-model.prl written by Dr. Zuckerman and updated/modified by August George. 

Last updated: September 13, 2017

!BUG! - Energy state not restored to last accepted state after rejection. Need to fix and then add unit test. 

TODO: METHODOLOGY,more assumptions?, fix "restore state after rejection" bug, add unit test to check energies have been restored 

TABLE OF CONTENTS:
1. INTRODUCTION 
2. REQUIREMENTS
3. INPUT/OUTPUT
4. ASSUMPTIONS (in progress)
5. USAGE
6. METHODOLGY/ALGORITHM (in progress)
7. THEORY
8. BUGS
9. CHANGELOG
10. ACKNOWLEDGEMENTS


1. INTRODUCTION:

This program package includes: proof-maker.prl, analyze-model.prl, README.txt, and a LICENSE. 
For legal/copyright concerns, please read LICENSE.

Proof-maker.prl uses monte carlo techniques to find a range of models which satisfy the given state space system constraints. 
Analyze-model.prl uses "successful" models from proof-maker.prl to calculate the flows with varying chemical potential. 

This software was developed to investigate the behavior of cell membrane transporter proteins. Specficially, the programs
creates an abastract state space representation of a transporter and then finds potential cycle pathways which optimize sugar 
flow in and toxin flow out. These models can then be analyzed to give the flows of sugar, sodium, and toxin flow for
a varying chemical potential. The overall goal is the use this data to better classify transporter cycle regimes. These regimes
can then be used to aid experimental work and the theoretical understanding of transport cycles. 


2. REQUIREMENTS:

All programs use Perl 5 and were tested/run on *nix machines.
*Proof-maker.prl uses Perl PDL to solve system of linear equations
*Analyze-model.prl uses Python/Numpy to "solve" system of linear equations (left over from proof-maker base code)


3. INPUT/OUTPUT:

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


4. ASSUMPTIONS(in progress):

*Na must bind first (based on experimental insight)
*Steady State Flow: constant source/drain of chemical potenial (i.e. sodium)
*Should toxin be strictly 2KT higher energy than sugar? 
*Should barrier difference = 0 for s & w?

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
*$dMu_N: chemical potential change of driving force (Mu_i - Mu_o) Na (default is 6). 
*$dMu_S: chemical potential change of sugar (default is 2)
*$dMu_W: chemical potential change of toxin (default is 2)
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


8. BUGS/Improvements:

!BUG! - Energy state not restored to last accepted state after rejection. Need to fix and then add unit test. 

Proof-maker:
*"use strict" errors out (currently commented out)
*Distance subroutine is a bottleneck

analyze model:
*still calls python/numpy to solve matrix (bottle neck)

Possible Improvements:
*Add unit tests ( sum of probabilites = 1, energy before trial move = energy after rejection/restore step)
*Optimize PDL matrix solving funtion
*Optimize code to run on cluster
*Port (or rewrite) programs to python/numpy
*Consolidate all programs into one program with streamlined input/output files
*Expand abstraction to include more complex systems


9. CHANGE LOG:

Update 2 (2017-08-13): Perl PDL used in place of Python/Numpy to elimate system IO call slow down in proof-maker. Uploaded as proof-maker2.prl (will keep original for historical purposes). Uploaded analyze-model.prl. Added more to README. 

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

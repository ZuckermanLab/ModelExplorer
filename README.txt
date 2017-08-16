Transporter State Space Modeling README

README.txt written by August George with help from Dr. Zuckerman.
Proof-maker.prl written by Dr. Zuckerman and updated/modified by August George.
Analyze-model.prl written by Dr. Zuckerman and updated/modified by August George. 
??Data-sort.prl (in progress) written by August George??

Last updated: August 14, 2017

TODO: rearrange sections and redo ToC, continue writing usage/assumptions, continue 
writing data in/out and files, METHODOLOGY, exapnd on introduction,  hdd requirements?, min specs?, graph script?

Table of Contents:
1. Introduction
2. Theory
3. Requirements/Configuration
4. Usage/Assumptions
5. Files
6. Data In/Out
7. Analysis
8. Bugs/Improvements
9. CHANGE LOG
10. Acknowledgements

INTRODUCTION:

This package includes: proof-maker, analyze-model, and data-sort(WIP). 

This software was developed to investigate the behavior of cell membrane transporter proteins. Specficially, the programs
creates an abastract state space representation of a transporter and then finds potential cycle pathways which optimize sugar 
flow in and toxin flow out. These pathways can be sorted and analyzed to give the flows of sugar, sodium, and toxin flow for
a different chemical potentials. The overall goal is the use this data to better classify transporter cycle regimes, which 
experimental scientists can then use. 

Proof-maker.prl uses monte carlo techniques to find a range of models which satisfy the given state space system constraints. 
Analyze-model.prl uses "successful" models from proof-maker.prl to calculate the flows with varying chemical potential. 
Data-sort.prl (WIP) bridges proof-maker.prl and analyze-model.prl by finding "successful" models given by proof-maker.prl output. 

REQUIREMENTS:

All programs use Perl 5. 
*Proof-maker.prl uses Perl PDL to solve system of linear equations
*Analyze-model uses Python/Numpy to "solve" system of linear equations (left over from proof-maker template)

USAGE:

Proof-maker usage:


METHODOLOGY: 

Proof-maker high level algorithm:
*

analyze-model high level algorithm:
*

?? Data-sort high level algorithm:
*Opens proof-maker output data (evolver_rates.dat) which contains MC number and MC Energies. 
*Searches for all negative values (negative = successful model) and writes them to a temp file.
*Searches for valleys in negative data values within a user determined variance
*Writes the bottom data values of each valley to a data file
??

ASSUMPTIONS:

State Space Assumptions:
*Na must bind first (based on experimental insight)
*Steady State Flow: constant source/drain of chemical potenial (i.e. sodium)

THEORY:

For a more in-depth discussion on cell membrane transporter proteins please see:
http://www.physicallensonthecell.org/ - Online textbook written by Dr. Zuckerman

For a more in-depth discussion on Transporter Cycles and Energies:
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill 

For a more in-depth discussion on statistcal mechanics methods applied to cell biology: 
"Statistical Physics of Biomolecules: An Introduction" by DM Zuckerman

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


BUGS/IMPROVEMENTS:

Proof-maker:
*"use strict" errors out (currently commented out)
*Distance subroutine is a bottleneck

analyze model:
*still calls python to solve matrix (bottle neck)

??data-sort:
*need algorithm/function to find valleys within certan variance
??

Possible Improvements:
*Optimize PDL matrix solving funtion
*Port (or rewrite) programs to python/numpy
*Consolidate all programs into one program with streamlined input/output files
*Expand abstraction to include more complex systems

CHANGE LOG:

Update 2 (2017-08-13): Perl PDL used in place of Python/Numpy to elimate system IO call slow down in proof-maker. Uploaded as proof-maker2.prl (will keep original for historical purposes). Uploaded analyze-model.prl. Added more to README. 

Update 1 (2017-07-13): README.txt created by August George. Added table of
contents and changelog. Minor changes to comments in proof-maker.prl.
Github repo setup.

Update 0 (2016-17): Proof-maker.prl program created by Dr. Dan Zuckerman

ACKNOWLEDGEMENTS:

Thanks to Dr. Zuckerman for his work writing the proof-maker and analyze-model programs, and also for his 
helpful advice and mentoring. 

Oregon Health and Sciences University (OHSU)

http://www.physicallensonthecell.org/ - Useful online resource written by Dr. Zuckerman
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill 


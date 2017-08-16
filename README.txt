Transporter State Space Modeling README

README.txt written by August George with help from Dr. Zuckerman.
Proof-maker.prl written by Dr. Zuckerman and updated/modified by August George.
Analyze-model.prl written by Dr. Zuckerman and updated/modified by August George. 

Last updated: August 14, 2017

Table of Contents:
1. Introduction
2. Background/Theory
3. Requirements/Configuration
4. Operation
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
flow in and toxin flow out. 

Proof-maker.prl uses monte carlo techniques to find a range of models which satisfy the given state space system constraints. 
Analyze-model.prl uses "successful" models from proof-maker.prl to calculate the flows with varying chemical potential. 
Data-sort.prl (WIP) bridges proof-maker.prl and analyze-model.prl by finding "successful" models given by proof-maker.prl output. 

THEORY:

For a more in depth discussion on cell membrane transporter proteins please see:
http://www.physicallensonthecell.org/ - Online resource written by Dr. Zuckerman

For a more in depth discussion on Transporter Cycles and Energies:
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill 

Transporters:
*These proteins consist of: uniporters, symporters, and antiporters. 
*Uniporters transport one molecule across the cell membrane. 
*Symporters transport two (or more) molecules across the cell membrane, in the same direction. 
*Antiporters transport two (or more) molecules across the cell mambrane, in opposite direction. 
*These transporter proteins are powered by the chemical potential differences (outside/inside cell) of each molecule. 
*By utilizing a large chemical potential difference of one molecule (ex Na), transporters are able to move another 
*molecule (ex Sugar) against its smaller (relative) chemical potential difference. 

Tranporter Cycles and State Space Representation:


CHANGE LOG:

Update 2 (2017-08-13): Perl PDL used in place of Python/Numpy to elimate system IO call slow down in proof-maker. Uploaded as proof-maker2.prl (will keep original for historical purposes). Uploaded analyze-model.prl. Added more to README. 

Update 1 (2017-07-13): README.txt created by August George. Added table of
contents and changelog. Minor changes to comments in proof-maker.prl.
Github repo setup.

Update 0 (2016-17): Proof-maker.prl program created by Dr. Dan Zuckerman

ACKNOWLEDGEMENTS:

README.txt written by August George with help from Dr. Zuckerman.
Proof-maker.prl written by Dr. Zuckerman and updated/modified by August George.
Analyze-model.prl written by Dr. Zuckerman and updated/modified by August George. 
Data-Sort.prl (in progress) written by August George. 

Oregon Health and Sciences University (OHSU)

http://www.physicallensonthecell.org/ - Useful online resource written by Dr. Zuckerman
"Free Energy Transduction and Biochemical Cycle Kinetics" by Terrel L. Hill 


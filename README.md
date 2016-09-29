# FastComp

This repository contains the FastComp algorithm for subcellular localisation prediction of metabolic reactions.

##Algorithm properties
* LP approximation of a MILP to determine reaction activity, the MILP is formulated in [1]
* Adapted FASTCC algorithm to allow simultaneous positive activation of all reactions.

##Repository structure
src contains the source files with CalculateCompartments the script that has to be run to get the prediction.
In additon the src folder contains the models used in the paper (iND750_Model.mat, Recon2ForFastComp.mat, ReconForFastComp.mat), 
along with scripts to calculate samples (CalcCompartStatisticsX.m, with X being either iND, Recon or Recon2).
These scripts take three arguments, the known percentage, the number of replicates and a count that can be used to split the computation into multiple runs.
All three scripts assume, that there is a FastCompResults directory in the users home folder and they will store their results in that folder.
The results can be analysed using the DisplayAllResults function from the Evaluation folder.


### Relevant Literature
[1] Mintz-Oron, S.; Aharoni, A.; Ruppin, E. & Shlomi, T. Network-based prediction of metabolic enzymes' subcellular localization Bioinformatics, 2009, 25, i247-i252
[2] Vlassis, N.; Pacheco, M. P. & Sauter, T. Fast reconstruction of compact context-specific metabolic network models. PLoS Comput Biol, 2014, 10, e1003424

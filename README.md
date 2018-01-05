# FastComp
This repository contains the FastComp algorithm for subcellular localisation prediction of metabolic reactions.

## Algorithm properties
* LP approximation of a MILP to determine reaction activity, the MILP is formulated in [1]

## Requisites
Two versions of fastComp exist:
A Version directly accessing cplex  and using the original fastcore algoirthm from [2] (used in the paper) and a version using the solve commands provided by the COBRA Toolbox[3].
If you use the former, make sure, that you have cplex installed and set up on the matlab path.

## Repository structure
* src contains the source files with CalculateCompartments the script that has to be run to assign reactions and create a model wih localisations.
* src/Comparison contains all scripts necessary to run the different models along with the corresponding data. 
* Decompartmentalisation contains the scripts used for decompartmentalisation.
* Evaluation contains the scritps used to create the Results figures (without the actual data).
All scripts in src are documented in the Standard COBRA Toolbox format. 

### Relevant Literature
[1] Mintz-Oron, S.; Aharoni, A.; Ruppin, E. & Shlomi, T. **Network-based prediction of metabolic enzymes' subcellular localization Bioinformatics**, 2009, 25, i247-i252  
[2] Vlassis, N.; Pacheco, M. P. & Sauter, T. **Fast reconstruction of compact context-specific metabolic network models.** PLoS Comput Biol, 2014, 10, e1003424
[3] Laurent Heirendt & Sylvain Arreckx, Thomas Pfau, Sebastian N. Mendoza, Anne Richelle, Almut Heinken, Hulda S. Haraldsdottir, et al. **Creation and analysis of biochemical constraint-based models: the COBRA Toolbox v3.0** (submitted), 2017, `arXiv:1710.04038 <https://arxiv.org/abs/1710.04038>`__.

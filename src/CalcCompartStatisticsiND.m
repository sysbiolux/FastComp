function CalcCompartStatisticsiND(Percentage,Replicates,filename,useMO,randomize)
%Calculates compartments for the given number of Replicates and the given
%percentage of known reaction localisations for the iND yeast model.
%It should always yield the same unknown reactions if randomiz is not set
%to 0. The function further assumes, that CPLEX is on the MATLAB path.
rng('shuffle')
if nargin < 4
    useMO = 0;
end
if nargin < 5
    seed = randi(intmax);
    disp('Using Random unknown sets')
else
    if(randomize)
        seed = randi(intmax);
        disp('Using Random unknown sets')
    else
        seed = 0;
    end
end
cpath = pwd;
addpath([pwd filesep 'CobraFunctions']);
addpath([pwd filesep])
addpath([pwd filesep 'FastCore' filesep])
addpath([pwd filesep 'Logic' filesep])
mkdir(['/tmp/FC_' num2str(Percentage) '_' num2str(seed)])
cd(['/tmp/FC_' num2str(Percentage) '_' num2str(seed)])
load('iND750_Model.mat')

[CompartResults,ResultFC,ResultMO,Predictions]= CalculateSampleForKnownPercentage(iND750_decomp_DirAdjust,1,'c',iND750_DirecAdjust_CompIDs,iND750_DirecAdjust_CompNames,0,iND750_DirecAdjust_OrigLocs,iND750_DirecAdjust_NonExt,Replicates,Percentage,seed,useMO)

cd(cpath);
save(filename,'CompartResults','ResultFC','ResultMO','Predictions');


function CalcCompartStatisticsiND(Percentage,Replicates,rngSeed,useMO,filename)
%Calculates compartments for the given number of Replicates and the given
%percentage of known reaction localisations for the iND yeast model.
%It should always yield the same unknown reactions if randomiz is not set
%to 0. The function further assumes, that CPLEX is on the MATLAB path.
rng('shuffle')
if nargin < 4
    useMO = 1; %By default this also runs the Mintz Oron algorithm.
end
if ~exist('rngSeed','var')
    %Do non random sampling
    seed = 0;
else
    seed = rngSeed;    
end
if ~exist('filename','var')
    filename = [getenv('HOME') filesep 'FastCompResults' filesep 'ResultiND750' num2str(Percentage) '-' num2str(seed)];
end

addpath([pwd filesep])
%Use the direct fastcore implementation.
addpath([pwd filesep 'FastCore' filesep])
mkdir([tempdir filesep  'FC_' num2str(Percentage) '_' num2str(seed)])
cd([tempdir filesep 'FC_' num2str(Percentage) '_' num2str(seed)])
load('iND750_Model.mat')

[CompartResults,ResultFC,ResultMO,Predictions]= CalculateSampleForKnownPercentage(iND750_decomp_DirAdjust,1,'c',iND750_DirecAdjust_CompIDs,0,iND750_DirecAdjust_OrigLocs,iND750_DirecAdjust_NonExt,Replicates,Percentage,seed,useMO)

save(filename,'CompartResults','ResultFC','ResultMO','Predictions');


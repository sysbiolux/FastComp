function CalcCompartStatisticsRecon(Percentage,Replicates,rngSeed,useMO,filename)
%Calculates compartments for the given number of Replicates and the given
%percentage of known reaction localisations for the human reconstruction Recon 1.
%It should always yield the same unknown reactions if randomiz is not set
%to 0. The function further assumes, that CPLEX is on the MATLAB path.
rng('shuffle')
if nargin < 4
    useMO = 0; %By default this does not run the Mintz oron algorithm
end
if ~exist('rngSeed','var')
    %Do non random sampling
    seed = 0;
else
    seed = rngSeed;    
end
if ~exist('filename','var')
    filename = [getenv('HOME') filesep 'FastCompResults' filesep 'ResultRecon' num2str(Percentage) '-' num2str(seed)];
end
addpath([pwd filesep])
addpath([pwd filesep 'FastCore' filesep])
mkdir([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
cd([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
load('Recon1ForFastComp.mat')

[CompartResults,ResultFC,ResultMO,Predictions]= CalculateSampleForKnownPercentage(Recon_Decomp,1      , 'c'      ,Recon_CompIDs , 0              , Recon_OrigLocs  , Recon_NonExtReacs, Replicates, Percentage, seed,useMO)

save(filename ,'CompartResults','ResultFC','ResultMO','Predictions');


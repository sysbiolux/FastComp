function CalcCompartStatisticsRecon(Percentage,Replicates,rngSeed,useMO,filename)
% Calculates a number of FastComp predictions for the Recon 1 human model.
%
% USAGE:
%
%    CalcCompartStatisticsiND(Percentage,Replicates,rngSeed,useMO,filename)
%
% INPUTS:
%
%    Percentage:        The percentage of unknown localisations.
%    Replicates:        The number of replicates
%    rngSeed:           Whether to use a specific RNG Seed. If non is
%                       if none is choosen 0 is used, and the results
%                       should be deterministic.
%    useMO:             Whether to use the Mintz-Oron algorithm for
%                       comparison. If this is set to 0, the results stored
%                       in the MO algorithm will be those returned from the
%                       single compartment based prediction.
%    filename:          The FileName to store the results in.
%
% .. Author:
%       - Thomas Pfau
%
% NOTE: To evaulate the results please use the evaluation scripts.
%
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
    if ~exist([getenv('HOME') filesep 'FastCompResults'],'dir')
        mkdir([getenv('HOME') filesep 'FastCompResults'])
    end
    filename = [getenv('HOME') filesep 'FastCompResults' filesep 'ResultRecon' num2str(Percentage) '-' num2str(seed)];
end
addpath([pwd filesep])
addpath([pwd filesep 'FastCore' filesep])
mkdir([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
cd([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
folder = fileparts(which(mfilename));

load([ folder filesep 'Data/Recon1ForFastComp.mat'])

[CompartResults,ResultFC,ResultMO,Predictions]= CalculateSampleForKnownPercentage(Recon_Decomp,1      , 'c'      ,Recon_CompIDs , 0              , Recon_OrigLocs  , Recon_NonExtReacs, Replicates, Percentage, seed,useMO)

save(filename ,'CompartResults','ResultFC','ResultMO','Predictions');


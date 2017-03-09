function CalcCompartStatisticsRecon2(Percentage,Replicates,filename,useMO,randomize)
%Calculates compartments for the given number of Replicates and the given
%percentage of known reaction localisations for the human reconstruction Recon 2.
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
mkdir([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
cd([tempdir filesep 'FC_Recon' num2str(Percentage) '_' num2str(seed)])
load('Recon2ForFastComp.mat')

[CompartResults,ResultFC,ResultMO,Predictions]= CalculateSampleForKnownPercentage(decompHumanWOExchangers,1,'c',[Recon_CompIDs],{},0,comps,1:numel(decompHumanWOExchangers.rxns),Replicates,Percentage,seed,useMO,Exchangers);

cd(cpath);
save(filename ,'CompartResults','ResultMO','ResultFC','Predictions');


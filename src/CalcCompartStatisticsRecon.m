function CalcCompartStatisticsRecon(Percentage,Replicates,Count, filename)
%try
    upath = userpath;
    addpath([pwd filesep 'CobraFunctions'])
    addpath([pwd filesep])
    addpath([pwd filesep 'FastCore' filesep])
    addpath([pwd filesep 'Logic' filesep])
    addpath([getenv('HOME') filesep 'CPLEX' filesep 'cplex' filesep 'matlab'])
    mkdir(['/tmp/FC_Recon' num2str(Percentage) '_' num2str(Count)])
    %cd(['/tmp/FC_Recon' num2str(Percentage) '_' num2str(Count)])
    load('ReconForFastComp.mat')
    
    [PrepTimes,ResultFCPure,ResultFC,ResultMO,ResultFCPlus, Predictions]= CalculateSampleForKnownPercentage(Recon_Decomp,1      , 'c'      ,Recon_CompIDs  , Recon_CompIDs   , 0              , Recon_OrigLocs  , Recon_NonExtReacs, Replicates, Percentage, Count)
                                                                                                           
    home = getenv('HOME')
    save([home filesep 'FastCompResults' filesep 'ResultsRecon_new' num2str(Percentage) '-' num2str(1+((Count-1)*Replicates))],'ResultFCPlus','PrepTimes','ResultFCPure','ResultFC','ResultMO','Predictions');


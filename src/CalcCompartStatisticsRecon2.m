function CalcCompartStatisticsRecon2(Percentage,Replicates,Count, filename)
%try
    upath = userpath;
    addpath([pwd filesep 'CobraFunctions']);
    addpath([pwd filesep])
    addpath([pwd filesep 'FastCore' filesep])
    addpath([pwd filesep 'Logic' filesep])

    mkdir(['/tmp/FC_Recon' num2str(Percentage) '_' num2str(Count)])
    %cd(['/tmp/FC_Recon' num2str(Percentage) '_' num2str(Count)])
    load('Recon2ForFastComp.mat')
    
    [ResultFCPure,ResultFC,Predictions]= CalculateSampleForKnownPercentage2(decompHumanWOExchangers,1,'c',[Recon_CompIDs],{},0,comps,1:numel(decompHumanWOExchangers.rxns),Replicates,Percentage,Count,Exchangers);
    home = getenv('HOME');
    save([home filesep 'FastCompResults' filesep 'ResultsRecon2_new' num2str(Percentage) '-'     num2str(1+((Count-1)*Replicates))] ,'ResultFCPure','ResultFC','Predictions');


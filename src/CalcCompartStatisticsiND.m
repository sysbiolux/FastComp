function CalcCompartStatisticsiND(Percentage,Replicates,Count, filename)
%try
    upath = userpath;
    addpath([pwd filesep 'CobraFunctions']);
    addpath([pwd filesep])
    addpath([pwd filesep 'FastCore' filesep])
    addpath([pwd filesep 'Logic' filesep])
    addpath([getenv('HOME') filesep 'CPLEX' filesep 'cplex' filesep 'matlab'])
    
    mkdir(['/tmp/FC_' num2str(Percentage) '_' num2str(Count)])
    cd(['/tmp/FC_' num2str(Percentage) '_' num2str(Count)])
    load('iND750_Model.mat')  
        
    [PrepTimes,ResultFCPure,ResultFC,ResultMO,ResultFCPlus, Predictions]= CalculateSampleForKnownPercentage(iND750_decomp_DirAdjust,1,'c',iND750_DirecAdjust_CompIDs,iND750_DirecAdjust_CompNames,0,iND750_DirecAdjust_OrigLocs,iND750_DirecAdjust_NonExt,Replicates,Percentage,Count)
    home = getenv('HOME')
    save([home filesep 'FastCompResults' filesep 'ResultiND750_new' num2str(Percentage) '-' num2str(1+((Count-1)*Replicates))],'ResultFCPlus','PrepTimes','ResultFCPure','ResultFC','ResultMO','Predictions');
    

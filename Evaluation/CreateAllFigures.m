%CreateAllFigures
algos = {'Consistency','FastComp','MILP','RandomResult'};
clearvars -except algos
load Recon2ForFastComp.mat
Recon2Fig = DisplayAllResults([40 60 80],50,'Results/ResultRecon2',['c', Recon_CompIDs],decompHumanWOExchangers,algos([1 2 4]),0);
print(Recon2Fig,'Recon2.svg','-dsvg')
clearvars -except algos
load Recon1ForFastComp.mat
ReconFigure = DisplayAllResults([40 60 80],50,'Results/ResultRecon',['c', Recon_CompIDs],Recon_Decomp,algos([1 2 4]));
print(ReconFigure,'Recon1.svg','-dsvg')

clearvars -except algos
load iND750_Model.mat
iNDFig = DisplayAllResults([40 60 80],50,'Results/ResultiND750',['c', iND750_DirecAdjust_CompIDs],iND750_decomp_DirAdjust,algos);
print(iNDFig,'yeast.svg','-dsvg')
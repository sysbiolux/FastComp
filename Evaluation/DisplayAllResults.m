function f = DisplayAllResults(Percentages,SampleCount,FileName,Comps,model,Algos,excludeC,onlypred)
%INPUT: Percentages: 	(an array of percentages generated) e.g. [40 60 80]
%	SampleCount: 	Identifiers for final samples e.g. 20
%	FileName:	The base fielname to look up. If 'Results' is provided, the evaluator will look up all files that match to 'ResultsX-Y.mat' where X is an element of Percentages and Y is an element of 1:samplecount
%	Comps:		The list of Compartments that were predicted (has to be the same as provided to the Compartmentalisation algorithm, but including the cytosol)
%	model:		the model used for prediction (the unlocalised one)
%	Algos:		The Algorithms that should eb displayed. A selection from: {'Compartment','FastComp','Mintz-Oron','RandomResult'}, only the choosen results will be displayed
%OPTIONAL Input:
%	excludeC:	an indicator whether to exclude the cytosol from the calculations (default false);

if nargin < 7
    excludeC = 0;
end

if nargin < 8
    onlypred = 0;
end
algonames = {'Consistency','FastComp','MILP','RandomResult'};

usedAlgos = ismember(algonames,Algos);
f = figure;

Colors = jet(4);
       
       
Colors =Colors(usedAlgos,:);
algocount = sum(usedAlgos);
Times = cell(0);
maxtime = -Inf;
for perc=1:numel(Percentages)
    
    Percentage = Percentages(perc);
    [FCResult, FCPureResult,MOResult,RandomResult] = DisplayResults(Percentage,SampleCount,FileName,Comps,model,excludeC,onlypred);    
    
    bars = [FCResult(1), FCPureResult(1), MOResult(1), RandomResult(1)];
    errors = [FCResult(2), FCPureResult(2), MOResult(2), RandomResult(2)];
    
    bars = plotBars(f,numel(Percentages),4,(perc-1)*4 + 1, errors(usedAlgos),bars(usedAlgos));    
    for i=1:algocount
        bars(i).FaceColor = Colors(i,:);
    end
    set(gca,'YLim',[0 1])
    ylabel([num2str(Percentages(perc)) '%'])
    if perc== numel(Percentages)
       xlabel('EMR')
    end
    bars = [FCResult(3), FCPureResult(3), MOResult(3), RandomResult(3)];
    errors = [FCResult(4), FCPureResult(4), MOResult(4), RandomResult(4)];
   bars = plotBars(f,numel(Percentages),4,(perc-1)*4 + 2, errors(usedAlgos),bars(usedAlgos));   
    for i=1:algocount
        bars(i).FaceColor = Colors(i,:);
    end
    set(gca,'YLim',[0 1])
    if perc== numel(Percentages)
       xlabel('L-Score')
    end
    bars = [FCResult(5), FCPureResult(5),  MOResult(5), RandomResult(5)];
    errors = [FCResult(6), FCPureResult(6),  MOResult(6), RandomResult(6)];
    bars = plotBars(f,numel(Percentages),4,(perc-1)*4 + 3, errors(usedAlgos),bars(usedAlgos));
    for i=1:algocount
        bars(i).FaceColor = Colors(i,:);
    end
    set(gca,'YLim',[0 1])
    if perc== numel(Percentages)
       xlabel('R-Score')
    end
    ctimes = [FCResult(7), FCPureResult(7), MOResult(7), RandomResult(7)] + [FCResult(8), FCPureResult(8), MOResult(8), RandomResult(8)];
    maxtime = max(union(maxtime, ctimes(usedAlgos)));
    bars = [FCResult(7), FCPureResult(7), MOResult(7), RandomResult(7)];
    errors = [FCResult(8), FCPureResult(8), MOResult(8), RandomResult(8)];
    bars = plotBars(f,numel(Percentages),4,(perc-1)*4 + 4, errors(usedAlgos),bars(usedAlgos)); 
    Times{perc} = gca;
    for i=1:algocount
        bars(i).FaceColor = Colors(i,:);
    end
    if perc== numel(Percentages)
       xlabel('Runtime')
    end
    
end

decimals = log10(maxtime);
divisor = 10^(fix(decimals)-1);
maxtime = (fix(maxtime / divisor) + 1)*divisor; %Round to the next highest
for t=1:numel(Times)
    Times{t}.YLim = [0 maxtime];
    try       
        Times{t}.YAxis.Exponent = floor(log10(maxtime)-0.1);
    catch
        %Happens if Matlab Version is prior to 2015a, in this case just
        %ignore it.
    end
end
leg = legend(algonames(usedAlgos));
f.PaperSize = [8 5.7];
f.PaperPosition = [-0.6, -0.5, 9.2 ,6.5];
leg.Position = [0.76, 0.425,0.16 0.15];

end
    
function pbars = plotBars(cfig,subplotx,subploty,subplotid,errors,bars)
set(0,'CurrentFigure', cfig);
currentpos = floor((subplotid-1)/subploty) * (subploty+1) + mod(subplotid-1,subploty) + 1;
subplot(subplotx,subploty+1,currentpos);
if size(errors,1) == 1
    errors(2,:,:) = 0;
    bars(2,:) = 0;    
    pbars = barwitherr(errors,bars);
    set(gca,'XLim',[0.5 1.5]);
    set(gca,'XTickLabels', [],'XTick',[]);
else
    pbars = barwitherr(errors,bars);
    
end
end

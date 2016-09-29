function [PrepTimes,ResultFCPure,ResultFC,ResultMO,ResultFCPlus,Predictions]= CalculateSampleForKnownPercentage(model, epsilon, cytosolID, CompartmentIDs, CompartmentNames, SingleCompModel, OrigLocalisation, NonExternals, replicates, percentage, rngseed, Exchangers )
%The results are double arrays with the following entries:
%Correct Predictions | False Predictions | No Prediction | RunTime
%OrigLocalisation is a list of localisations for all reactions, that are
%not external and not transporters from external into the cytosol
%NonExternals is the list corresponding to the position of those reactions
%in the model.

%initialize the random number generator with the given seed (so that we can
%produce comparable results in different runs
rngseed
rng(rngseed)
%if no exchangers were defined, we create an empty exchanger set.
if nargin < 12
    Exchangers = {};
end
%initialize the return struct
Predictions = struct();
%Initialize the empty compartmentalisation data array to be filled later
ReactionCompartmentalisationData = cell(numel(model.rxns),1);
ReactionCompartmentalisationData(:) = {{}};
%The same for the gene compartmentalisation array
GeneComp = {};
%Initialize the algorithm specific Results vectors
PrepTimes = [];
ResultFCPure = [];
ResultFCPlus = [];
ResultFC = [];
ResultMO = [];
%Remove exchange reactions from the Original localisation, as those will be
%automatically added to the exchange reactions and could otherwise cause
%issues.
[imps,exps] = findImpsAndExpRxns(model);
exchange = find(imps|exps);
toRemove = ismember(NonExternals,exchange);
OrigLocalisation(toRemove) = [];
NonExternals(toRemove) = [];

%Now perform the prediction for the given number of replicates
for i=1:replicates
    %Determine the known and unknown set of reactions
    ReactionCompartmentalisationData(:) = {{}};    
    %Set the rng to the current sample of the given rngseed.
    knowns = datasample(1:numel(OrigLocalisation),floor(numel(OrigLocalisation)*percentage/100),'Replace',false);  %nonexchange,floor(numel(nonexchange)*percentage/100),'Replace',false);
    unknowns = setdiff(1:numel(OrigLocalisation),knowns); 
    ReactionCompartmentalisationData(NonExternals(knowns)) = OrigLocalisation(knowns);
    ToPredict = OrigLocalisation(unknowns);
    %save('Target','ToPredict');
    %Start the prediction for this set
    [FCComps,FCtime,FullComps,Fulltime,FullPlus,FullPlusTime,MOComps,MOTime] =  determineCompartments_compare( model, epsilon, {},...
                                                       cytosolID, CompartmentIDs, CompartmentNames,...
                                                       ReactionCompartmentalisationData, GeneComp,...
                                                       SingleCompModel, 1,...
                                                       '[e]',Exchangers);
    Predictions.(['Replicate' num2str(i)]) = {FCComps,FullComps,FullPlus,MOComps,OrigLocalisation(unknowns),model.rxns(NonExternals(unknowns))};
    PrepTimes(i,1) = 0;    
    %Determine some statistics.
    [corr, wrong, missing, unpred] =  checkprediction(FCComps,OrigLocalisation(unknowns));
    ResultFCPure(i,:) = [corr, wrong, missing, unpred, FCtime];
    [corr, wrong, missing, unpred] =  checkprediction(FullComps,OrigLocalisation(unknowns));    
    ResultFC(i,:) = [corr, wrong, missing, unpred, Fulltime];
    [corr, wrong, missing, unpred] =  checkprediction(FullPlus,OrigLocalisation(unknowns));        
    ResultFCPlus(i,:) = [corr, wrong, missing, unpred, FullPlusTime];
    [corr, wrong, missing, unpred] =  checkprediction(MOComps,OrigLocalisation(unknowns));     
    ResultMO(i,:) = [corr, wrong, missing, unpred, MOTime];    
    fprintf('New RNG Seed is %i\n',(rngseed + i));    
    rng(rngseed+i)
end

end

function [correct, wrong, missing, unpred] = checkprediction(prediction,original)
%input are two cell arrays of cell arrays. we need to check whether each
%cell array has the items of the other cell array contained if yes, it is 
%a correct prediction, otherwise it is a wrong prediction. if prediction
%of one entry is empty, there is no prediction
correct = 0;
wrong = 0;
unpred = 0;
missing = 0;
for entry = 1:size(prediction,1)
    predicted = prediction{entry};
    orig = original{entry};
    if numel(predicted) == 0
        unpred = unpred + 1;
        missing = missing + numel(orig);
    else
        correct = correct + numel(intersect(predicted,orig));
        wrong = wrong + numel(setdiff(predicted,orig));
        missing = missing + numel(setdiff(orig,predicted));
    end
end
end 


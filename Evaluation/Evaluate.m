function [EMR, L_Score, R_Score,origin,pred] = Evaluate(origdata,predictiondata,Comps,excludec,randomsample,onlypredicted)
%Get the evaluation based on Retrieval Score Labelling Score and Exact Match Ratio
%Godbole and Sunita - Discriminative Methods for Multi-Labeled Classification
%
%input are 2 cell arrays of cell arrays

if nargin < 4
    excludec = 0;
end

if nargin < 5
    randomsample = 0;
end

if nargin < 6
    onlypredicted = false;
end

if onlypredicted
    nonpred = cellfun(@isempty, predictiondata);
    predictiondata = predictiondata(~nonpred);
    origdata = origdata(~nonpred);
end

[origin, pred] =generateROCData(origdata,predictiondata,excludec,randomsample,Comps);
EMR = ExactMatchRatio(origin,pred);
L_Score = LabellingFScore(origin,pred);
R_Score = retrievalScore(origin,pred);
end



function [origbin,predbin] = generateROCData(origdata,predictiondata,excludec,randomsample,Comps)
%First, obtain the different possibilities
classes = Comps;

if excludec
    classes = setdiff(classes,'c');
end
origbin = zeros(size(origdata,1),numel(classes));
predbin = zeros(size(origdata,1),numel(classes));
for entry = 1:numel(origdata)    
    origbin(entry,:) = ismember(classes,origdata{entry});    
end

if randomsample
    predbin = createRandomSample(origbin');
else
    for entry = 1:numel(origdata)    
        predbin(entry,:) = ismember(classes,predictiondata{entry});    
    end
    predbin = predbin';    
end


origbin = origbin';


end


function randdata = createRandomSample(origdata)

randdata = zeros(size(origdata));
for i=1:size(origdata,2)
    pos = numel(find(origdata(:,i)));
    sample = randsample(1:size(origdata,1),pos);
    randdata(sample,i) = 1;
end
end

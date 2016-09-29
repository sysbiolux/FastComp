function fscore = LabellingFScore(origdata,preddata)
    
matches = (origdata == 1 & origdata == preddata);
%Eliminate cases with no label in origdata

%totalclassdata = origdata | preddata;
itemscores = sum(matches,1);
totalscores = sum(origdata,1) + sum(preddata,1);
%Remove empty stuff
itemscores(totalscores == 0) = [];
totalscores(totalscores==0) = [];
fscore = sum(2* itemscores./totalscores) / size(origdata,2);


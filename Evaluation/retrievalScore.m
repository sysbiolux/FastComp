function score = retrievalScore(origdata,preddata)

matches = (origdata == 1 & origdata == preddata);
classscores = sum(matches,2);
totalscores = sum(origdata,2) + sum(preddata,2);
classscores(totalscores == 0) = [];
totalscores(totalscores == 0) = [];
score = sum(2 * classscores./totalscores) / size(origdata,1);



function ratio = ExactMatchRatio(origdata,preddata)


exactmatches = all(origdata == preddata,1);

ratio = sum(exactmatches / (size(origdata,2)));
end
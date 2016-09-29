function dist = subsetAcc(origdata,targetdata)

dist = 0;
for i=1:size(origdata,2)
    dist = dist + all(origdata(:,i) == targetdata(:,i));
end
dist = dist / size(origdata,2);
       
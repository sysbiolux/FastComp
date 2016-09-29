function acc = overallAcc(origdata,targetdata)

acc = 0;
for i = 1:size(origdata,1)
    acc = acc + (sum(origdata(i,:) == targetdata(i,:)) /size(origdata,2));
end
acc = acc/size(origdata,1);
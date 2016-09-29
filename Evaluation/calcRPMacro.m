function [Precision,Recall,f1] = calcRPMacro(origdata,preddata)

sizes = sum(origdata,2);
TP = origdata == 1 & origdata == preddata;
FP = origdata == 1 & origdata ~= preddata;
TN = origdata == 0 & origdata == preddata;
FN = origdata == 0 & origdata ~= preddata;

Precision = 0;
Recall = 0;

for i = 1:size(sizes,1)
    if(sum(TP(i,:)) > 0 | (sum(FP(i,:) > 0 & sum(FN(i,:) > 0))))
        Recall = Recall + ( sum(TP(i,:)) / (sum(TP(i,:)) + sum(FN(i,:))));
        Precision = Precision + ( sum(TP(i,:)) / (sum(TP(i,:)) + sum(FP(i,:))));
    end
end
Recall = Recall / size(sizes,1);
Precision = Precision/ size(sizes,1);
f1 = (2*Recall*Precision) / (Recall + Precision);
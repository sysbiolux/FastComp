function [Precision,Recall,f1] = calcRPMicro(origdata,preddata)

sizes = sum(origdata,2);
TP = origdata == 1 & origdata == preddata;
FP = origdata == 1 & origdata ~= preddata;
TN = origdata == 0 & origdata == preddata;
FN = origdata == 0 & origdata ~= preddata;

Precision = 0;
Recall = 0;

Recall = sum(TP) / (sum(TP) + sum(FN));
Precision = sum(TP) / (sum(TP) + sum(FP));
%Recall = Jaccard_similarity(origdata,preddata);
f1 = (2*Recall*Precision) / (Recall + Precision);
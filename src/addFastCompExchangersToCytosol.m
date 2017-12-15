function [ComparisonModel] = addFastCompExchangersToCytosol(model,Exchangers,cytosolID)

%First: Add all exchangers (to the cytosol)
ComparisonModel = model;
metIDs = strcat(Exchangers(:,1),'[',cytosolID,']');
Coefs = diag([Exchangers{:,3}]);
%make the metIDs unique
[umets,order,source] = unique(metIDs);
reacCoefs = sparse(numel(umets),size(Exchangers,1));
for i = 1:numel(umets)
    reacCoefs(i,:) = sum(Coefs(source==i,:),1);
end
lbs = cell2mat(Exchangers(:,4));
ubs = cell2mat(Exchangers(:,5));
rxnIDs = strcat('Ex_', Exchangers(:,1), '[', Exchangers(:,2), ']');
ComparisonModel = addReactionBatch(model,rxnIDs,umets,reacCoefs,'lb',lbs,'ub',ubs);

end

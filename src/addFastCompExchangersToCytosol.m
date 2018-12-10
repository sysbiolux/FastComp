function [ComparisonModel] = addFastCompExchangersToCytosol(model,Exchangers,cytosolID)
% Add all Exchangers in the list to the cytosol, i.e. use the cytosol ID
% for all metabolite localisations
% 
% USAGE:
%    [ComparisonModel] = addFastCompExchangersToCytosol(model,Exchangers,cytosolID)
% 
% INPUTS: 
%    model:         A COBRA model structure. 
%    Exchangers:    List of Exchangers for the model including Demand and uptake reactions.  
%                   This is a n x 5 cell array with the first element being the
%                   metabolite(without compartment) the second element is the compartment ID
%                   The third element is the stoichiometric coefficient 
%                   The fourth and fifth elementes are the lower and upper bound respectively.
%    cytosolID:     The ID for all metabolites. commonly 'c'
%
% OUTPUTS:
% 
%    CoparisonModel:   The MOdel with all exchangers added to the cytosol
%
% Author:
%       Thomas Pfau


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
ComparisonModel = addMultipleReactions(model,rxnIDs,umets,reacCoefs,'lb',lbs,'ub',ubs);

end

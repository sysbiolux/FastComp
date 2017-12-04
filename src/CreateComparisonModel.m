function [ComparisonModel] = CreateComparisonModel(model,Exchangers,cytosolID)
%Add the exchangers to the model. 
%The model has to be able to carry positive flux through all reactions
%including the specified exchange reactions.
ComparisonModel = model;
for reac = 1:size(Exchangers)
    cmet = Exchangers{reac,1};    
    ComparisonModel = addReaction(ComparisonModel,['Ex_' cmet],'metaboliteList',{[cmet '[' cytosolID ']']},'stoichCoeffList',Exchangers{reac,3},'lowerBound',Exchangers{reac,4},'upperBound',Exchangers{reac,5});
end
end

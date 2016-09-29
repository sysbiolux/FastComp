function [updatedmodel,UpdatedExchangers,ComparisonModel] = CreateComparisonModelAndUpdateExchangersOrder(model,Exchangers,cytosolID)

%First: Add all exchangers (to the cytosol)
ComparisonModel = model;
for reac = 1:size(Exchangers,1)
    cmet = [Exchangers{reac,1} '[' cytosolID ']'];
    if isempty(find(ismember(model.mets,cmet)))
        %if the metabolite is not part of the model, we skip it.
        fprintf('%s was not part of the model, skipping the addition of the requested exchange reaction\n',cmet);
        continue
    end
    name = [ 'Ex_' cmet ];
    if isempty(find(ismember(ComparisonModel.rxns,name)))
        %Either add it
        ComparisonModel = addReaction(ComparisonModel,['Ex_' cmet],{cmet},-Exchangers{reac,3},0,0,max(max(ComparisonModel.ub),1000),0);
    else
        %Or adjust the bounds.
        pos = find(ismember(ComparisonModel.rxns,name));
        direc = full(ComparisonModel.S(ismember(ComparisonModel.mets,cmet),pos));
        if sign(Exchangers{reac,3}) ==  sign(direc)
            %we have already added the reaction, but with the opposite
            %directionality
            if ComparisonModel.lb(pos) == 0
                ComparisonModel.rev(pos) = 1;
                ComparisonModel.lb(pos) = min(min(ComparisonModel.lb),-1000);
            end
        end
    end

end

[~,ComparisonModel] = fastcc_ChangeModel(ComparisonModel,1);
%update the original model directionalities
updatedmodel = model;
updatedmodel.S = ComparisonModel.S(:,1:numel(model.rxns));

%And now, update the exchangers.
[imps,exps] = findImpsAndExpRxns(ComparisonModel);
imps = find(imps);
exps = find(exps);
ExchangeRs = unique([imps; exps]);

NewExchangers = {};
for i=1:numel(ExchangeRs)
    met = find(ComparisonModel.S(:,ExchangeRs(i)));
    commet = regexprep(ComparisonModel.mets(met),'\[[a-z]\]$','');
    comp = regexprep(ComparisonModel.mets(met),'.*\[([a-z])\]$','$1');
    NewExchangers(end+1,1:3) = {commet{1},comp{1},full(-sign(ComparisonModel.S(met,ExchangeRs(i))))};    
    if ComparisonModel.lb(ExchangeRs(i)) < 0
        NewExchangers(end+1,1:3) = {commet{1},comp{1},full(sign(ComparisonModel.S(met,ExchangeRs(i))))};
    end
end

UpdatedExchangers = cell(0,3);

exmets = unique(NewExchangers(:,1));


for i=1:numel(exmets)
    origpos = ismember(Exchangers(:,1),exmets{i});
    newpos = find(ismember(NewExchangers(:,1),exmets{i}));
    for cpos = 1:numel(newpos)
        dir = NewExchangers{newpos(cpos),3};
        matches = origpos & cellfun(@(x) x == dir, Exchangers(:,3));       
        UpdatedExchangers = [UpdatedExchangers ; Exchangers(matches,:)];
    end
end
        

end

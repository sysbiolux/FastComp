function [selUpt,selExp] = findImpsAndExpRxns(model)
    % Find exchange rxns
    selExc = (sum(model.S~=0) ==1)';

    %Uptake needs a reversible reaction with a lower bound 
    selUpt = full((model.lb < 0 & selExc & (model.rev == 1) & (sum(model.S) < 0)') | (model.ub > 0 & selExc & (sum(model.S) > 0)'));
    selExp = full((model.ub > 0 & selExc & (sum(model.S) < 0)') | (model.lb < 0 & (model.rev == 1) & selExc & (sum(model.S) > 0)'));
    
end
    
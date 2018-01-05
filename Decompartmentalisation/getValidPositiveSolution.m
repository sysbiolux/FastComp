function [modelMod,Problematic] = getValidPositiveSolution(model)

%first, adjust directionality for reactions that have an ub <= 0;
modelMod = flipRxnOrientation(model,model.rxns(model.ub <= 0));
modelMod.c(:) = 1;
sol = LP7(1:numel(model.rxns),model,1);
if sum(sol >=1) > numel(model.rxns)/3
    full = sol;
    sol = struct();
    sol.full = full;
else
    sol = optimizeCbModel(modelMod,'max');    
end
P = find(sol.full >= 1 | sol.full <= -1);
NP = setdiff(1:numel(modelMod.rxns),P);
flippedReactions = find(sol.full<= -1);
modelMod = flipRxnOrientation(modelMod,modelMod.rxns(sol.full <= -1));
modelTemp = modelMod;
modelTemp.lb(P) = 1;
modelTemp.c(:) = 0;
sol = optimizeCbModel(modelTemp,'max','one');
P = union(P,find(sol.v >= 1 | sol.v <= -1));
NP = setdiff(NP,P);
flippedReactions = union(flippedReactions, find(sol.v <= -1));
modelMod = flipRxnOrientation(modelMod,model.rxns(sol.v <= -1));
Problematic = [];
NP = setdiff(NP,P);
%now, lets add each reaction

while(numel(NP) ~= 0)    
    creac = NP(1);    
    sol = getSolForAllP(modelMod,P);
    P = find(sol.v >= 1 | sol.v <= -1);    
    flippedReactions = union(flippedReactions, find(sol.v <= -1));
    modelMod = flipRxnOrientation(modelMod,modelMod.rxns(sol.v <= -1));
    sol.v(sol.v < -1) = - sol.v(sol.v < -1);    
    if any(ismember(P,creac))
        NP = setdiff(NP,P);
        continue
    end
    if sol.v(creac) < 0
        modelModTemp = flipRxnOrientation(modelMod,modelMod.rxns(creac));
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelModTemp,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);
            modelMod = flipRxnOrientation(modelMod,modelMod.rxns(creac));
            flippedReactions = union(flippedReactions,creac);            
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end
    elseif sol.v(creac) > 0 
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelMod,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);                        
            continue;
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end
    else %Neither larger, nor smaller , i.e. its 0. %We should be able to add the solution.
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelMod,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);                        
            continue;
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end        
    end    
end

    
%Now if we have problematic reactions, than thats due to the 
%lets retry the problematic reactions

sol = LP7(1:numel(modelMod.rxns),modelMod,1);
P = find(sol >= 1);
NP = setdiff(1:numel(modelMod.rxns),P);
Problematic = [];

while(numel(NP) ~= 0)    
    creac = NP(1);    
    sol = getSolForAllP(modelMod,P);
    P = find(sol.v >= 1 | sol.v <= -1);    
    flippedReactions = union(flippedReactions, find(sol.v <= -1));
    modelMod = flipRxnOrientation(modelMod,modelMod.rxns(sol.v <= -1));
    sol.v(sol.v < -1) = - sol.v(sol.v < -1);    
    if any(ismember(P,creac))
        NP = setdiff(NP,P);
        continue
    end
    if sol.v(creac) < 0
        modelModTemp = flipRxnOrientation(modelMod,modelMod.rxns(creac));
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelModTemp,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);
            modelMod = flipRxnOrientation(modelMod,modelMod.rxns(creac));
            flippedReactions = union(flippedReactions,creac);            
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end
    elseif sol.v(creac) > 0 
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelMod,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);                        
            continue;
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end
    else %Neither larger, nor smaller , i.e. its 0. %We should be able to add the solution.
        Ptemp = union(P,creac);
        solTemp = getSolForAllP(modelMod,Ptemp);
        if solTemp.stat == 1
            P = union(P,creac);
            NP = setdiff(NP,creac);                        
            continue;
        else
            Problematic = union(Problematic,creac);
            NP = setdiff(NP,creac);
        end        
    end    
end

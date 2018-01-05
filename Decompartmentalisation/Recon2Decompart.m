%Decompartmentalisation
clearvars
Model = readCbModel('~/Work/Models/Human/Recon2.v04.mat');
%Remove constraints
Model.lb(Model.lb > 0) = 0;
Model.ub(Model.ub < 0) = 0;

%Increase bounds
Model.lb = Model.lb*1e3;
Model.ub= Model.ub*1e3;

FastCCModel = Model;
FastCCResult = fastcc(Model,1);

%Reduce the Model to a model that is internally consistent, i.e. all
%reactions can carry proper flux.
while ~ ( numel(FastCCResult) == numel(FastCCModel.rxns))
    fprintf('Removing %i reactions\n', numel(setdiff(FastCCModel.rxns,FastCCModel.rxns(FastCCResult))));
    FastCCModel = removeRxns(FastCCModel,setdiff(FastCCModel.rxns,FastCCModel.rxns(FastCCResult)));    
    FastCCResult = fastcc(FastCCModel,1);
end
FastCCModel = removeRxns(FastCCModel,setdiff(FastCCModel.rxns,FastCCModel.rxns(FastCCResult)));    

%Make all reactions able to carry positive flux, this is necessary, in
%order to be able to extract the Exchangers properly.
[PosModel,problems] = getValidPositiveSolution(FastCCModel);

%%
%clear everything except the fastcc results. 
clearvars -except FastCCResult Model FastCCModel PosModel
Model = PosModel;

%% First, we have to get the compartment information for all reactions
metComps = regexprep(Model.mets, '.*\[([a-z])\]$','$1');
reacComps = cellfun(@(x) unique(metComps(find(Model.S(:,ismember(Model.rxns,x))))),Model.rxns,'UniformOutput',false);

%%And Find the exchange reactions
[imps,exps] = findImpsAndExpRxns(Model);

ExchangerList = setupExchangers(Model,find(imps|exps),reacComps(imps|exps));

%% Now, we will have to set the localisation of the ETM and correct the ETM.
% There are a few SPECIAL reactions (i.e. the mitochondrial electron
% transfer chain) For these reactions we move every cytosolic hydrogen to a
% special hydrogen.
ETC = {'CYOR_u10m','L_LACDcm','NADH2_u10m','THD1m','CYOOm3','SULFOX','ATPS4m'};
Model = addMetabolite(Model,'h_ext[c]','External Proton','H');
for i = 1:numel(ETC)
    Model = changeRxnMets(Model,{'h[c]'},{'h_ext[c]'},ETC{i});
    reacComps{ismember(Model.rxns,ETC{i})} = {'m'};
end

%% Then we merge all metabolites.
NewModel = Model;

mets = regexprep(NewModel.mets,'\[.\]$','[c]'); %We merge all.
%replace all met names
NewModel.mets = mets;
[umets,ic,ia] = unique(mets);
for i= 1:numel(umets)    
    [pres,pos] = ismember(NewModel.mets,umets{i});       
    NewModel = mergeModelFieldPositions(NewModel,'mets',pres,{'metCharges',@(x) x(1); 'metFormulas',@(x) x(1)});
end

%% And finally, we A) remove all empty reactions and B) merge all non unique reactions (we don't care about the genes for now...)

%% A: 
toRemove = sum(abs(NewModel.S),1) == 0;
NewModel = removeRxns(NewModel,NewModel.rxns(toRemove));
reacComps(toRemove) = [];

%% B:
ubs = NewModel.ub;
lbs = NewModel.lb;
cs = NewModel.c;
i = 1;
%checkedReacs = {};
OrigReacs = containers.Map;
[imps,exps] = findImpsAndExpRxns(NewModel);
ExReactions = NewModel.rxns(imps|exps);
while i <= numel(NewModel.rxns)
    %Skip any exchange reaction
    if any(ismember(ExReactions,NewModel.rxns(i)))
        i = i+1;
        continue
    end
    %find all forward/reverse reactions that are the same.
    %checkedReacs{end+1} = NewModel.rxns{i};      
    reacMat = repmat(NewModel.S(:,i),1,numel(NewModel.rxns));
    sameForward = all(reacMat == NewModel.S,1);
    sameBackward = all(reacMat == -NewModel.S,1);
    lb = min([min(NewModel.lb(sameForward)),min(-NewModel.ub(sameBackward))]);
    ub = max([max(NewModel.ub(sameForward)),max(-NewModel.lb(sameBackward))]);
    NewModel.lb(i) = lb;
    NewModel.ub(i) = ub;
    reacComps(i) = {unique(vertcat(reacComps{sameForward},reacComps{sameBackward}))};
    NewModel = mergeModelFieldPositions(NewModel,'rxns',sameBackward | sameForward,{'S',@(x) x(:,1);...
                                        'grRules',@(x) {strjoin(setdiff(x,''), ' or ')};'rules',@(x) {strjoin(setdiff(x,''), ' | ')};...
                                        'rxnNames', @(x) x(1); 'rxns', @(x) x(1);...
                                        'rxnGeneMat',@(x) any(x,1) });
    
    %NewModel = removeRxns(NewModel,NewModel.rxns(setdiff(find(sameForward | sameBackward),i)));
    reacComps(setdiff(find(sameForward | sameBackward),i)) = [];        
    i = i+1;
end


%% Now, we have representative reactions, so we need to redo the metabolite merger.
A = fastcc(NewModel,1);
[Nimps,Nexps] = findImpsAndExpRxns(NewModel);

DecompModel = removeRxns(NewModel,setdiff(NewModel.rxns,NewModel.rxns(A)));

OrigComps= reacComps(A);
%% Set up the output

decompHuman = getValidPositiveSolution(DecompModel);
%And now, update the exchange reactions.
[imps,exps] = findImpsAndExpRxns(decompHuman); 
Exch = find(imps|exps);
Exchangers = ExchangerList;
for i = 1:numel(Exch)
    Exchangers(i,3:5) = {decompHuman.S(find(decompHuman.S(:,Exch(i))),Exch(i)),decompHuman.lb(Exch(i)),decompHuman.ub(Exch(i))};
end
decompHumanWOExchangers = removeRxns(decompHuman,decompHuman.rxns(imps|exps),'metFlag',false);
OrigComps(imps|exps) = [];
Recon_CompIDs = setdiff(unique(metComps),'c')';
comps = OrigComps;
Exchangers = ExchangerList;
save('Recon2ForFastComp','Recon_CompIDs', 'decompHuman', 'decompHumanWOExchangers', 'comps','Exchangers');

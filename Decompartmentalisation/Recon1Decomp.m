%Decompartmentalisation

Model = RECON1;
%Remove constraints
Model.lb(Model.lb > 0) = 0;
Model.ub(Model.ub < 0) = 0;

%Increase bounds
Model.lb(Model.lb == -999999) = -1000000;
Model.ub(Model.ub == 999999) = 1000000;

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

%%
%clear everything except the fastcc results. 
clearvars -except FastCCResult RECON1 Model FastCCModel
Model = FastCCModel;

%% First, we have to get the compartment information for all reactions
metComps = regexprep(Model.mets, '.*\[([a-z])\]$','$1');
reacComps = cellfun(@(x) unique(metComps(find(Model.S(:,ismember(Model.rxns,x))))),Model.rxns,'UniformOutput',false);

%% Now, we will have to set the localisation of the ETM and correct the ETM.
Comps = {'c';'e';'g';'l';'m';'n';'r';'x'};
% There are a few SPECIAL reactions (i.e. the mitochondrial electron
% transfer chain) For these reactions we move every cytosolic hydrogen to a
% special hydrogen.
ETC = {'ATPS4m','CYOOm3','CYOR_u10m','NADH2_u10m','THD1m'};
Model = addMetabolite(Model,'h_ext[c]','External Proton','H');
for i = 1:numel(ETC)
    Model = changeRxnMets(Model,{'h[c]'},{'h_ext[c]'},ETC{i});
    reacComps{ismember(Model.rxns,ETC{i})} = {'m'};
end
%There are also two replicate dead reactions.


%% Then we merge all metabolites.
NewModel = Model;

mets = regexprep(NewModel.mets,'\[[cglmnrx]\]$','[c]'); %All Except external
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
while i <= numel(NewModel.rxns)
    %find all forward/reverse reactions that are the same.
    %checkedReacs{end+1} = NewModel.rxns{i};      
    reacMat = repmat(NewModel.S(:,i),1,numel(NewModel.rxns));
    sameForward = all(reacMat == NewModel.S,1);
    sameBackward = all(reacMat == -NewModel.S,1);
    lb = min([min(NewModel.lb(sameForward)),min(-NewModel.ub(sameBackward))]);
    ub = max([max(NewModel.ub(sameForward)),max(-NewModel.lb(sameBackward))]);
    NewModel.lb(i) = lb;
    NewModel.ub(i) = ub;
    reacComps(i) = {unique([reacComps{sameForward},reacComps{sameBackward}])};
    
    NewModel = mergeModelFieldPositions(NewModel,'rxns',sameBackward | sameForward,{'S',@(x) x(:,1);...
                                        'grRules',@(x) {strjoin(setdiff(x,''), ' or ')};'rules',@(x) {strjoin(setdiff(x,''), ' | ')};...
                                        'rxnNames', @(x) x(1); 'rxns', @(x) x(1);...
                                        'rxnGeneMat',@(x) any(x,1) });
    
    %NewModel = removeRxns(NewModel,NewModel.rxns(setdiff(find(sameForward | sameBackward),i)));
    reacComps(setdiff(find(sameForward | sameBackward),i)) = [];        
    if any(ismember(reacComps{i},{'e'}))
        %THis ensures, that we do remove duplicate reactions. 
        reacComps{i} = {};            
    end  
    i = i+1;
end


%% Now, we have representative reactions, so we need to redo the metabolite merger.


A = fastcc(NewModel,1);
DecompModel = removeRxns(NewModel,setdiff(NewModel.rxns,NewModel.rxns(A)));
%% now, replace the localisation 
OrigComps= reacComps(A);
%% And then check again, that all reaction localisation makes sense.

Recon_Decomp = getValidPositiveSolution(DecompModel);
%Recon_Decomp = DecompModel;
Recon_CompIDs = {'g','l','m','n','r','x'};
Recon_NonExtReacs = find(~cellfun(@isempty ,OrigComps));
Recon_OrigLocs = OrigComps(Recon_NonExtReacs);
save('Recon1ForFastComp','Recon_CompIDs', 'Recon_Decomp', 'Recon_NonExtReacs', 'Recon_OrigLocs');
function [compmodel,CoreReactions,Transporters,nonLocReacSets,nonLocReacNames] = ...
MakeCompartmentsOneTransport(model,Complist,IntMets,ExtMets,ExtRxns,NonCompMetList,...
CompartmentIDs,GeneCompData,cytosolID,ExchangeReactions)
% Generate an extended model with unlocalised reactions in all compartments and localised reactions 
% in their specific compartment. 
%
% USAGE:
%    [compmodel,CoreReactions,Transporters,nonLocReacSets,nonLocReacNames] = ...
%               MakeCompartmentsOneTransport(model,Complist,IntMets,ExtMets,ExtRxns,NonCompMetList,...
%               CompartmentIDs,GeneCompartmentalisationData,cytosolID,ExchangeReactions)
%
% INPUTS: 
%    model:                 The compartmentilised model (i.e. the model with all
%                           non localised reactions in all compartments).
%    Complist:              Indicators of reaction localisations, cell array of cell arrays of strings indicating the compartments..
%    IntMets:               Indices of all Metabolites which are internal
%                           (i.e. not in a fixed external compartment)
%    ExtMets:               Indices of all Metabolites which are external
%                           (i.e. in a fixed external compartment)
%    NonCompMetList:        Names of all non Compartmentalised metabolites.
%                           includes both IntMets and ExtMets. 
%    CompartmentIDs:        The compartment ids (cell of strings),
%                           excluding the cytosol
%    GeneCompData:          As Complist but for Genes (can be empty).
%    SingleCompModel:       Boolean indicator, whether this is a single compartment
%                           model (true -> only cytosol), or has a fixed
%                           external compartment (false).
%    cytosolID:             Id for the cytosol.
%    ExchangeReactions:     List of Exchangers for the model including Demand and uptake reactions.  
%                           This is a n x 5 cell array with the first element being the
%                           metabolite(without compartment) the second element is the compartment ID
%                           The third element is the stoichiometric coefficient 
%                           The fourth and fifth elementes are the lower and upper bound respectively.
%
% OUTPUTS:
%
%    compmodel:             The Model with all transporters, and
%                           compartments.
%    CoreReactions:         Positions of the localised reactions in the
%                           Extended Model.                           
%    Transporters:          Positions of the transporters in the extended
%                           model.
%    nonLocReacSets:        The sets of non localised reactions. A double
%                           array of indices, with one row per non localised reaction indicating
%                           all positions of the reaction in the different
%                           compartments.
%    nonLocReacNames:       NAmes of all non localised reactions.
%
% .. Authors:
%       - Thomas Pfau 
%                               
             

warn = warning();%Deactivate Warnings for this function call, there will be plenty
%otherwise.
warning('off');

CoreReactions = [];
Transporters = [];
CompTransporters = {};
CompartmentIDs = [cytosolID CompartmentIDs];
NonLocReacs = unique(setdiff(find(cellfun(@length,Complist) == 0),ExtRxns));
nonLocReacNames = model.rxns(NonLocReacs);
%Make a copy of the original model
OrigModel = model;
nonLocReacSets = zeros(numel(NonLocReacs),numel(CompartmentIDs));
%create a model with only the external Reactions, all available genes and
%the respective GPRS etc pp.
%fprintf('The following reactions have been assigned as external in MakeCompartments\n');

%For a singleCompModel this should remove ALL reactions and metabolites, so
%we will have to rebuild them. 
%For an external + cytosol model, the external Reactions should be kept.
compmodel = removeRxns(model,setdiff(model.rxns,model.rxns(ExtRxns)));
ExchangedMetabolites = {};
for reac = 1:size(ExchangeReactions,1)
    cmet = [ExchangeReactions{reac,1} '[' ExchangeReactions{reac,2} ']'];
    ExchangedMetabolites = union(ExchangedMetabolites,cmet);
end
cytosolmetoffset = numel(compmodel.mets);
% find(ismember(compmodel.rxns,'6HTSTSTERONEte'))
%Create a list of Exchanged Metabolites


%fprintf('After Removing Exchange reactions, the model contains %i reactions and there are %i external reactions\n',numel(compmodel.rxns), numel(ExtRxns));
% save('Inmodel2')
emptyreac = 0;
metFields = getModelFieldsForType(model,'mets');
relMetFields = setdiff(metFields,{'S','mets','b','csense'}); % remove FBA fields.
metpars = {};
for metField = 1: numel(relMetFields)
    metpars{end+1} = relMetFields{metField};
    metpars{end+1} = model.(relMetFields{metField})(IntMets);
end

rxnFields = getModelFieldsForType(model,'rxns');
relRxnFields = setdiff(rxnFields,{'S','rxns','C', 'rxnGeneMat'}); % remove FBA fields.
rxnspars = {};
for rxnField = 1: numel(relRxnFields)
    rxnspars{end+1} = relRxnFields{rxnField};
    rxnspars{end+1} = model.(relRxnFields{rxnField})(NonLocReacs);
end

NonLocReacMatrix = OrigModel.S(IntMets,NonLocReacs);

%For each compartment:
%1. add all internal metabolites (Except those already present).
%2. Add All non Localised reactions.
%3. Add all localised Reactions

for c = 1:numel(CompartmentIDs)       
    %To do this, we will first add all metabolites to the new compartment            
    %Extend the S Matrix by IntMets Metabolites and NonLocReacs reactions    
    compMetNames = strcat(NonCompMetList(IntMets),['[' CompartmentIDs{c} ']']);    
    %check for presence
    pres = ismember(compMetNames,compmodel.mets);
    if any(pres)
        cmetpars = metpars;
        for i = 2:2:numel(metpars)
            cmetpars{i} = metpars{i}(~pres,:);
        end
    else
        cmetpars = metpars;
    end    
    compmodel = addMetaboliteBatch(compmodel,compMetNames(~pres),...
                cmetpars{:});                
    if numel(NonLocReacs) > 0
        compRxnNames = strcat(regexprep(model.rxns(NonLocReacs),'\([a-z]+\)',''),['(' CompartmentIDs{c} ')']);
        compmodel = addReactionBatch(compmodel,compRxnNames,compMetNames,NonLocReacMatrix, rxnspars{:});
    end
    %Now add the localised Reactions of this compartment    
    ExclusiveRxns = find(cellfun(@(x) not(isempty(find(ismember(x,CompartmentIDs{c})))),Complist));
    ExclusiveStoich = OrigModel.S(IntMets,ExclusiveRxns);
    exrxnspars = {};
    for rxnField = 1: numel(relRxnFields)
        exrxnspars{end+1} = relRxnFields{rxnField};
        exrxnspars{end+1} = model.(relRxnFields{rxnField})(ExclusiveRxns);
    end
    nRxnsBefore = numel(compmodel.rxns);
    if numel(ExclusiveRxns) > 0
        ExRxnNames = strcat(regexprep(model.rxns(ExclusiveRxns),'\([a-z]+\)',''),['(' CompartmentIDs{c} ')']);
        compmodel = addReactionBatch(compmodel,ExRxnNames,compMetNames,ExclusiveStoich, exrxnspars{:});
    end
    nRxnsAdded = numel(compmodel.rxns);
    CoreReactions = [CoreReactions, (nRxnsBefore+1):nRxnsAdded];
    fp = FormulaParser();
    for r=1:numel(ExclusiveRxns)
        %We might need to update the GPR rules.
        if ~all(cellfun(@isempty,GeneCompData)) && ~isempty(compmodel.rules{nRxnsBefore+r})
             CurrentRule = fp.parseFormula(compmodel.rules{nRxnsBefore+r});
             %We convert this to a DNF form and then we remove all clauses,
             %which have members not localised to this compartment.
             %And we will add all clauses, which have either localised
             %genes, or no genes localised to different compartments.
             dnfNode = CurrentRule.convertToDNF(); %the DNF node can either be a single AND nod, a literal node, or an OR node. 
             if isa(dnfNode, 'AndNode') || isa(dnfNode,'LiteralNode')
                genes = cellfun(@str2num, dnfNode.getLiterals());
                correctLoc = any(ismember( CompartmentIDs{c}, [GeneCompData{genes}])) || all(cellfun(@isempty, GeneCompData(genes)));
                if ~correctLoc % Otherwise everything is well.
                    compmodel.rules{nRxnsBefore+r} = '';
                    if isfield(compmodel,'grRules')
                        compmodel.grRules{nRxnsBefore+r} = '';
                    end
                    if isfield(compmodel,'rxnGeneMat')
                        compmodel.rxnGeneMat(nRxnsBefore+r,:) = false;
                    end
                end
             end
             if isa(dnfNode,'OrNode')
                 %Now, we have to walk over all children and add them as
                 %above, if they fit.
                 rule = '';
                 NewNode = OrNode();
                 for cchild = 1:numel(dnfNode.children)
                     childNode = dnfNode.children(cchild);
                     genes = cellfun(@str2num, childNode.getLiterals());
                     correctLoc = any(ismember( CompartmentIDs{c}, [GeneCompData{genes}])) || all(cellfun(@isempty, GeneCompData(genes)));
                     if correctLoc % extend the rule
                         rule = strcat(rule,' | ',childNode.toString());
                         NewNode.addChild(childNode);
                     end
                 end
                 compmodel.rules{nRxnsBefore+r} = rule;
                 if isfield(compmodel,'grRules')
                    compmodel = creategrRulesField(compmodel,nRxnsBefore+r);
                 end
                 if isfield(compmodel,'rxnGeneMat')
                     
                     compmodel.rxnGeneMat(nRxnsBefore+r,:) = false;
                     literals = NewNode.getLiterals();
                     if ~isempty(literals)
                         literals = cellfun(@str2num, literals);
                         compmodel.rxnGeneMat(nRxnsBefore+r,literals) = true;
                     end
                 end
             end             
        end
       
    end
    %Add Transporters for all non external Metabolites
    %skip this step if we are in the cytosol
    if not(strcmp(CompartmentIDs{c},cytosolID))
        nRxnsStart = numel(compmodel.rxns);
        %Look up the positions of the metabolites 
        metabPositions = ismember(compmodel.mets,compMetNames);
        ExMetsForComp = sum(compmodel.S(:,:)~=0,2)> 0 & metabPositions;
        compMetNames = compmodel.mets(ExMetsForComp);
        cytMetNames = regexprep(compMetNames,'\[[^\[]+\]$',['[' cytosolID ']']);
        metIDs = [cytMetNames;compMetNames];
        nMets = numel(cytMetNames);
        stoich = [speye(nMets); -speye(nMets)];
        lbs = -abs(max([model.lb ; model.ub])) * ones(nMets,1);
        ubs = -lbs;
        metNames = regexprep(compMetNames,'\[[^[]+\]$','');
        rxnNames = strcat(metNames, {' Transport'});
        rxnIDs = strcat('R_T', upper(CompartmentIDs{c}), 'C_', metNames);
        compmodel = addReactionBatch(compmodel,rxnIDs,metIDs,stoich,'lb',lbs,'ub',ubs,'rxnNames',rxnNames);
        nRxnsEnd = numel(compmodel.rxns);
        CompTransporters{end+1} = (nRxnsStart+1):nRxnsEnd;
        Transporters = [Transporters, (nRxnsStart+1):nRxnsEnd];                    
    end
  
end
%Now we need to perform some cleanup:
%combine all Metabolites with the same name - There should be none in a 
%singlecompmodel but in an external/cytosol model there are those involved in transporters.
%to do this: We need to: 
% -> Combine all respective rows (delete the latter)
% -> Remove the respective mets, metNames, metFormula, b rows
[uni, idx_last, idx] = unique(compmodel.mets);

nonunique = find(accumarray(idx,1) > 1);
if numel(nonunique) > 0
    fprintf('There are %i non unique metabolites\n', numel(nonunique))
end

for i = 1:numel(nonunique)
    %get all positions of these metabolite
    
    nonuniquemets = find(ismember(compmodel.mets,uni(nonunique(i))));    
    for metpos = 2:numel(nonuniquemets)       
        compmodel.S(nonuniquemets(1),:) = compmodel.S(nonuniquemets(1),:)+compmodel.S(nonuniquemets(metpos),:);
    end
    %We adjusted the S matrix, lets remove the surplus stuff.
    compmodel = removeMetabolites(compmodel,nonuniquemets(2:end));       
end
%     disp('After removing non unique metabolites the size of S is')
%     size(compmodel.S)
% %now clean up transporters that are dead ends.
% save PreExchange
%Also add the Exchange reactions. These reactions will also be part of the
%Core!
% size(compmodel.rxns);
nRxnsStart = numel(compmodel.rxns);
metIDs = strcat(ExchangeReactions(:,1),'[',ExchangeReactions(:,2) ,']');
Coefs = diag([ExchangeReactions{:,3}]);
%make the metIDs unique
[umets,order,source] = unique(metIDs);
reacCoefs = sparse(numel(umets),size(ExchangeReactions,1));
for i = 1:numel(umets)
    reacCoefs(i,:) = sum(Coefs(source==i,:),1);
end
lbs = cell2mat(ExchangeReactions(:,4));
ubs = cell2mat(ExchangeReactions(:,5));
rxnIDs = strcat('Ex_', ExchangeReactions(:,1), '[', ExchangeReactions(:,2), ']');
compmodel = addReactionBatch(compmodel,rxnIDs,umets,reacCoefs,'lb',lbs,'ub',ubs);

nRxnsEnd = numel(compmodel.rxns);
CoreReactions = [CoreReactions, (nRxnsStart+1):nRxnsEnd];


CoreReactions = unique(CoreReactions);
% remove all Metabolites which are not present in any reaction (those were
% added temporarily).
compmodel = removeMetabolites(compmodel,compmodel.mets(sum(compmodel.S~=0,2) == 0));

deadendmets = {};
% Get all metabolites only involved in exactly one reaction.
% Those metabolites are obviously dead.
% These metabolites can arise by metabolites allocated to the different 
% compartments, which we can't rremove until all compartments are set up,
% as we otherwise wouldn't know whether they might be transported through
% the cytosol to a different compartment.
deadendmets = compmodel.mets(sum(compmodel.S~=0,2) <= 1);


while numel(deadendmets) > 0
    rs_to_remove = findRxnsFromMets(compmodel,deadendmets);
    rs_pos = find(ismember(compmodel.rxns,rs_to_remove));
    transp_names = compmodel.rxns(Transporters);
    core_names = compmodel.rxns(CoreReactions);
    compmodel = removeRxns(compmodel,rs_to_remove,false);
    Transporters = find(ismember(compmodel.rxns,transp_names));
    CoreReactions = find(ismember(compmodel.rxns,core_names));
    deadendmets = compmodel.mets(sum(compmodel.S~=0,2) <= 1);    
end

if numel(NonLocReacs) > 0
    %We will change this to a cell array of double arrays, one for each
    %reaction.
    for c = 1:numel(CompartmentIDs)
        compRxnNames = strcat(regexprep(model.rxns(NonLocReacs),'\([a-z]+\)',''),['(' CompartmentIDs{c} ')']);
       [~,pos] = ismember(compRxnNames,compmodel.rxns);
        nonLocReacSets(:,c) = pos;
    end    
end

CoreReactions = columnVector(CoreReactions);
%Reactivate warnings.
warning(warn);
end


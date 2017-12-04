function [compmodel,CoreReactions,Transporters,nonLocReacSets,nonLocReacNames] = ...
MakeCompartmentsOneTransport(model,Complist,IntMets,ExtMets,ExtRxns,NonCompMetList,...
CompartmentIDs,GeneCompartmentalisationData,exclusivegenes,cytosolID,ExchangeReactions)

%Input: 
% model                         without any exchange Reactions 
% Complist                      A list of compartmentalisation information; the length
%                               has to be the same as the length of the model.rxns, field. Empty cells
%                               indicate that there is no information and the reactions are supposed to
%                               be available in all compartments and part of the nonLocReacSet
% IntMets                       A double array indicating all internal metabolites 
% ExtMets                       A double array indicating all External metabolites (for
%                               models with cytosol and external compartment
% ExtRxns                       A double array indicating all External reactions (for
%                               models with cytosol and external compartment
% NonCompMetList                A List of metabolite names without localisation
%                               information
% CompartmentIDs                All compartments that shall eb created, without the
%                               cytosol
% GeneCompartmentalisationData  Compartmentalisation Information of genes,
%                               if any
% exclusivegene                 Indicator, whether Genes should only be
%                               available in their respective compartments
% cytosolID                     ID of the Cytosol
% ExchangeReactions             List of ExchangeReactions to add.


CoreReactions = [];
Transporters = [];
CompTransporters = {};
CompartmentIDs = [cytosolID CompartmentIDs];
NonLocReacs = unique(setdiff(find(cellfun(@length,Complist) == 0),ExtRxns));
nonLocReacNames = model.rxns(NonLocReacs);
%Make a copy of the original model
OrigModel = model;
nonLocReacSets = [];
%create a model with only the external Reactions, all available genes and
%the respective GPRS etc pp.
%fprintf('The following reactions have been assigned as external in MakeCompartments\n');

%save('Temp.mat','ExtRxns');

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
    metpars{end+1} = compmodel.(relMetFields{metField})(IntMets);
end

rxnFields = getModelFieldsForType(model,'rxns');
relRxnFields = setdiff(metFields,{'S','rxns','C', 'rxnGeneMat'}); % remove FBA fields.
rxnspars = {};
for rxnField = 1: numel(relRxnFields)
    rxnspars{end+1} = relRxnFields{rxnField};
    rxnspars{end+1} = compmodel.(relMetFields{metField})(NonLocReacs);
end

NonLocReacMatrix = OrigModel.S(IntMets,NonLocReacs);

for c = 1:numel(CompartmentIDs)       
%     fprintf('There are currently %i reactions and %i gprs\n',length(compmodel.rxns),length(compmodel.grRules));
%     disp('The Current size of S is')
%     size(compmodel.S)
%     disp('Extending S')
%     fprintf('Adding %i metabolite rows and %i non localised reaction columns for compartment %s\n',size(model.S,1)-numel(ExtMets),numel(NonLocReacs),CompartmentIDs{c});    
%     fprintf('There are a total of %i metabolites with %i external metabolites, while S has a size of %i\n',numel(compmodel.mets),numel(ExtMets),size(compmodel.S,1))
%     fprintf('We add a total of %i rows and %i columns\n',size(compmodel.S,1)+size(model.S,1)-numel(ExtMets)-size(compmodel.S,1)+1,size(compmodel.S,2)+numel(NonLocReacs)-size(compmodel.S,2)+1)
%     fprintf('We add From Row %i to %i and col %i to %i\n',size(compmodel.S,1)+size(model.S,1)-numel(ExtMets),size(compmodel.S,1)+1,size(compmodel.S,2)+numel(NonLocReacs),size(compmodel.S,2)+1)
    %duplicate the non localised part of the S matrix
    %To do this, we will first add all metabolites to the new compartment    
    [compmetcount,compreaccount] = size(compmodel.S);    
    %Extend the S Matrix by IntMets Metabolites and NonLocReacs reactions
    compMetNames = strcat(NonCompMetList(IntMets),['[' CompartmentIDs{c} ']']);
    compmodel = addMetaboliteBatch(compmodel,compMetNames,...
                metpars{:});
    %And set the new bits to the old S-Matrix for those entries
    
    compmodel.S(compmetcount+1:compmetcount+numel(IntMets),compreaccount+1:compreaccount+numel(NonLocReacs)) = ...
        OrigModel.S(IntMets,NonLocReacs);       
    
    %if we have Non ocalised Reactions
    if numel(NonLocReacs) > 0
        compRxnNames = strcat(regexprep(model.rxns(NonLocReacs),'\([a-z]+\)',''),['(' CompartmentIDs{c} ')'])]
        compmodel = addReactionBatch(compmodel,compRxnNames,compMetNames,NonLocReacMatrix, rxnspars{:});
    end
    %Now add the localised Reactions of this compartment
    ExclusiveRxns = find(cellfun(@(x) not(isempty(find(ismember(x,CompartmentIDs{c})))),Complist));
    for r=1:numel(ExclusiveRxns)
        %add this as a core reaction.         
        CoreReactions(end+1,1) = numel(compmodel.rxns)+1;
        compmodel = addReaction(compmodel,strcat(regexprep(model.rxns{ExclusiveRxns(r)},'\([a-z]+\)',''),['(' CompartmentIDs{c} ')']),...
                                'lowerBound',model.lb(ExclusiveRxns(r)), 'upperBounds', model.ub(ExclusiveRxns(r)),...
                                'stoichCoeffList',model.S(IntMets,ExclusiveRxns(r)),...
                                'metaboliteList',compmodel.mets([1:numel(IntMets)] + compmetcount),...
                                'rxnName', [model.rxnNames{ExclusiveRxns(r)} '(' CompartmentIDs{c} ')']);                                    
        %now this gets interesting.... Since here we can assign gene
        %localisation.
        %First get the genes involved with this reaction        
        genelist = findGenesFromRxns(model,model.rxns{ExclusiveRxns(r)});
        genes = genelist{1};
        %disp('Assigning genes')
        if not(isempty(genes)) && ~isempty(GeneCompartmentalisationData)
            genepos = find(ismember(model.genes,genes));
            %Now, check which ones are assigned to the current compartment.
            %and which ones to other
            correctloc = {};
            differentloc = {};
            if exclusivegenes
                %This creates Gene Associations which only contain clauses
                %that have a relation to the assigned localisation
               for gene = 1:numel(genes)
                    if numel(GeneCompartmentalisationData{genepos(gene)}) > 0
                        if not(isempty(find(ismember(GeneCompartmentalisationData{genepos(gene)},CompartmentIDs(c)),1)))
                            correctloc(end+1) = genes(gene);
                        else
                            differentloc(end+1) = genes(gene);
                        end
                    end
                end
            else    
                %otherwise just add all non localised as correct as well.
                for gene = 1:numel(genes)
                    
                    if numel(GeneCompartmentalisationData{genepos(gene)}) > 0
                        if not(isempty(find(ismember(GeneCompartmentalisationData{genepos(gene)},CompartmentIDs(c)),1)))
                            correctloc(end+1) = genes(gene);
                        else
                            differentloc(end+1) = genes(gene);
                        end
                    else
                        correctloc(end+1) = genes(gene);
                    end
                end
            end

            %now, parse the GPR and extract the correct GPRs
            FP = FormulaParser();
            Formula = FP.parseFormula(model.rules{ExclusiveRxns(r)});
            Formula = Formula.convertToDNF();
            FormString = Formula.toString(0);
            Formula.reduce();
            %reduce the Formula.            
            while not(strcmp(FormString,Formula.toString(0)))
                FormString = Formula.toString(0);
                Formula.reduce();
            end
            Formula.toString(1);
            %The New Formula starts with an OR Node
            %disp('Creating Compartment specific Form')
            NewFormula = OrNode();
            %Add all clauses containing a correctly localised gene
            for gene = 1:numel(correctloc)
                if numel(Formula.children) > 0
                    for cid = 1:numel(Formula.children)
                        child = Formula.children(cid);
                        if child.contains(correctloc{gene})
                            NewFormula.addChild(child);
                        end
                    end
                else
                    if Formula.contains(correctloc{gene})
                        NewFormula.addChild(Formula);
                    end
                end
            end
            childstodel = [];
            %and remove from these clauses all which show differently
            %localised genes
            for gene = 1:numel(differentloc)
                for cid = 1:numel(NewFormula.children)
                    child = NewFormula.children(cid);
                    if child.contains(differentloc{gene})
                        childstodel(end+1) = cid;
                    end
                end
            end
            NewFormula.children(childstodel) = [];
            %Create the new GPR String
            GPRString = NewFormula.toString(0);
            NewFormula.reduce();
            %And reduce it to a DNF Form String
            while not(strcmp(GPRString,NewFormula.toString(0)))
                GPRString = NewFormula.toString(0);
                NewFormula.reduce();
            end
            NewFormula.removeDNFduplicates();
            GPRString = NewFormula.toString(1);            
            %Extend the matrix and set all involved genes to 1.
            % Also create the rules string by replacining the gene by its
            % position.            
            compmodel.rules{end+1} = GPRString;            
            %At the end we have to update the grRules and rxnGeneMat
            %fields.
        else
            %there are no assigned genes with this reaction so just
            compmodel.rules{end+1} = '';
            compmodel.grRules{end+1} = '';
            compmodel.rxnGeneMat(end+1,:) = 0;                        
        end
       
    end
    %Add Transporters for all non external Metabolites
    %skip this step if we are in the cytosol
    if not(strcmp(CompartmentIDs{c},cytosolID))
        nRxnsStart = numel(compmodel.rxns);
        
        CompTransporters{end+1} = [];
        ExMetsForComp = sum(compmodel.S(compmetcount+(1:numel(IntMets)),:)~=0,2)> 0;
        cytMetNames = compmodel.mets(cytosolmetoffset + find(ExMetsForComp));
        compMetNames = compmodel.mets(compmetcount+ find(ExMetsForComp));
        metIDs = [cytMetNames;compMetNames];
        nMets = numel(cytMetNames);
        stoich = [speye(nMets); -speye(nMets)];
        lbs = -abs(max([model.lb ; model.ub])) * ones(nMets,1);
        ubs = -lbs;
        metNames = NonCompMetList(IntMets(ExMetsForComp));
        rxnNames = strcat(metNames, {' Transport'});
        rxnIDs = strcat('R_T', upper(CompartmentIDs{c}), 'C_', metNames);
        compmodel = addReactionBatch(compmodel,rxnIDs,metIDs,stoich,'lb',lbs,'ub',ubs,'rxnNames',rxnNames);
        nRxnsEnd = numel(compmodel.rxns);
        CompTransporters{end} = (nRxnsStart+1):nRxnsEnd;
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
for reac = 1:size(ExchangeReactions,1)
    cmet = [ExchangeReactions{reac,1} '[' ExchangeReactions{reac,2} ']'];
    name = [ 'Ex_' cmet ];
    compmodel = addReaction(compmodel,['Ex_' cmet],'metaboliteList',{cmet},...
        'stoichCoeffList',ExchangeReactions{reac,3},'lowerBound',ExchangeReactions{reac,4},'upperBound',ExchangeReactions{reac,5});
    CoreReactions(end+1) = numel(compmodel.rxns);
end

CoreReactions = unique(CoreReactions);
nonlocreacs = compmodel.rxns(nonLocReacSets);

deadendmets = {};
%Get all metabolites only involved in exactly one reaction.
%Those metabolites are obviously dead.
for i = 1:numel(compmodel.mets)
met = i;
involved = sum(compmodel.S(met,:) ~= 0);
if involved <= 1
deadendmets{end+1} = compmodel.mets{met};
end
end
% save('Compmodel','compmodel')

while numel(deadendmets) > 0
    rs_to_remove = findRxnsFromMets(compmodel,deadendmets);
    rs_pos = find(ismember(compmodel.rxns,rs_to_remove));
    transp_names = compmodel.rxns(Transporters);
    core_names = compmodel.rxns(CoreReactions);
    compmodel = removeRxns(compmodel,rs_to_remove);
    Transporters = find(ismember(compmodel.rxns,transp_names));
    CoreReactions = find(ismember(compmodel.rxns,core_names));
    deadendmets = {};
    for i = 1:numel(compmodel.mets)
        met = i;
        involved = sum(compmodel.S(met,:) ~= 0);
        if involved <= 1
            deadendmets{end+1} = compmodel.mets{met};
        end
    end    
end

%And now, produce a Consistent Model - Since this is only necessary for FC
%runs we do this at another place and provide MO with a full network.
% [A,compmodel,toRemove] = fastcc_ChangeModel(compmodel,1);
% 
% rs_to_remove = compmodel.rxns(toRemove);    
% transp_names = compmodel.rxns(Transporters);
% core_names = compmodel.rxns(CoreReactions);
% compmodel = removeRxns(compmodel,rs_to_remove);
% Transporters = find(ismember(compmodel.rxns,transp_names));
% CoreReactions = find(ismember(compmodel.rxns,core_names));

if numel(nonlocreacs) > 0
    %We will change this to a cell array of double arrays, one for each
    %reaction.
    nonLocReacCellArrays = cell(numel(NonLocReacs),1);
    for i = 1:numel(nonLocReacNames)
        reacpositions = zeros(1,numel(CompartmentIDs));
        for c = 1:numel(CompartmentIDs)
            reacname = [nonLocReacNames{i} , '(' CompartmentIDs{c} ')'];
            reacpos = find(ismember(compmodel.rxns,reacname));
            if ~isempty(reacpos)
                reacpositions(c) = reacpos;
            else
                reacpositions(c) = 0;
            end
        end
        nonLocReacCellArrays{i} = reacpositions;
    end
    %save('Test','nonLocReacCellArrays')
    nonLocReacSets = cell2mat(nonLocReacCellArrays);
end

%save('compafterRemove.mat','compmodel','nonLocReacSets');
end


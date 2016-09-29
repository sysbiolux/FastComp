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
% find(ismember(model.rxns,'6HTSTSTERONEte'))
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
    compmodel.S(compmetcount+1:compmetcount+numel(IntMets),:) = 0;
    compmodel.S(:,compreaccount+1:compreaccount+numel(NonLocReacs)) = 0;
    %And set the new bits to the old S-Matrix for those entries
    compmodel.S(compmetcount+1:compmetcount+numel(IntMets),compreaccount+1:compreaccount+numel(NonLocReacs)) = ...
        OrigModel.S(IntMets,NonLocReacs);
    %also duplicate the corresponding metabolite rows.
    %disp('Extending Metabolites')
    compmodel.mets = [compmodel.mets ; strcat(NonCompMetList(IntMets),['[' CompartmentIDs{c} ']'])];
    compmodel.metNames = [compmodel.metNames ; model.metNames(IntMets)];
    compmodel.metFormulas = [compmodel.metFormulas ; model.metFormulas(IntMets)];
    compmodel.b = [compmodel.b ; model.b(IntMets)];
    if (isfield(model,'metChEBIID'))
        compmodel.metChEBIID = [compmodel.metChEBIID ; model.metChEBIID(IntMets)];
    end
    if isfield(model,'metKEGGID')
        compmodel.metKEGGID = [compmodel.metKEGGID ; model.metKEGGID(IntMets)];
    end
    if isfield(model,'metPubChemID')
        compmodel.metPubChemID = [compmodel.metPubChemID ; model.metPubChemID(IntMets)];  
    end
    if isfield(model,'metInChIString')
        compmodel.metInChIString = [compmodel.metInChIString ; model.metInChIString(IntMets)];
    end
    if isfield(model,'metCharge')
        compmodel.metCharge =   [compmodel.metCharge ; model.metCharge(IntMets)];
    end
    
    
    
    %if we have Non ocalised Reactions
    if numel(NonLocReacs) > 0
        % and the gene relationships (which are duplicates for each of this, as
        % localisation is unknown)
        %disp('Extending GPRs')
        compmodel.rules = [compmodel.rules ; model.rules(NonLocReacs)];
        compmodel.grRules = [compmodel.grRules; model.grRules(NonLocReacs)];
        compmodel.rxnGeneMat(size(compmodel.rxnGeneMat,1)+1:size(compmodel.rxnGeneMat,1)+numel(NonLocReacs),:) = model.rxnGeneMat(NonLocReacs,:);
        
        
        %and the reaction rows
        %disp('Extending Reactions')
        crxns = numel(compmodel.rxns);
        %add the positions of the nonLocalised Reactions
        nonLocReacSets(:,end+1) = [1:numel(NonLocReacs)] + crxns;
        %And update the remainig fields
        compmodel.rxns = [compmodel.rxns ; strcat(regexprep(model.rxns(NonLocReacs),'\([a-z]+\)',''),['(' CompartmentIDs{c} ')'])];
        compmodel.lb = [compmodel.lb ; model.lb(NonLocReacs)];
        compmodel.ub = [compmodel.ub ; model.ub(NonLocReacs)];
        compmodel.rev = [compmodel.rev ; model.rev(NonLocReacs)];
        compmodel.rxnNames = [compmodel.rxnNames ; model.rxnNames(NonLocReacs)];
        compmodel.c = [compmodel.c ; model.c(NonLocReacs)];
        if isfield(model,'subSystems')
            compmodel.subSystems = [compmodel.subSystems ; model.subSystems(NonLocReacs)];
        end
        if isfield(model,'confidenceScores')
            compmodel.confidenceScores = [compmodel.confidenceScores ; model.confidenceScores(NonLocReacs)];
        end
        if isfield(model,'rxnNotes')
            compmodel.rxnNotes = [compmodel.rxnNotes ; model.rxnNotes(NonLocReacs)];
        end
        if isfield(model,'rxnECNumbers')
            compmodel.rxnECNumbers = [compmodel.rxnECNumbers ; model.rxnECNumbers(NonLocReacs)];
        end
        if isfield(model,'rxnReferences')
            compmodel.rxnReferences = [compmodel.rxnReferences ; model.rxnReferences(NonLocReacs)];
        end
        if (isfield(model, 'rxnConfidenceScores'))
            compmodel.rxnConfidenceScores = [compmodel.rxnConfidenceScores ; model.rxnConfidenceScores(NonLocReacs)];;
        end
        if (isfield(model, 'rxnKeggID'))
            compmodel.rxnKeggID = [compmodel.rxnKeggID ; model.rxnKeggID(NonLocReacs)];
        end
        if (isfield(model, 'rxnsboTerm'))
            compmodel.rxnsboTerm = [compmodel.rxnsboTerm ; model.rxnsboTerm(NonLocReacs)];
        end
    else
        if (strcmp(CompartmentIDs{c},cytosolID)) && (size(compmodel.S,2) == 1)
            %Now we got a problem! we have created an empty reaction and we
            %need to fix this.
            emptyreac = 1;
        end
            
                
    end
    %Now add the localised Reactions of this compartment
    ExclusiveRxns = find(cellfun(@(x) not(isempty(find(ismember(x,CompartmentIDs{c})))),Complist));
    for r=1:numel(ExclusiveRxns)
        %add this as a core reaction.         
        CoreReactions(end+1,1) = numel(compmodel.rxns)+1;
        compmodel.rxns{end+1,1} = strcat(regexprep(model.rxns{ExclusiveRxns(r)},'\([a-z]+\)',''),['(' CompartmentIDs{c} ')']);
        compmodel.lb(end+1,1) = model.lb(ExclusiveRxns(r));
        compmodel.ub(end+1,1) = model.ub(ExclusiveRxns(r));
        compmodel.rev(end+1,1) = model.rev(ExclusiveRxns(r));
        compmodel.rxnNames{end+1,1} = [model.rxnNames{ExclusiveRxns(r)} '(' CompartmentIDs{c} ')'];
        compmodel.c(end+1,1) = model.c(ExclusiveRxns(r));        
        if emptyreac
            compmodel.S([1:numel(IntMets)] + compmetcount,end) = model.S(IntMets,ExclusiveRxns(r));        
            emptyreac = 0;
        else
            compmodel.S([1:numel(IntMets)] + compmetcount,end+1) = model.S(IntMets,ExclusiveRxns(r));        
        end
        if isfield(model,'subSystems')
            compmodel.subSystems{end+1,1} = model.subSystems{ExclusiveRxns(r)};
        end
        if isfield(model,'confidenceScores')
            compmodel.confidenceScores{end+1,1} = model.confidenceScores{ExclusiveRxns(r)};
        end        
        if isfield(model,'rxnNotes')
            compmodel.rxnNotes{end+1,1} = model.rxnNotes{ExclusiveRxns(r)};
        end
        if isfield(model,'rxnECNumbers')
            compmodel.rxnECNumbers{end+1,1} = model.rxnECNumbers{ExclusiveRxns(r)};
        end               
        if isfield(model,'rxnReferences')
            compmodel.rxnReferences{end+1,1} = model.rxnReferences{ExclusiveRxns(r)};
        end
        if (isfield(model, 'rxnConfidenceScores'))
            compmodel.rxnConfidenceScores{end+1,1} = model.rxnConfidenceScores{ExclusiveRxns(r)};
        end
        if (isfield(model, 'rxnKeggID'))
            compmodel.rxnKeggID{end+1,1} = model.rxnKeggID{ExclusiveRxns(r)};
        end
        if (isfield(model, 'rxnsboTerm'))
            compmodel.rxnsboTerm{end+1,1} = model.rxnsboTerm{ExclusiveRxns(r)};
        end
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
            %disp('Creating Formula')
            FP = FormulaParser();
            Formula = FP.parseFormula(model.grRules{ExclusiveRxns(r)});
            %disp('Creating DNNF Form')
            Formula = Formula.convertToDNF();
            FormString = Formula.toString(0);
            Formula.reduce();
            %reduce the Formula.
            %disp('Reducing Formula')
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
            NewFormula.toString(1);
            %disp('Removing non compartment genes')
            for gene = 1:numel(differentloc)
                for cid = 1:numel(NewFormula.children)
                    child = NewFormula.children(cid);
                    if child.contains(differentloc{gene})
                        childstodel(end+1) = cid;
                    end
                end
            end
            NewFormula.children(childstodel) = [];
            NewFormula.toString(1);
            %Create the new GPR String
            GPRString = NewFormula.toString(0);
            NewFormula.reduce();
            %And reduce it to a DNF Form String
            while not(strcmp(GPRString,NewFormula.toString(0)))
                GPRString = NewFormula.toString(0);
                NewFormula.reduce();
            end
            NewFormula.removeDNFduplicates();
            GPRString = NewFormula.toString(0);
            rulestring = NewFormula.toString(1);
            %Extend the matrix and set all involved genes to 1.
            % Also create the rules string by replacining the gene by its
            % position.
            compmodel.rxnGeneMat(end+1,:) = 0;
            for gid = 1:numel(genes)
                gene = genes{gid};
                oldrulesstring = rulestring;
                rulestring = strrep(rulestring,['(' gene ')'],['(x(' num2str(find(ismember(model.genes,gene))) '))']);
                if not(strcmp(oldrulesstring,rulestring))
                    %only if we changed anything it is necessary to create
                    %a gene reaction link
                    compmodel.rxnGeneMat(end,find(ismember(model.genes,gene))) = 1;
                end
            end
            compmodel.rules{end+1} = rulestring;
            compmodel.grRules{end+1} = GPRString;
            
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
        CompTransporters{end+1} = [];
        for m=1:numel(IntMets)
            %Check whether this metabolite is used in any Reaction,
            %otherwise we can simply "not add" the transporter.
            present = find(compmodel.S(compmetcount + m,:)~= 0);
            %If its not yet in the compartment and not exchanged, than we
            %can ignore it!
            if isempty(present) && isempty(intersect(ExchangedMetabolites,compmodel.mets(compmetcount + m)))
                %There are no reactions which could carry any flux so just
                %skip this metabolite
                
                continue;
            end
            %fprintf('Compartment %s and Metabolite %s\n',CompartmentIDs{c},model.mets{m})
            %disp('Adding info')
            compmodel.rxns{end+1} = ['R_T' upper(CompartmentIDs{c}) 'C_' NonCompMetList{IntMets(m)}];
            compmodel.rules{end+1} = '';
            compmodel.grRules{end+1} = '';
            %Take the absolute largest values as upper and lower bounds
            compmodel.lb(end+1) = -abs(max([model.lb ; model.ub]));
            compmodel.ub(end+1) = abs(max([model.lb ; model.ub]));
            compmodel.rev(end+1) = 1;
            compmodel.rxnNames{end+1} = [NonCompMetList{IntMets(m)} ' Transport'];
            compmodel.c(end+1) = 0;
            if isfield(model,'subSystems')
                compmodel.subSystems{end+1} = 'Transport internal for FastCompart';
            end
            if isfield(model,'confidenceScores')
                compmodel.confidenceScores{end+1} = '';
            end
            
            if isfield(model,'rxnNotes')
                compmodel.rxnNotes{end+1} = 'Automatically added  Transporter';
            end
            if isfield(model,'rxnECNumbers')
                compmodel.rxnECNumbers{end+1} = '';
            end                        
            if isfield(model,'rxnReferences')
                compmodel.rxnReferences{end+1} = '';
            end

            if (isfield(model, 'rxnConfidenceScores'))
                compmodel.rxnConfidenceScores{end+1} = '';
            end            
            if (isfield(model, 'rxnKeggID'))
                compmodel.rxnKeggID{end+1} = '';
            end
            if (isfield(model, 'rxnsboTerm'))
                compmodel.rxnsboTerm{end+1} = '';
            end

            %and now comes the hard part, updating the Stoichiometric matrix
            %disp('Modifingy matrices')
            compmodel.rxnGeneMat(end+1,:) = 0;
            compmodel.S(:,end+1) = 0;
            compmodel.S(cytosolmetoffset + m,end) = 1; % consumes one cytosolic compound
            compmodel.S([compmetcount + m],end) = -1;
            % produces one compound in the target compartment this metabolite
            % is at position
            % #totalmets + (#added comparts -1) * #Internal mets + currentmet
            CompTransporters{end}(end+1) = numel(compmodel.rxns);
            Transporters(end+1) = numel(compmodel.rxns);            
        end
        
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
    %Now we adjusted the S matrix, lets remove the surplus stuff.
    compmodel.mets(nonuniquemets(2:end)) = [];
    compmodel.metNames(nonuniquemets(2:end)) = [];
    compmodel.b(nonuniquemets(2:end)) = [];
    compmodel.metFormulas(nonuniquemets(2:end)) = [];        
    compmodel.S(nonuniquemets(2:end),:) = [];
    if (isfield(model,'metChEBIID'))
        compmodel.metChEBIID(nonuniquemets(2:end)) =  [];
    end
    if isfield(model,'metKEGGID')
        compmodel.metKEGGID(nonuniquemets(2:end)) =  [];
    end
    if isfield(model,'metPubChemID')
        compmodel.metPubChemID(nonuniquemets(2:end)) =  [];
    end
    if isfield(model,'metInChIString')
        compmodel.metInChIString(nonuniquemets(2:end)) =  [];
    end
    if isfield(model,'metCharge')
        compmodel.metCharge(nonuniquemets(2:end)) =  [];
    end
    
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
    if isempty(find(ismember(compmodel.rxns,name)))
        %Either add it
        compmodel = addReaction(compmodel,['Ex_' cmet],{cmet},...
            -ExchangeReactions{reac,3},0,0,max(max(compmodel.ub),1000),0);
%        compmodel = addReaction(compmodel,{cmet},...
%            -min(0,ExchangeReactions{reac,3})*min(-1000,min(compmodel.lb)),...
%            max(0,ExchangeReactions{reac,3})*max(max(compmodel.ub),1000));
        %add it to the core.
        CoreReactions(end+1) = numel(compmodel.rxns);
    else
        %Or adjust the bounds.
        pos = find(ismember(compmodel.rxns,name));
        %pos = find(ismember(ComparisonModel.rxns,name));
        direc = full(compmodel.S(ismember(compmodel.mets,cmet),pos));
        CoreReactions(end+1) = pos;
        if sign(ExchangeReactions{reac,3}) ==  sign(direc)
            %we have already added the reaction, but with the opposite
            %directionality
            if compmodel.lb(pos) == 0
                compmodel.rev(pos) = 1;
                compmodel.lb(pos) = min(min(compmodel.lb),-1000);
            end
        end
%         if ExchangeReactions{reac,3} > 0
%             if compmodel.ub(pos) == 0
%                 compmodel.rev(pos) = 1;
%                 compmodel.ub(pos) = max(max(compmodel.ub),1000);
%             end
%         else
%             if compmodel.lb(pos) == 0
%                 compmodel.rev(pos) = 1;
%                 compmodel.lb(pos) = min(min(compmodel.lb),-1000);
%             end
%         end
    end

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


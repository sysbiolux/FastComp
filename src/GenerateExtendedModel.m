function [ExtendedModel,core,transporters,nonLocReacSets,nonLocReacNames] = GenerateExtendedModel(model, ...
                        cytosolID, CompartmentIDs, CompartmentNames, ReactionCompartmentalisationData, GeneCompData, ...
                        SingleCompModel, ExclusiveGenePos, ExternalID,ExchangeReactions)
% Model                                 the original uncompartmentalised Model, 
% epsilon                               the epsilon used for fastcore
% cytosolID                             The ID of the Cytosol
% CompartmentIDs                        a cell aray of compartment identifiers excluding the cytosol and external compartments ('x','m'...)
% CompartmentNames                      a cell array of Compartment Names (in the same order
%                                       as ComprtmentIDs) again excluding the cytosol and external compartments ('Peroxisome','Mitochondrion'....)
% ReactionCompartmentalisationData      a cell array of cell arrays with one  entry
%                                       per reaction in model, assigning reactions to specific compartments.
%                                       {{},{'c','x','p'},{'m'},......}, this includes the demand reactions which
%                                       have to be non specific initially It also includes uptake and export reactions which have to be assigned to their respective compartments.
% GeneCompartmentalisationData          a cell array of cell arrays with one
%                                       entry per gene in the model similar to ReactionCompartmentalisationData
% SingleCompModel                       Whether this is a Model with one single compartment or whether there are more 
%                                       in the original model (i.e. if it has a well defined external compartment
% ExternalMets                          a logical array of the size of model.mets containing all
%                                       external metabolites. This is used to extract all reactions containing
%                                       external metabolites which are to be excluded from duplication.
%                                       The model metabolite ids are not allowed to have [] at any position that
%                                       is not enclosing the compartment id.
% ExclusiveGenePos                      indicates, whether only localised genes will be assigned
%                                       (and only a minimal number of non localised added to fulfill the GPRs) or
%                                       whether all non localised genes will be considered appropriate to all
%                                       compartments.
% ReactionToAllCompartments             indicates whether reactions with evidence are
%                                       only allowed to be present in their assigned compartment, or whether they
%                                       are allowed to be present in multiple compartments - NOT YET IMPLEMENTED
% ExternalID                            ID of the external compartment 
% ExchangeReactions                     A cell array with metabolites and compartment ids and directionality. {'co2' , 'c', -1 ;'co2','e', 1} would 
%                                       indicate that there is an importer for co2 in the cytosol and an exporter in the external compartment
% 
%Adjust the ReacCompartmentalisation Sets
if strcmp(class(ReactionCompartmentalisationData),class(containers.Map()))
    temp = cell(1,numel(model.rxns));
    temp(:) = {{}};
    for k = ReactionCompartmentalisationData.keys
        key = k{1};
        temp{find(ismember(model.rxns,key))} = ReactionCompartmentalisationData(key);        
    end
    ReactionCompartmentalisationData = temp;
end
%And the GeneCompartmentalisation Sets
if strcmp(class(GeneCompData) ,class(containers.Map()))
    temp = cell(1,numel(model.genes));
    temp(:) = {{}};
    for k = GeneCompData.keys
        key = k{1};
        gene = find(ismember(model.genes,key));
        if not(isempty(gene))
            temp{find(ismember(model.genes,key))} = GeneCompData(key);        
        end
    end
    GeneCompData = temp;
end
GeneCompartmentalisationData = GeneCompData;

%Find exchange reactions. There should be none if we also have a
%ExchangeReactions array defined, if, we will simply add those to the
%exchangers for the cytosol in the end.
[exps,imps] =  findImpsAndExpRxns(model);
imps = find(imps);
exps = find(exps);
ExchangeRs = unique([imps; exps]);

%Add all Exchangers to the ExchangeReactions array
%Always add the "upper bound" exchanger first (as this is the one that
%carries positive flux.
for i=1:numel(ExchangeRs)
    met = find(model.S(:,ExchangeRs(i)));
    commet = regexprep(model.mets(met),'\[[a-z]\]$','');
    comp = regexprep(model.mets(met),'.*\[([a-z])\]$','$1');
    if model.ub(ExchangeRs(i)) > 0
        ExchangeReactions(end+1,1:3) = {commet{1},comp{1},full(-sign(model.S(met,ExchangeRs(i))))};
    end
    if model.lb(ExchangeRs(i)) < 0
        ExchangeReactions(end+1,1:3) = {commet{1},comp{1},full(sign(model.S(met,ExchangeRs(i))))};
    end
end

%Now, if we have a non single compartment model we will need to define external metabolites
%if those are not defined by an ID, we assume all metabolites involved with
%Exchange reactions as external. This can include 
if ~SingleCompModel
    if isempty(ExternalID)
        %If no external ID is given, then there are no external
        %metabolites!
        ExternalMets = {};        
    else
        idx =strfind(model.mets,ExternalID);
        ExternalMets = model.mets(find(not(cellfun('isempty',idx))));
    end
else
    %If we have a single comp model, there are no external mets (initially)
    ExternalMets = {};
end

%% now, remove all exchange Reactions!, those reactions will be added later

model = removeRxns(model,model.rxns(union(imps,exps)),0,0); 
%Adjust the ReactionCompData!
ReactionCompartmentalisationData(union(imps,exps)) = [];
%But do not remove the metabolites, as those could otherwise cause
%problems if there is only the exchanger!

%% Get some important pieces of data, e.g. a list of reactions with external metabolites to distinguish reactions to compartmentalise from import/export reactions
% The initial model is uncompartmentalised, with only 2 compartments
% (external, cytosol), or only 1 compartment (singlecompmodel)
% To be able to multiply the model, we will first collect all reactions including external metabolites
% i.e. including reactions which have the "externalID" as ID.
% This is to exclude them from the multiplication
ExtRxns = find(ismember(model.rxns,findRxnsFromMets(model,ExternalMets)));
%fprintf('There are %f external Reactions\n',numel(ExtRxns));
ExtMets = find(ismember(model.mets,ExternalMets));
% size(ExtMets);
%In a singlecompmodel these are ALL Metabolites
IntMets = setdiff(1:numel(model.mets),ExtMets);
% size(IntMets);
%now get all non external reactions  
NonExtRxns = model.rxns;
NonExtRxns(ExtRxns) = [];
NonExtRxns = find(ismember(model.rxns,NonExtRxns));
%The assumption is, that all metabolites contain [c] as an indictaor for
%cytosolic localisation we will have to remove this from the base
%metabolite names
comppattern = '\[[a-zA-Z_0-9]+\]$';
%remove the items enclosed in [] from the metabolite name
NonCompMetList = regexprep(model.mets,comppattern,'');    


%% now prepare the list of assignments for reactions. This includes 2 steps: First the direct reaction associations and then the gene associations
%disp('Assigning Genes');
Complist = ReactionCompartmentalisationData;

% Check each gene, whether it is assigned
% this will contain all reactions which were assigned such that the gene
% responsible for the assignment was essential
%and this is the array of responsible genes
if ~isempty(GeneCompartmentalisationData)
    for g = 1:numel(model.genes)
        %fprintf('Computing gene %s\n',model.genes{g});
        %if there is no assignment ignore the gene, it will later be added to
        %all duplicate reactions in all localisations
        if not(isempty(GeneCompartmentalisationData{g}))
            GeneComparts = GeneCompartmentalisationData{g};
            %fprintf('Computing gene %s\n',model.genes{g});
            %check what effect the deletion of this gene has
            rxns = findRxnsFromGene(model,model.genes{g});
            if not(isempty(rxns))
                %for each reaction check for all valid assignments based on this
                %gene. i.e. double check for invalid assignments due to other
                %genes...
                for ri = 1:numel(rxns)
                    r = rxns(ri);
                    if r > length(model.rxns)
                        disp(model.genes{g})
                        disp(ri)
                    end
                    relgenes = findGenesFromRxns(model,model.rxns(r));
                    relgenes = relgenes{1};
                    FP = FormulaParser();
                    model.grRules{r}
                    Formula = FP.parseFormula(model.grRules{r});
                    %disp('Creating DNNF Form')
                    Formula.convertToDNF();
                    FormString = Formula.toString(0);
                    Formula.reduce();
                    Formula.toString(0);
                    %reduce the Formula.
                    %disp('Reducing Formula')
                    while not(strcmp(FormString,Formula.toString(0)))
                        FormString = Formula.toString(0);
                        Formula.reduce();
                    end
                    % obtain all clauses relevant for the current gene
                    NewFormula = OrNode();
                    if strcmp(class(Formula),class(LiteralNode('a')))
                        NewFormula.addChild(Formula);
                    else
                        
                        
                        for cid = 1:numel(Formula.children)
                            child = Formula.children(cid);
                            if child.contains(model.genes{g})
                                NewFormula.addChild(child);
                            end
                        end
                    end
                    NewFormula.toString(0);
                    %and now get all literals from this new formula and check,
                    %whether there are localisation contradictions
                    locs = getLocalisation(NewFormula,GeneCompartmentalisationData,GeneComparts,model);
                    notaccept = true;
                    while notaccept && isempty(locs)
                        [notaccept,locs,GeneCompartmentalisationData] = CorrectGeneLoc(NewFormula,model,...
                            GeneCompartmentalisationData,g,r);
                    end
                    if not(isempty(locs))
                        % if strcmp(model.genes{g},'YPR140W')
                        %fprintf('Assigning reaction %s to the following compartments due to gene %s',model.rxns{r},model.genes{g});
                        %disp(locs);
                        %end
                        Complist{r} = [Complist{r} locs];
                    end
                end
            end
        end
    end
end

%If we have a non single comp model, we should define, that there are external reactions 
if ~SingleCompModel    
    for r=1:numel(ExtRxns)
        if isempty(Complist{ExtRxns(r)})
            Complist(ExtRxns(r)) = {{ExternalID}};
        end
    end
end

[ExtendedModel,core,transporters,nonLocReacSets,nonLocReacNames] = ...
    MakeCompartmentsOneTransport(model, Complist, IntMets, ExtMets,...
    ExtRxns, NonCompMetList, CompartmentIDs, GeneCompartmentalisationData,...
    ExclusiveGenePos, cytosolID,ExchangeReactions);
%save('Data.mat','ExtendedModel','core','transporters','nonLocReacSets','nonLocReacNames')
end
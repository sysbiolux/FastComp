function [CompartModel,core,transporters,nonLocReacSets] = GenereateExtendedModel(model, epsilon, Demandreactions, ...
                        cytosolID, CompartmentIDs, CompartmentNames, ReactionCompartmentalisationData, GeneCompData, ...
                        SingleCompModel, ExclusiveGenePos, ReactionToAllCompartments,ExternalID)
% Model                                 the original uncompartmentalised Model, 
% Demandreactions                       a map of metabolite names to compartment ids ('GLC' : {'c,'m'}) 
% epsilon                               the epsilon used for fastcore
% CompartmentIDs                        a cell aray of compartment identifiers excluding the cytosol and external compartments ('x','m'...)
% CompartmentNames                      a cell array of Compartment Names (in the same order
%                                       as ComprtmentIDs) again excluding the cytosol and external compartments ('Peroxisome','Mitochondrion'....)
% ReactionCompartmentalisationData      a cell array of cell arrays with one  entry
%                                       per reaction in model, assigning reactions to specific compartments.
%                                       {{},{'c','x','p'},{'m'},......}, this includes the demand reactions which
%                                       have to be non specific initially It also includes uptake and export reactions which have to be assigned to their respective compartments.
% GeneCompartmentalisationData          a cell array of cell arrays with one
%                                       entry per gene in the model similar to ReactionCompartmentalisationData
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
% 
%tic
%disp('Preprocessing')
save('GenerateExtModel.mat');
if strcmp(class(ReactionCompartmentalisationData),class(containers.Map()))
    temp = cell(1,numel(model.rxns));
    temp(:) = {{}};
    for k = ReactionCompartmentalisationData.keys
        key = k{1};
        temp{find(ismember(model.rxns,key))} = ReactionCompartmentalisationData(key);        
    end
    ReactionCompartmentalisationData = temp;
end
%fprintf('There are %f Reactions with localisation information\n', numel(find(not(cellfun(@isempty, ReactionCompartmentalisationData)))));
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

%model.rxns(find(~cellfun(@isempty,regexp(model.rxns,'^DM_'))))
[exps,imps] =  findImpsAndExpRxns(model);
imps = find(imps);
exps = find(exps);
ExchangeRs = unique([imps; exps]);
Internalmets = find(~cellfun(@isempty ,regexp(model.mets,['\[' cytosolID '\]'])));
[~,Internalexchange] = find(model.S(Internalmets,ExchangeRs));
ExchangeRs = setdiff(ExchangeRs,ExchangeRs(Internalexchange));
%fprintf('The following reactions are Exchange Reactions:\n');


%and exclude any internal Demand reactions


if ~SingleCompModel
    if nargin < 12       
        [ExternalMets,reacs] = find(model.S(:,ExchangeRs));
        ExternalMets = model.mets(ExternalMets);
    else
        idx =strfind(model.mets,ExternalID);
        ExternalMets = model.mets(find(not(cellfun('isempty',idx))));
    end
else
    ExternalMets = {};
end
%ExternalMets

if strcmp(class(Demandreactions) ,class({'a'}))
    DReacs = find(ismember(model.rxns,Demandreactions));    
else
    DReacs = Demandreactions;    
end

if size(DReacs,2) > size(DReacs,1)
        DReacs = DReacs';
end        

GeneCompartmentalisationData = GeneCompData;
%fprintf('After preprocessing the model contains %i reactions\n',numel(model.rxns));

%% Get some important pieces of data, e.g. a list of reactions with external metabolites to distinguish reactions to compartmentalise from import/export reactions
% The initial model is uncompartmentalised, with only 2 compartments (external, cytosol).
% To be able to multiply the model, we will first collect all reactions including external metabolites;
% This is to exclude them from the multiplication
ExtRxns = find(ismember(model.rxns,findRxnsFromMets(model,ExternalMets)));
%model.rxns(ExtRxns)
%fprintf('There are %f external Reactions\n',numel(ExtRxns));
ExtMets = find(ismember(model.mets,ExternalMets));
IntMets = setdiff(1:numel(model.mets),ExtMets);
%now get all non external reactions and 
NonExtRxns = model.rxns;
NonExtRxns([ExtRxns ;DReacs]) = [];
NonExtRxns = find(ismember(model.rxns,NonExtRxns));
%The assumption is, that all metabolites contain [c] as an indictaor for
%cytosolic localisation we will have to remove this from the base
%metabolite names
comppattern = '\[[a-zA-Z_0-9]+\]';
NonCompMetList = {};
for m = 1:numel(model.mets)
    %remove the items enclosed in [] from the metabolite name
    NonCompMetList{end+1} = regexprep(model.mets{m},comppattern,'');    
end
%toc
%fprintf('After External determination the model contains %i reactions',numel(model.rxns));

%% now prepare the list of assignments for reactions. This includes 2 steps: First the direct reaction associations and then the gene associations
%disp('Assigning Genes');
Complist = ReactionCompartmentalisationData;

% Check each gene, whether it is assigned
% this will contain all reactions which were assigned such that the gene
% responsible for the assignment was essential
%and this is the array of responsible genes
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
%toc
%fprintf('After Gene Compartment assignment determination the model contains %i reactions\n',numel(model.rxns));
if SingleCompModel    
    for r = 1: numel(ExchangeRs)
        if isempty(Complist{ExchangeRs(r)})
%            fprintf('Assigning exchange reaction %s to the cytosol\n',model.rxns{ExchangeRs(r)});
            Complist(ExchangeRs(r)) = {{cytosolID}};
        end
    end
else
    for r=1:numel(ExtRxns)
        if isempty(Complist{ExtRxns(r)})
%            fprintf('Assigning exchange reaction %s to the External\n',model.rxns{ExchangeRs(r)});
            Complist(ExtRxns(r)) = {{ExternalID}};
        end
    end
    for r = 1: numel(ExchangeRs)
        if isempty(Complist{ExchangeRs(r)})
%            fprintf('Assigning exchange reaction %s to the External\n',model.rxns{ExchangeRs(r)});
            Complist(ExchangeRs(r)) = {{ExternalID}};
        end
    end
end


%% Create the model with all reactions distributed.
% %and now extend the initial model adding a copy of all reactions at each
% %position, each time expanding it by the according assigned reactions,
% %remembering those.
%disp('Creating Compartments')
%save('Complist.mat','Complist');
%fprintf('There are %f localised reactions',numel(find(not(cellfun(@isempty,Complist)))))
[CompartModel,core,transporters,CompTransporters,nonLocReacSets] = MakeCompartmentsOneTransport(model, Complist, IntMets, ExtMets, ExtRxns, NonCompMetList, CompartmentIDs, GeneCompartmentalisationData,ExclusiveGenePos, ReactionToAllCompartments,cytosolID);
%CompartModel.rxns(find(~cellfun(@isempty,regexp(CompartModel.rxns,'^DM_'))))
CompartmentedModel = CompartModel;
%fprintf('After Compartment creation the model contains %i reactions\n',numel(CompartModel.rxns));
%return
%% Define demand reactions for each compartment (i.e. what can be produced in this compartment (and thus does not need to be exported)
if isempty(Demandreactions)
    demandreactions = [];
    ExtendedModel = CompartModel;
else
    [demandreactions, ExtendedModel] = createDemandRxns(Demandreactions,CompartModel);
end
Exchange = find(ismember(ExtendedModel.rxns,findRxnsFromMets(ExtendedModel,ExternalMets)));
core = unique([core]);
coreRxns = ExtendedModel.rxns(core);
transrxns = ExtendedModel.rxns(transporters);
%save('IPExtendedModel.mat','ExtendedModel','core','transporters','model', 'Complist', 'IntMets', 'ExtMets', 'ExtRxns', 'nonLocReacSets', 'CompartmentIDs', 'GeneCompartmentalisationData','ExclusiveGenePos', 'ReactionToAllCompartments','cytosolID');
%toc

%update the compartment specific transporter list 


%[A,s,IP7] = fastcoreIP(core,ExtendedModel,transporters,3,nonLocReacSets);

%compmodel = removeRxns(ExtendedModel,setdiff(ExtendedModel.rxns,ExtendedModel.rxns(A)));

end


function [rxns] = findRxnsFromGene(model,gene)
gpos = find(ismember(model.genes,gene));
rxns = [];
if not(isempty(gpos))
    rxns = find(model.rxnGeneMat(:,gpos)~= 0);    
end
end
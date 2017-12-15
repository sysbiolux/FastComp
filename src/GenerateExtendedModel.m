function [ExtendedModel,core,transporters,nonLocReacSets] = GenerateExtendedModel(model, ...
                        cytosolID, CompartmentIDs, ReactionCompartmentalisationData, GeneCompData, ...
                        SingleCompModel, ExternalID,ExchangeReactions)
% Model                                 the original uncompartmentalised Model, 
% epsilon                               the epsilon used for fastcore
% cytosolID                             The ID of the Cytosol
% CompartmentIDs                        a cell aray of compartment identifiers excluding the cytosol and external compartments ('x','m'...)
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
% ExternalID                            ID of the external compartment 
% ExchangeReactions                     A cell array with metabolites and compartment ids and directionality. {'co2' , 'c', -1 ;'co2','e', 1} would 
%                                       indicate that there is an importer for co2 in the cytosol and an exporter in the external compartment
% 

if isempty(regexp(ExternalID,'^\[.*\]$')) %External ID is assumed to have [];
    ExternalID = ['[' ExternalID ']'];
end


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
LocalExRs = cell(numel(ExchangeRs),5);
for i=1:numel(ExchangeRs)    
    met = find(model.S(:,ExchangeRs(i)));
    commet = regexprep(model.mets(met),'\[[a-z]\]$','');
    comp = regexprep(model.mets(met),'.*\[([a-z])\]$','$1');
    LocalExRs(i,:) = {commet{1},comp{1},full(sign(model.S(met,ExchangeRs(i)))),model.lb(ExchangeRs(i)),model.ub(ExchangeRs(i))};    
end
ExchangeReactions = [ExchangeReactions;LocalExRs];

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
        ExternalMets = model.mets(find(~(cellfun(@isempty,idx))));
    end
else
    %If we have a single comp model, there are no external mets (initially)
    ExternalMets = {};
end

%% now, remove all exchange Reactions!, those reactions will be added later

model = removeRxns(model,model.rxns(union(imps,exps)),'metFlag',false); 
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


Complist = ReactionCompartmentalisationData;

%If we have a non single comp model, we should define, that there are external reactions 
if ~SingleCompModel    
    for r=1:numel(ExtRxns)
        if isempty(Complist{ExtRxns(r)})
            Complist(ExtRxns(r)) = {{regexprep(ExternalID,'^\[(.*)\]$','$1')}};
        end
    end
end

[ExtendedModel,core,transporters,nonLocReacSets] = ...
    MakeCompartmentsOneTransport(model, Complist, IntMets, ExtMets,...
    ExtRxns, NonCompMetList, CompartmentIDs, GeneCompartmentalisationData,...
    cytosolID,ExchangeReactions);
%save('Data.mat','ExtendedModel','core','transporters','nonLocReacSets','nonLocReacNames')
end
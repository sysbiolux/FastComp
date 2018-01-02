function [ExtendedModel,core,transporters,nonLocReacSets] = GenerateExtendedModel(model, ...
                        cytosolID, CompartmentIDs, ReacCompData, GeneCompData, ...
                        SingleCompModel, ExternalID,ExchangeReactions)
% Generate an extended model with unlocalised reactions in all compartments and localised reactions 
% in their specific compartment. 
%
% USAGE:
%    [ExtendedModel,core,transporters,nonLocReacSets] = GenerateExtendedModel(model, ...
%                        cytosolID, CompartmentIDs, ReactionCompartmentalisationData, GeneCompData, ...
%                        SingleCompModel, ExternalID,ExchangeReactions)
%
% INPUTS: 
%    model:                 The compartmentilised model (i.e. the model with all
%                           non localised reactions in all compartments).
%    cytosolID:             The cytosol ID (char)
%    CompartmentIDs:        The compartment ids (cell of strings),
%                           excluding the cytosol and a potential fixed
%                           external compartment.
%    ReacCompData:          Reaction Compartmentalisation data for the
%                           uncompartmentalised model, cell array of cell arrays of strings indicating the compartments.
%    GeneCompData:          As ReacCompData but for Genes (can be empty).
%    SingleCompModel:       Boolean indicator, whether this is a single compartment
%                           model (true -> only cytosol), or has a fixed
%                           external compartment (false).
%    ExternalID:            IDentifier of the external compartment.
%    ExchangeReactions:     List of Exchangers for the model including Demand and uptake reactions.  
%                           This is a n x 5 cell array with the first element being the
%                           metabolite(without compartment) the second element is the compartment ID
%                           The third element is the stoichiometric coefficient 
%                           The fourth and fifth elementes are the lower and upper bound respectively.
%
% OUTPUTS:
%
%    ExtendedModel:         The Model with all transporters, and
%                           compartments.
%    core:                  Positions of the localised reactions in the
%                           Extended Model.                           
%    transporters:          Positions of the transporters in the extended
%                           model.
%    nonLocReacSets:        The sets of non localised reactions. A double
%                           array of indices, with one row per non localised reaction indicating
%                           all positions of the reaction in the different
%                           compartments.%
% .. Authors:
%       - Thomas Pfau 
%                               
                        

if isempty(regexp(ExternalID,'^\[.*\]$')) %External ID is assumed to have [];
    ExternalID = ['[' ExternalID ']'];
end


%Adjust the ReacCompartmentalisation Sets
if strcmp(class(ReacCompData),class(containers.Map()))
    temp = cell(1,numel(model.rxns));
    temp(:) = {{}};
    for k = ReacCompData.keys
        key = k{1};
        temp{find(ismember(model.rxns,key))} = ReacCompData(key);        
    end
    ReacCompData = temp;
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
ReacCompData(union(imps,exps)) = [];
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


Complist = ReacCompData;

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
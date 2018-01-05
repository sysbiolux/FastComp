function [CompModel,Unlocalised, Localisations] = fastComp( model, varargin)
% Network based localisation of reactions using both Gene
% Compartmentalisation data and Reaction compartmentalisation data. Input
% can be provided either as a single compartment model (with all reactions
% in one compartment), or a two compartment model with one external and one
% internal compartment.
%
% USAGE:
%    [CompModel,Unlocalised] = fastComp( model, cytosolID, CompartmentIDs,...
%                                        ReactionCompartmentalisationData, GeneCompartmentalisationData,...
%                                        SingleCompModel,ExternalID,Exchangers)
%
% INPUTS:
%    model:                            the original uncompartmentalised model (or split between cytosol and external). 
%                                      The model must (with any added Exchanger) must be able to carry positive flux in all reactions. 
% OPTIONAL INPUTS:
%    varargin:  variable parameter/Value pairs. The following Parameters
%               are possible:
%               * ReactionLocalisation:  A Cell Array of Cell Arrays containing one or more 
%                                        Strings indicating compartments one entry for all reactions in the model. (e.g. {{'c','m'};{'c'};{'p','x'}})
%                                        Either ReactionCompartmentalisationData or GeneCompartmentalisationData 
%                                        need to be provided.
%               * GeneLocalisation:      A Cell array of cell arrays with localisations for each gene. 
%                                        The following assumption is made wrt activity of a reaction:
%                                        all unlocalised genes are assumed to be present in all compartments.
%                                        all localised genes are assumed to only be present in the indicated compartment.
%                                        only reactions for which the gpr rule evaluates to true under these conditions are assumed to be in the respective compartment. 
%               * Exchangers:            List of Exchangers for the model including Demand and uptake reactions.  
%                                        This is a n x 5 cell array with the first element being the
%                                        metabolite(without compartment) the second element is the compartment ID
%                                        The third element is the stoichiometric coefficient 
%                                        The fourth and fifth elementes are the lower and upper bound respectively.
%               * cytosolID:             the string for the cytosol id (just the letter no enclosing brackets) (default 'c')
%               * CompartmentIDs:        a cell array of letters for the compartments that shall be predicted. Default is computed from geen or reaction localisation.
%               * ExternalID:            ID of the external Compartment (Default 'e')
%               * SingleCompModel:       Indicator whether this is a single compartment model (only cytosol) or if it has external reactions
%                                        Default Computed from the metabolites in the model
% OUTPUTS:
%    CompModel:     A consistent model with reactions localised to
%                   their predicted compartments (i.e. all reactions can carry flux).
%                   Transporters are added as necessary. Unlocalised reactions are added to the cytosol.
%    Unocalised:    A list of reaction IDs for which a compartment could not be predicted.
%    Localisations: A Cell array containing the ID, and a cell array of localisations for
%                   all unlocalised reactions from the input.
%   

parser = inputParser();

parser.addRequired('model',@isstruct);
parser.addParamValue('ReactionLocalisation',repmat({{}},size(model.rxns)),@iscell);
parser.addParamValue('GeneLocalisation',repmat({{}},size(model.genes)),@iscell);
parser.addParamValue('cytosolID','c',@ischar);
parser.addParamValue('ExternalID','e',@ischar);
parser.addParamValue('Exchangers',{},@iscell);
parser.KeepUnmatched = true;


parser.parse(model,varargin{:});
ReactionCompartmentalisationData = parser.Results.ReactionLocalisation;
GeneCompartmentalisationData = parser.Results.GeneLocalisation;
Exchangers = parser.Results.Exchangers;
if all(cellfun(@isempty,ReactionCompartmentalisationData)) && all(cellfun(@isempty,GeneCompartmentalisationData))
    error('Either Gene or Reaction localisation data is required for this function.');
end
cytosolID = parser.Results.cytosolID;

defaultComps = unique(union([ReactionCompartmentalisationData{:}],[GeneCompartmentalisationData{:}]));
if ~isempty(Exchangers)
    defaultComps = union(Exchangers(:,2),defaultComps);
end
defaultComps = setdiff(defaultComps,cytosolID);
metComps = unique(regexprep(model.mets,'.*\[([^\[]+)\]$','$1'));
if ~isempty(setxor(metComps,{cytosolID}))
    defaultSingle = false;
else
    defaultSingle = true;
end
parser.addParamValue('SingleCompModel',defaultSingle,@(x) isnumeric(x) || islogical(x));
parser.addParamValue('CompartmentIDs',defaultComps,@iscell);
parser.parse(model,varargin{:});
SingleCompModel = parser.Results.SingleCompModel;
ExternalID = parser.Results.ExternalID;
CompartmentIDs = parser.Results.CompartmentIDs;

%Make sure, that CompartmentIDs is a row vector
CompartmentIDs = columnVector(CompartmentIDs)';
if ~all(cellfun(@isempty,GeneCompartmentalisationData))
    %Get localisation of each reaction, based on the genes
    %We need to update the localisation. 
    if ~SingleCompModel
        %We need to exclude all external reactions
        idx =regexp(model.mets,['\[' ExternalID '\]$']);
        ExtMets = model.mets(~(cellfun(@isempty,idx)));
        IntReactions = ~ismember(model.rxns,findRxnsFromMets(model,ExtMets));
    else
        IntReactions = true(size(model.rxns));        
    end    
    if ~isfield(model,'rxnGeneMat')
        model = buildRxnGeneMat(model);    
    end
    %Now, for each compartment lets have a look
    unlocalised = cellfun(@isempty, GeneCompartmentalisationData);
    for comp = 1:numel(CompartmentIDs)
        ccomp = CompartmentIDs{comp};
        cloc = cellfun(@(x) any(ismember(x,ccomp)),GeneCompartmentalisationData);
        %Now, the associated reactions are:
        assocreacs = any(model.rxnGeneMat(:,cloc),2);
        actReacs = ismember(model.rxns,findRxnsActiveWithGenes(model,model.genes(cloc | unlocalised)));
        %Now, all reactions active and associated with the genes will get a
        %localisation to this compartment.
        ReactionCompartmentalisationData(actReacs & assocreacs) = cellfun(@(x) union(x,{ccomp}),ReactionCompartmentalisationData(actReacs & assocreacs),'UniformOutput',false);
    end
end

disp('Generating Extended Model');
fastCompModel = model;
fastCompExchangers = Exchangers;
if ~isempty(Exchangers)
    [ComparisonModelFC] = addFastCompExchangersToCytosol(model,Exchangers,cytosolID);
else    
    ComparisonModelFC = fastCompModel;    
end
[CompartModelFC,coreFC,transportersFC,nonLocReacSetsumodFC] = ...
    GenerateExtendedModel(fastCompModel, cytosolID, CompartmentIDs,...
    ReactionCompartmentalisationData, GeneCompartmentalisationData, SingleCompModel,...
    ExternalID, fastCompExchangers);
%Create the mapping.
mappingFC = getOrigReacs(nonLocReacSetsumodFC,CompartModelFC,ComparisonModelFC);
disp('Predicting localisations')
[~,locs] = FastCompartPrediction(CompartModelFC , transportersFC, nonLocReacSetsumodFC, coreFC, [cytosolID,CompartmentIDs], fastCompModel, mappingFC,ReactionCompartmentalisationData, ComparisonModelFC );
%Now, we want all those reactions to be part of the model, and the minimal
%number of transport reactions necessary.
localised = ~cellfun(@isempty, locs);
Unlocalised = model.rxns(mappingFC(~localised));
Localisations = [model.rxns(mappingFC(localised)), locs(localised)];
%Now build the final model.
%We localise all non localised reactions to the cytosol
locs(~localised) = {{'c'}};
CompIDs = [cytosolID CompartmentIDs];
toRemove = false(size(nonLocReacSetsumodFC));
for comp = 1:numel(CompIDs)
    ccomp = CompIDs{comp};
    toRemove(:,comp) = cellfun(@(x) ~any(ismember(x,ccomp)),locs);
end
rxnsToRemove = CompartModelFC.rxns(nonLocReacSetsumodFC(toRemove));
transpReacs = CompartModelFC.rxns(transportersFC);
modelWithLoc = removeRxns(CompartModelFC,rxnsToRemove);
%Do a fastcore run with everything but transporters.
A = fastcore(setdiff(1:numel(modelWithLoc.rxns),find(ismember(modelWithLoc.rxns,transpReacs))),modelWithLoc,1);
CompModel = removeRxns(modelWithLoc,setdiff(modelWithLoc.rxns,modelWithLoc.rxns(A)));    

end

function mapping = getOrigReacs(nonLocReacSets,CompartModel,ComparisonModel)
%get the first non zero entry in each row

[pos,firstNonZero] = max(nonLocReacSets~=0,[],2);
%This gives you a logical array (pos) detailiing the respective row a 1 for
%all entries) and a double index firstNonZero, that indicates the
%respective position.
%So, when taking these entries, it is necessary to get the diagonal matrix.
reacs = CompartModel.rxns(diag(nonLocReacSets(pos,firstNonZero)));
origreacs = regexprep(ComparisonModel.rxns,'\([a-z]\)','');
%This is for exchange reactions.
origreacs = regexprep(origreacs,'\[[a-z]\]','');
newreacs =  regexprep(reacs,'\([a-z]\)','');
%This is for exchange reactions.
newreacs =  regexprep(newreacs,'\[[a-z]\]','');
[Orig,Target] = ismember(newreacs,origreacs);
mapping = Target(Orig);

end

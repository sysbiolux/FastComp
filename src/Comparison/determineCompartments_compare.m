function [FCComps,FCtime,FullComps,Fulltime,MOComps,MOTime] = determineCompartments_Compare( model, epsilon,...
                                                       cytosolID, CompartmentIDs,...
                                                       ReactionLocalisation, GeneLocalisation,...
                                                       SingleCompModel,...
                                                       ExternalID,Exchangers,useMO)
% This function performs a run of fastComp, and the Mintz-ORON MILP
% algorithm. It provides the predicted localisations along with the run
% time for each algorithm.
%
% USAGE:
%    [FCComps,FCtime,FullComps,Fulltime,MOComps,MOTime] = determineCompartments_Compare( model, epsilon,...
%                                                       cytosolID, CompartmentIDs,...
%                                                       ReactionCompartmentalisationData, GeneLocalisation,...
%                                                       SingleCompModel,...
%                                                       ExternalID,Exchangers,useMO)
%
% INPUTS:
%    model:                 the original uncompartmentalised model (or split between cytosol and external). 
%                           The model must (with any added Exchanger) must be able to carry positive flux in all reactions. 
%    cytosolID:             the string for the cytosol id (just the letter no enclosing brackets) 
%    CompartmentIDs:        a cell array of letters for the compartments that shall be predicted. 
%    ReactionLocalisation:  A Cell Array of Cell Arrays containing one or more 
%                           Strings indicating compartments one entry for all reactions in the model. (e.g. {{'c','m'};{'c'};{'p','x'}})
%                           Either ReactionCompartmentalisationData or GeneCompartmentalisationData 
%                           need to be provided.
%    GeneLocalisation:      A Cell array of cell arrays with localisations for each gene. 
%                           The following assumption is made wrt activity of a reaction:
%                           all unlocalised genes are assumed to be present in all compartments.
%                           all localised genes are assumed to only be present in the indicated compartment.
%                           only reactions for which the gpr rule evaluates to true under these conditions are assumed to be in the respective compartment. 
%    SingleCompModel:       Indicator whether this is a single compartment model (only cytosol) or if it has external reactions
%    ExternalID:            ID of the external Compartment 
%    Exchangers:            List of Exchangers for the model including Demand and uptake reactions.  
%                           This is a n x 5 cell array with the first element being the
%                           metabolite(without compartment) the second element is the compartment ID
%                           The third element is the stoichiometric coefficient 
%                           The fourth and fifth elementes are the lower and upper bound respectively.
%    useMO:                 Whether to run the Mintz-Oron algorithm or not.
%    
% OUTPUTS:
%    FCComps:               The Predicted compartments for the consistency
%                           based prediction
%    FCtime:                The runtime for the consistency based prediction
%    FullComps:             The Predicted compartments for fastComp
%    Fulltime:              The runtime for fastcomp
%    MOComps:               The Predicted compartments for the Mintz Oron MILP
%    MOTime:                The runtime of the MIntz-Oron MILP


starttime = clock;
fcpreptime = 0;
disp('Generating Extended Model');
fastCompModel = model;
fastCompExchangers = Exchangers;
if ~isempty(Exchangers)
    [ComparisonModelFC] = addFastCompExchangersToCytosol(model,Exchangers,cytosolID);
else    
    ComparisonModelFC = fastCompModel;    
end

fcpreptime = fcpreptime + etime(clock,starttime)

MOStart = clock;
[CompartModel,core,transporters,nonLocReacSets] = ...
    GenerateExtendedModel(model, cytosolID, CompartmentIDs,...
    ReactionLocalisation, GeneLocalisation, SingleCompModel, ...
    ExternalID, Exchangers);
MOPrepTime = etime(clock,MOStart);
ctime = clock;
[CompartModelFC,coreFC,transportersFC,nonLocReacSetsumodFC] = ...
    GenerateExtendedModel(fastCompModel, cytosolID, CompartmentIDs,...
    ReactionLocalisation, GeneLocalisation, SingleCompModel,...
    ExternalID, fastCompExchangers);
fcmodelgentime = etime(clock,ctime)

%Set Timer
ctime = clock;

FCmodtime = etime(clock,ctime);

mappingFC = getOrigReacs(nonLocReacSetsumodFC,CompartModelFC,ComparisonModelFC);
FCumodtime = etime(clock,ctime) - FCmodtime
ctime = clock;
disp('Doing FastComp')
[FullComps,FCComps,SingleCompTime] = FastCompartPrediction( CompartModelFC , transportersFC, nonLocReacSetsumodFC, coreFC, [cytosolID,CompartmentIDs], fastCompModel, mappingFC,ReactionLocalisation, ComparisonModelFC );
FastCompTime = etime(clock,ctime)
Fulltime =  FastCompTime + FCumodtime + fcpreptime + fcmodelgentime; 
SingleCompTime

FCtime = SingleCompTime + fcmodelgentime + fcpreptime + FCumodtime ;

ctime = clock;
%FullPlusTime =  newelapsed - elapsed;
disp('Doing Mintz-Oron Algorithm')
if useMO
    MOComps = OronMILP(CompartModel,transporters,nonLocReacSets,core,[cytosolID,CompartmentIDs],epsilon, model );
else
    MOComps = FullComps;
end
MOTime = etime(clock,ctime) + MOPrepTime;
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


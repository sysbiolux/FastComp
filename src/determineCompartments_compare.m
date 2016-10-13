function [FCComps,FCtime,FullComps,Fulltime,MOComps,MOTime] = determineCompartments_Compare( model, epsilon, Demandreactions,...
                                                       cytosolID, CompartmentIDs, CompartmentNames,...
                                                       ReactionCompartmentalisationData, GeneCompData,...
                                                       SingleCompModel, ExclusiveGenePos,...
                                                       ExternalID,Exchangers,useMO)
%determineCompartments This function Runs multiple network based
%compartment prediction algorithms on a input set.
% Parameters:
% model                         - the original uncompartmentalised model (or split between cytosol and external) 
% epsilon                       - Epsilon for FastComp, should be 1
% Demandreactions               - no function, will be removed in a final version
% cytosolID                     - the string for the cytosol id (just the letter no enclosing brackets)
% CompartmentIDs                - a cell array of letters for the compartments that shall be predicted.
% CompartmentNames              - Names for the Compartments, (optional)
% ReactionCompartmentalisationData - A Cell Array of Cell Arrays ontaining one or more 
%                                    Strings contained in CompartmentIDs (e.g. {{'c','m'};{'c'};{'p','x'}})

% SingleCompModel               - Indicator whether this is a single compartment model (only cytosol) or if it has external reactions
% ExclusiveGenePos              - Only has an effect 
% ExternalID                    - ID of the external Compartment
% Exchangers                    - List of Exchangers for the model.  This is a Cell Array of 
%                                 Cell Arrays with 3 entries, indicating metabolite, compartment and directionality of the exchanger.
%                                 e.g. {'atp','c',-1; 'glc','e',1}, 
% 
                                                   
starttime = clock;
fcpreptime = 0;
disp('Generating Extended Model');
fastCompmodel = model;
fastCompExchangers = Exchangers;
if ~isempty(Exchangers)
    [fastCompModel,fastCompExchangers,ComparisonModelFC] = CreateComparisonModelAndUpdateExchangersOrder(model,Exchangers,cytosolID);
else
    [~,fastCompModel,~] = fastcc_ChangeModel(model,epsilon);
    ComparisonModelFC = fastCompModel;
    
end
fcpreptime = fcpreptime + etime(clock,starttime);

MOStart = clock;
[CompartModel,core,transporters,nonLocReacSets,nonLocReacNames] = ...
    GenerateExtendedModel(model, cytosolID, CompartmentIDs, CompartmentNames,...
    ReactionCompartmentalisationData, GeneCompData, SingleCompModel, ExclusiveGenePos,...
    ExternalID, Exchangers);
MOPrepTime = etime(clock,MOStart);
ctime = clock;
[CompartModelFC,coreFC,transportersFC,nonLocReacSetsumodFC,nonLocReacNamesumodFC] = ...
    GenerateExtendedModel(fastCompmodel, cytosolID, CompartmentIDs, CompartmentNames,...
    ReactionCompartmentalisationData, GeneCompData, SingleCompModel, ExclusiveGenePos,...
    ExternalID, fastCompExchangers);
fcmodelgentime = etime(clock,ctime);

%Set Timer
ctime = clock;

FCmodtime = etime(clock,ctime);

mappingFC = getOrigReacs(nonLocReacSetsumodFC,CompartModelFC,ComparisonModelFC);
FCumodtime = etime(clock,ctime) - FCmodtime;
ctime = clock;
disp('Doing FastComp')
[FCComps,FullComps,FCtime] = FastCompartPrediction( CompartModelFC , transportersFC, nonLocReacSetsumodFC, coreFC, [cytosolID,CompartmentIDs], epsilon, fastCompModel, mappingFC,ReactionCompartmentalisationData, ComparisonModelFC ,'umod');
Fulltime = etime(clock,ctime) + FCumodtime + fcpreptime + fcmodelgentime; 

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


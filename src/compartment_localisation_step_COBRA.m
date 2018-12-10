function [localisedReactions] = compartment_localisation_step_COBRA( model,CompIDs, localisedReactions, noncompmodel,...
    nonLocSets,ReacCompData, epsilon,...
    TransportedMetabs, TransPresenceInComp,TransDirections)
% Assign localisations based on model consistency, i.e. use all reactions
% localised to a compartment along with a set of exchange reactions for the
% compartments to indicate additional reactions in the compartment.
% This function exclusively uses COBRA Toolbox functions, while the
% original relies on the original implementation of fastcore.
%
% USAGE:
%    [localisedReactions] = compartment_localisation_step_COBRA( model,CompIDs, localisedReactions, noncompmodel,...
%    nonLocSets,ReacCompData, epsilon,...
%    TransportedMetabs, TransPresenceInComp,TransDirections)
%
% INPUTS: 
%    model:                 The compartmentilised model (i.e. the model with all
%                           non localised reactions in all compartments).
%    CompIDs:               The IDs of the compartment (a Cell Array of strings).
%    localisedReactions:    Indices of all reactions which are localised in the compartmentalised model.
%    noncompmodel:          The model without compartmentalisation (e.g. only
%                           the cytosol with exchangers for 
%    nonLocSets:            The sets of non localised reactions. A double
%                           array of indices, with one row per non localised reaction indicating
%                           all positions of the reaction in the different
%                           compartments.
%    ReacCompData:          Reaction Compartmentalisation data for the
%                           uncompartmentalised model.
%    epsilon:               activity epsilon.
%    TransportedMetabs:     Indicators which MEtabolites are transported.
%    TransPresenceInComp:   Indicator, which metabolites are transported in
%                           which compartment.
%    TransDirections:       Indication of the directionality of transport
%                           to/from a compartment.
%
% OUTPUTS:
%
%    localisedReactions:    Indices of all reactions localised in the
%                           compartmentalised model. I.e. all indices not
%                           returned for reactions with at least on
%                           localised reaction can be deactivated.
%
%
% .. Authors:
%       - Maria Pires Pacheco - initial implementation
%       - Thomas Pfau - addition of Transport restrictions
%

comp_model=model;
mother_model=noncompmodel;

%first determine the reactions, which are localised but have no Compartment
%Assignment. These reactions are the external reactions (and should be
%assigned to the cytosol.

%Set up the transporters
%To do so, we have to determine which metabolites are actually exchanged.
ActiveTransport = sum(TransPresenceInComp,2) > 0;
ActiveTransMetabs = TransportedMetabs(ActiveTransport);
ActiveDirections = TransDirections(ActiveTransport,:);
ActivePresence = TransPresenceInComp(ActiveTransport,:);

%Now, we don't care about duplication. for this test.
[TransMetabPresence,TransMetabPositions] = ismember(mother_model.mets,strcat(ActiveTransMetabs,['[' CompIDs{1} ']']));
ActiveDirections = ActiveDirections(TransMetabPositions(TransMetabPresence),:);
ActivePresence = ActivePresence(TransMetabPositions(TransMetabPresence),:);
transporterpositions = zeros(numel(ActiveTransMetabs),1);
for i = 1:numel(ActiveTransMetabs)
    
    [mother_model,reacID] = addReaction(mother_model,['FC_Exchange_' ActiveTransMetabs{i}],{[ActiveTransMetabs{i} '[' CompIDs{1} ']']},[-1],...
        1);
    mother_model.lb(reacID) = 0;
    mother_model.ub(reacID) = 0;
    if(isempty(reacID))
        reacID = numel(mother_model.rxns);
    end
    transporterpositions(i) = reacID;
end

%Obtain a list of localised reactions
Localised = find(~cellfun(@isempty,ReacCompData));
%Initialise the localisation information
for j=1:numel(Localised)
    presence_tag(j,:)=ismember(CompIDs, ReacCompData{Localised(j)});
end
%and attach the Transporter Localisation
presence_tag = [presence_tag;ActivePresence];
Localised = [Localised; transporterpositions];

%Initialise A_keep as a #Compartments * #noncompreactions cell array
A_keep=cell(numel(CompIDs), numel(mother_model.rxns));
A_sets=cell(numel(CompIDs), 1);
%Exchangers are set up as "Ex_metaboliteID[compID]\d+"
Exchangers = ~cellfun(@isempty, regexp(mother_model.rxns,'^Ex_.*\[[^\]]+\][0-9]+$'));
ExComps = regexprep(mother_model.rxns,'^Ex_.*\[([^\]]+)\][0-9]+$','$1');

for i=1:numel(CompIDs)
    model=mother_model;
    %Define the reactions that are to be removed from the model. (for wrong
    %localisation
    to_remove=Localised( presence_tag(:,i)==0);
    %And define the core as the reactions localised to this compartment.
    C=setdiff(Localised( presence_tag(:,i)==1),transporterpositions);
    model.ub(transporterpositions) = ((ActivePresence(:,i) == 1) & (ActiveDirections(:,i) == 1)) * max(model.ub);
    model.lb(transporterpositions) = ((ActivePresence(:,i) == 1) & (ActiveDirections(:,i) == -1)) * min(model.lb);
    %Further shut down all Exchangers which are not in this compartment.    
    CompEx = cellfun(@(x) isequal(x,CompIDs{i}),ExComps);
    %Shut down all exchangers which are not in this compartment.
    model.lb(Exchangers & ~CompEx) = 0;
    model.ub(Exchangers & ~CompEx) = 0;    
    %reduce the model (i.e. remove all "localised" reactions with the wrong localisation).
    model=removeRxns(model,model.rxns(to_remove));
    %obtain the consistent model (for numerical reasons, we repeat this
    %process until it converges. normally, one run is sufficient, but there
    %are situations where the numerics of the system allow a flux over the
    %threshold by distributing it over many reactions and make it
    %infeasible in the next step.
    reacset = 1:numel(model.rxns);
    A=fastcc(model, epsilon,0);
    while numel(reacset) > numel(A)
        %Reduce the model to only contain reactions with this
        model=removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A)));
        reacset = 1:numel(model.rxns);
        A = fastcc(model,epsilon,0);
    end
    %Define the new core as all localised reactions in the original model and the
    %reactions of the current model.
    %i.e. remove anything from the core, that is not present in the current
    %reduced model (i.e. that cannot be active in the current submodel).
    C=find(ismember(model.rxns, mother_model.rxns(C)));
    
    %Determine the fastcore of the new model.
    newmodel = fastcore(model,C,epsilon,0 );
    A = find(ismember(model.rxns,newmodel.rxns));
    if numel(A) > 0
        A_keep(i,1:numel(A))=model.rxns(A)';
        A_sets{i} = model.rxns(A)';
    end
end

A_keep2=A_keep(:);
not_empty =find(~cellfun('isempty',A_keep2));
[A_keep2_uni, IA]=unique(A_keep2(not_empty));
%These are all detected reactions. NOw, see in which compartments they are
%present
found_tags = zeros(size(A_keep2_uni'));
if ~isempty(A_sets{1})
    found_tags(1,:) = ismember(A_keep2_uni,A_sets{1});
else
    found_tags(1,:) = 0;
end
A_keep2_uni_name = {};
for i = 2:numel(CompIDs)
    if ~isempty(A_sets{i})
        presence = ismember(A_keep2_uni,A_sets{i})';
    else
        presence = ismember(A_keep2_uni,{})';
    end
    found_tags(i,:) = presence;
    found_tags(1,:) = found_tags(1,:) & ~presence;
    A_keep2_uni_name = [A_keep2_uni_name; strcat(A_keep2_uni(presence),'(',CompIDs(i),')')];
end
A_keep2_uni_name = [A_keep2_uni_name; strcat(A_keep2_uni(logical(found_tags(1,:))),'(',CompIDs(1),')')];



new_localized=find(ismember(comp_model.rxns,A_keep2_uni_name ));
localisedReactions=unique([localisedReactions; new_localized]);

end




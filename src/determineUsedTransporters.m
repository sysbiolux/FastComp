function [Metabolites,AvailableInComp,Direction] = determineUsedTransporters(model , weightedRxns, nonLocSets,...
                                   localisedReactions, CompIDs, epsilon,...
                                   noncompmodel)
% Determine the transporters used in each compartment when activating all
% localised reactions.
%
% USAGE:
%    [Metabolites,AvailableInComp,Direction] = determineUsedTransporters(model , weightedRxns, nonLocSets,...
%                                   localisedReactions, CompIDs, epsilon,...
%                                   noncompmodel)
%
% INPUTS: 
%    model:                 The compartmentilised model (i.e. the model with all
%                           non localised reactions in all compartments).
%    weightedRxns:          The indices of all transporters which should receive a penalty.
%    nonLocSets:            The sets of non localised reactions. A double
%                           array of indices, with one row per non localised reaction indicating
%                           all positions of the reaction in the different
%                           compartments.
%    localisedReactions:    Indices of all reactions which are localised in the compartmentalised model.
%    CompIDs:               The compartment ids (cell of strings)
%    epsilon:               activity epsilon.
%    noncompmodel:          A non com,partmentalised model.
%
% OUTPUTS:
%
%    Metabolites:           A list of transported metabolites.
%    AvailableInComp:       A indicator matrix which metabolites are
%                           transported into/from which compartment.
%    Direction:             The direction of the transport, per
%                           compartment.
%
% .. Authors:
%       - Thomas Pfau 
%                               
                              
t = numel(weightedRxns);
 
transpenalty = ones(t,1);
%get the sums of involved reactions for each original metabolite
%(otherwise this is too dependent on the number of compartments)
% these 
MetOcc = noncompmodel.S ~= 0;
MetSums = sum(MetOcc,2);
%Find the most common metabolites to exclude them from penalisation
relmets = find(MetSums > 0.05 * numel(noncompmodel.mets));
%Get the positions of these Metabolites in the Compartmentalised model
CompModMets = find(ismember(model.mets,noncompmodel.mets(relmets)));
%Get the positions of the reactions these metabolites are involved with
[~,InvolvedReacs] = find(model.S(CompModMets,:));
InvolvedReacs = unique(InvolvedReacs);
% These transports will not be punished.
InvolvedTrans = intersect(weightedRxns,InvolvedReacs);
transpenalty(ismember(weightedRxns,InvolvedReacs)) = 0;

%if everything is done, return the solution.
if numel(nonLocSets) == 0
    return
end

lp = Cplex('LP');
m = numel(model.mets);
n = numel(model.rxns) ;
k = numel(nonLocSets(nonLocSets~=0));
q = numel(localisedReactions);
nonLocReacsVector = reshape(nonLocSets',1,[]);
nonLocReacs = nonLocReacsVector(nonLocReacsVector~=0);
NP(1:k) = nonLocReacs(1:k);
NL = zeros(n,1);  NL(NP) = 1;   NL = diag(NL);  NL(find(all(NL==0,2)),:) = [];
L = zeros(n,1); L(localisedReactions) = 1; L = diag(L); L(find(all(L==0,2)),:) = [];
TP = zeros(n,1); TP(weightedRxns) = 1; TP = diag(TP); TP(find(all(TP==0,2)),:) = [];

%Create the A matrix
 A =   [model.S  , sparse(m,q) , sparse(m,t) ; % S*v = 0
  L , -speye(q,q), sparse(q,t) %  z <= v for all Localised Reactions
  TP, sparse(t,q), -speye(t,t) ; % v - z <= 0
  TP, sparse(t,q), speye(t,t) ]; % 0 <= v + z

rhs = [zeros(m,1) ; % S*v = 0;       
       inf * ones(q,1);% lb <= (lb - eps)*y+ + v for all 'Localised Reactions'                      
       zeros(t,1) ; % v - z <= 0
       inf * ones(t,1) ]; % 0 <= v + z

lhs =  [zeros(m,1) ; % S*v = 0;
       zeros(q,1) ; % z <= v for all Localised Reactions
       -inf * ones(t,1) ; % v - z <= 0
       zeros(t,1) ]; % 0 <= v + z


lbs = [model.lb;zeros(q,1);zeros(t,1)]; % lb <= v; 0 <= z ; 0 <= t
ubs = [model.ub;epsilon*ones(q,1);inf * ones(t,1)]; % v <= ub; z <= 1 ; t <= inf

%set names for debugging (this could be removed for performance)
colnames = [model.rxns ;
    strcat(model.rxns(localisedReactions),'_u');
    strcat(model.rxns(weightedRxns),'_z')] ;
    
rownames = [model.mets;
    strcat(model.rxns(localisedReactions),'_vgtu');
    strcat(model.rxns(weightedRxns),'_zpos') ; 
    strcat(model.rxns(weightedRxns),'_zneg') ];

%The Objective for LP7 is to dependent on the current "core" (even though
%we do not have to distingish between forward and backward)
%Initially however it is simply all a's and all localised zs

%The objective is Sum z_i - eps * 1/(#transports) * Sum (flux through transports)
obj = [zeros(n,1);ones(q,1);epsilon*-1/t*transpenalty];

%We only have continous variables
coltype = '';
coltype(1,1:n) = 'C';
coltype(1,(end+1):(end+q+t)) = 'C';

%Set up the cplex model.
lp.Model.A = A;
%turn of the Display
lp.Param.output.clonelog.Cur = -1;
lp.Param.simplex.display.Cur = 0;
lp.Param.barrier.display.Cur = 0;
lp.Param.mip.display.Cur = 0;
lp.Model.lb = lbs;
lp.Model.ub = ubs;
lp.Model.rhs = rhs;
lp.Model.lhs = lhs;
lp.Model.colname = colnames;
lp.Model.rowname = rownames;
lp.Model.sense = 'maximize';
lp.Model.obj = obj;
%lp.Model.ctype = coltype;

sol = lp.solve();
reactionVals = sol.x(1:n);
activetranspos = reactionVals(weightedRxns) > 1e-6;
activetransneg = reactionVals(weightedRxns) < 1e-6;

Transporterspos = model.rxns(weightedRxns(activetranspos));
Transportersneg = model.rxns(weightedRxns(activetransneg));
CompIDRegexp = '';
for i = 1:numel(CompIDs)
    CompIDRegexp = [CompIDRegexp ,'(' , CompIDs{i}  ')'];
end

Metabolites = unique(regexprep(model.mets,['\[[' CompIDRegexp ']\]$'],''));

Compspos = cellfun(@(x) x(4), Transporterspos,'UniformOutput',false);
Compsneg = cellfun(@(x) x(4), Transportersneg,'UniformOutput',false);
transmetabolitespos = cellfun(@(x) x(7:end), Transporterspos,'UniformOutput',0);
transmetabolitesneg = cellfun(@(x) x(7:end), Transportersneg,'UniformOutput',0);

AvailableInComp = zeros(numel(Metabolites),numel(CompIDs));
Direction = zeros(numel(Metabolites),numel(CompIDs));
for i = 1:numel(transmetabolitespos)
    metpos = ismember(Metabolites,transmetabolitespos{i});
    comppos = ismember(upper(CompIDs),Compspos{i});    
    AvailableInComp(metpos,comppos) = 1;
    Direction(metpos,comppos) = 1;    
    %And make it available in the cytosol!!
    %Don't make it available in the cytosol, to avoid overrepresentation.
    %AvailableInComp(metpos,1) = 1;    
    %Direction(metpos,1) = 1;
end

for i = 1:numel(transmetabolitesneg)
    metpos = ismember(Metabolites,transmetabolitesneg{i});
    comppos = ismember(upper(CompIDs),Compsneg{i});    
    AvailableInComp(metpos,comppos) = 1;
    Direction(metpos,comppos) = -1;    
    %And make it available in the cytosol!!
    %Don't make it available in the cytosol, to avoid overrepresentation.
    %AvailableInComp(metpos,1) = 1;
    %Direction(metpos,1) = -1;    
end


    
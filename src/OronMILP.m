function comps = OronMILP( model , weightedRxns, nonLocSets , localisedReactions ,CompIDs,epsilon,noncompmodel)
% This function assumes the following:
% Irreversible reactions have a lower bound of 0. 
% Otherwise bounds are +/- 10^9 for upper/lower bounds respectively
% The structure of the IP will be as follows:
% With 
% m = numel(model.mets);
% n = numel(model.rxns); 
% k = numel(nonLocSets);
% q = numel(localisedReactions)
% t = numel(weightedRxns)
% NP = unique(nonLocSets);
% NL = zeros(n,1); NL(NP) = 1; NL = diag(NL); NL(find(all(NL==0,2))) = [];
% L = zeros(n,1); L(localisedReaction) = 1; L = diag(L); L(find(all(L==0,2))) = [];
% TP = zeros(n,1); TP(weightedRxns) = 1; TP = diag(TP); TP(find(all(TP==0,2))) = [];
% 
% A = 
%  [model.S  , sparse(m,2*q) , sparse(m,2*k) , sparse(m,t) ; % S*v = 0
%  speye(q,n) , sparse(q,q), diag(model.ub(localisedReactions)+eps), sparse(q,t) %  v + y- * (ub + eps)  <= ub for all 'Localised Reactions'
%  speye(q,n) , diag(model.lb(localisedReactions)-eps) , sparse(q,q) , sparse(q,t) ; % lb <= (lb - eps)*y+ + v for all 'Localised Reactions'
%  speye(k,n) , sparse(k,k), diag(model.ub(NP)+eps), sparse(k,t) %  v + y- * (ub + eps)  <= ub for all 'non Localised Reactions'
%  speye(k,n) , diag(model.lb(NP)-eps) , sparse(k,k) , sparse(k,t) ; % lb <= (lb - eps)*y+ + v for all 'non Localised Reactions'
%  speye(k,n) , NL , NL , sparse(k,t) ; % y+ + y- = ? for all 'non Localised Reactions'
%  TP, sparse(k,q) , sparse(k,q)  , speye(k,k) ; % v - z <= 0
%  TP, sparse(k,q) , sparse(k,q)  , speye(k,k) ; % 0 <= v + z
%]
% Order of variables:
% R1 .. RN , y+1...y+(q+k) , y-1...y-(q+k), z1 .. zk 
% 
%save('OronMILP.mat')
ip = Cplex('MILP');
 

%%
m = numel(model.mets);
n = numel(model.rxns) ;
k = numel(nonLocSets);
q = numel(localisedReactions);
t = numel(weightedRxns);
NP(1:k) = nonLocSets(1:k);
NL = zeros(n,1); NL(NP) = 1; NL = diag(NL); NL(find(all(NL==0,2)),:) = [];
L = zeros(n,1); L(localisedReactions) = 1; L = diag(L); L(find(all(L==0,2)),:) = [];
TP = zeros(n,1); TP(weightedRxns) = 1; TP = diag(TP); TP(find(all(TP==0,2)),:) = [];
%%
%Create the A matrix
 A=    [model.S  , sparse(m,2*q) , sparse(m,2*k) , sparse(m,t) ; % S*v = 0
  L , sparse(q,q), diag(model.ub(localisedReactions)+epsilon), sparse(q,2*k), sparse(q,t); %  v + y- * (ub + eps)  <= ub for all 'Localised Reactions'
  L , diag(model.lb(localisedReactions)-epsilon) , sparse(q,q), sparse(q,2*k) , sparse(q,t) ; % lb <= (lb - eps)*y+ + v for all 'Localised Reactions'
  NL,  sparse(k,2*q), sparse(k,k), diag(model.ub(NP)+epsilon), sparse(k,t); %  v + y- * (ub + eps)  <= ub for all 'non Localised Reactions'
  NL ,  sparse(k,2*q), diag(model.lb(NP)-epsilon) , sparse(k,k) , sparse(k,t) ; % lb <= (lb - eps)*y+ + v for all 'non Localised Reactions'
  sparse(k,n) ,  sparse(k,2*q), speye(k) , speye(k) , sparse(k,t) ; % y+ + y- = ? for all 'non Localised Reactions'
  TP, sparse(t,2*q) , sparse(t,2*k)  , -speye(t,t) ; % v - z <= 0
  TP, sparse(t,2*q) , sparse(t,2*k)  , speye(t,t) ]; % 0 <= v + z

%initialize rhs
rhs = [zeros(m,1) ; % S*v = 0;
       model.ub(localisedReactions) ; %  v + y- * (ub + eps)  <= ub for all 'Localised Reactions'
       inf * ones(q,1);% lb <= (lb - eps)*y+ + v for all 'Localised Reactions' 
       model.ub(NP); %  v + y- * (ub + eps)  <= ub for all 'non Localised Reactions'
       inf * ones(k,1); %  % lb <= (lb - eps)*y+ + v for all 'non Localised Reactions'
       ones(k,1); %  y+ + y- = ? for all 'non Localised Reactions' initially
       zeros(t,1) ; % v - z <= 0
       inf * ones(t,1) ]; % 0 <= v + z

lhs =  [zeros(m,1) ; % S*v = 0;
       -inf * ones(q,1) ; % v + y- * (ub + eps)  <= ub for all 'Localised Reactions'
       model.lb(localisedReactions); % lb <= (lb - eps)*y+ + v for all 'Localised Reactions'
       -inf * ones(k,1) %  v + y- * (ub + eps)  <= ub for all 'non Localised Reactions'
       model.lb(NP) ; % lb <= (lb - eps)*y+ + v for all 'non Localised Reactions'
       zeros(k,1) ; % y+ + y- = ? for all 'non Localised Reactions' initially
       -inf * ones(t,1) ; % v - z <= 0
       zeros(t,1) ]; % 0 <= v + z


lbs = [model.lb;zeros(q,1);zeros(q,1); zeros(2 * k + t,1)]; % essentially all the other variables are strictly positive
ubs = [model.ub;ones(2*q + 2*k,1);inf * ones(t,1)]; % except for the u's all other helper variables are binary
%set names for debugging (this will be removed later on for performance)
colnames = [model.rxns ;
    strcat(model.rxns(localisedReactions),'_ypos');
    strcat(model.rxns(localisedReactions),'_yneg');    
    strcat(model.rxns(NP),'_ypos');
    strcat(model.rxns(NP),'_yneg');
    strcat(model.rxns(weightedRxns),'_z')] ;
    
rownames = [model.mets;
    strcat(model.rxns(localisedReactions),'_ypos');
    strcat(model.rxns(localisedReactions),'_yneg');    
    strcat(model.rxns(NP),'_ypos');
    strcat(model.rxns(NP),'_yneg');    
    strcat(model.rxns(NP),'_nonzero');
    strcat(model.rxns(weightedRxns),'_zpos') ; 
    strcat(model.rxns(weightedRxns),'_zneg') ];

%Determine transporters of reactions with multiple entries (>20);
transpenalty = ones(t,1);
%get the sums of involved reactions for each original metabolite
MetOcc = noncompmodel.S ~= 0;
MetSums = sum(MetOcc,2);
%Find the most common metabolites to exclude them from penalisation
relmets = find(MetSums > 20);
%Get the positions of these Metabolites in the Compartmentalised model
CompModMets = find(ismember(model.mets,noncompmodel.mets(relmets)));
%Get the positions of the reactions these metabolites are involved with
[~,InvolvedReacs] = find(model.S(CompModMets,:));
InvolvedReacs = unique(InvolvedReacs);

%printRxnFormula(model,model.rxns(weightedRxns(ismember(weightedRxns,InvolvedReacs))));
transpenalty(ismember(weightedRxns,InvolvedReacs)) = 0;
    
obj = [zeros(n,1);ones(2*q,1);zeros(2*k,1);-1/(100*t)*transpenalty];

%Set up the MIP
coltype = '';
coltype(1,1:n) = 'C';
coltype(1,n+1:n+2*q+2*k) = 'I';
coltype(1,end+1:end+t) = 'C';

ip.Model.A = A;
ip.DisplayFunc = [];
ip.Model.lb = lbs;
ip.Model.ub = ubs;
ip.Model.rhs = rhs;
ip.Model.lhs = lhs;
ip.Model.colname = colnames;
ip.Model.rowname = rownames;
ip.Param.output.clonelog.Cur = -1;
ip.Model.sense = 'maximize';
ip.Model.obj = obj;
ip.Model.ctype = coltype;
%Set the integrality tolerance to zero! If we have an integer it has to be
%integer!
ip.Param.mip.tolerances.integrality.Cur = 0;
%use the highest possible precision
ip.Param.simplex.tolerances.feasibility.Cur = 1e-9;
sol = ip.solve();
%Set a time limit
%ip.Param.timelimit.Cur = 7200;
scores = calcScores(model,ip,NP,m+2*q+2*k,nonLocSets);
scores = vec2mat(scores,size(nonLocSets,1))';
comps = GetCompartmentsFromScore(scores, CompIDs, ip.Param.simplex.tolerances.feasibility.Cur);
save('MO.mat');

end

function [MILP,score] = calcScore(MILP,pos,altpos)
MILP.Param.output.clonelog.Cur = -1;
MILP.Param.mip.tolerances.integrality.Cur = 0;
MILP.Param.mip.tolerances.absmipgap.Cur = 0;
MILP.Param.mip.tolerances.mipgap.Cur = 0;
MILP.Param.read.scale.Cur = -1;

    try    
        %force flux through the reaction by setting the respective activity
        %constraint to and inactivating the alternatepositions.
        MILP.Model.lhs(altpos) = 0;
        MILP.Model.rhs(altpos) = 0;        
        MILP.Model.lhs(pos) = 1;
        MILP.Model.rhs(pos) = 1;
        sol = MILP.solve();
        score = sol.objval;        
        MILP.Model.lhs(altpos) = 0;
        MILP.Model.rhs(altpos) = 1;        
        MILP.Model.lhs(pos) = 0;
        MILP.Model.rhs(pos) = 1;        
    catch        
        %if this fails, we could not calculatee a score, so the score
        %becomes -inf
        score = -inf;
        MILP.Model.lhs(pos) = 0;
        MILP.Model.rhs(pos) = 1;
        %Deactivate the activity constraint again.
        MILP.Model.lhs(altpos) = 0;
        MILP.Model.rhs(altpos) = 1; 
    end
    
end

function alternatepos = getAltPositions(pos,nonLocSets)
entry = nonLocSets(pos);
[row,~] = find(nonLocSets == entry);
alts = nonLocSets(row,:);
alternatepos = find(ismember(nonLocSets(:),alts));

end

function scores = calcScores(model,MILP,nonLocReacs, offset, nonLocSets)
scores = -inf * ones(size(nonLocReacs));
MILP.Param.output.clonelog.Cur = -1;
starttime = toc;
for i=1:numel(nonLocReacs)
    altpos = getAltPositions(i,nonLocSets);
    %if we take longer than 4 days we stop the computation.
    if toc - starttime > 345600
        break
        disp('Mintz Oron reached time limit')
    end
    %calculate the score for reaction i
    [MILP,score] = calcScore(MILP,offset+i,offset+altpos);
    scores(i) = score;
    if scores(i) == -inf
        fprintf('ORONMIL (using P: No Solution found for reactions', model.rxns{nonLocReacs(i)});
    end
    %When the score is calculated, set the others as possible again.
   % MILP.Model.rhs([altpos;i]) = 1;
end
end


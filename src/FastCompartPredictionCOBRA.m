function [FinalComps,FCComps,FCTime] = FastCompartPrediction( model , weightedRxns, nonLocSets,...
                                                              localisedReactions, CompIDs, ...
                                                              noncompmodel, nonLocToOrigMap, ReacLocal,...
                                                              ComparisonModel)
% Network based localisation of reactions using both Gene
% Compartmentalisation data and Reaction compartmentalisation data. Input
% can be provided either as a single compartment model (with all reactions
% in one compartment), or a two compartment model with one external and one
% internal compartment.
%
% USAGE:
%    [FCComps,FinalComps,FCTime] = FastCompartPrediction( model , weightedRxns, nonLocSets,...
%                                                              localisedReactions, CompIDs, ...
%                                                              noncompmodel, nonLocToOrigMap, ReactionCompartmentalisationData,...
%                                                              ComparisonModel)%
% INPUTS:
%    model:                 The Extended Model in which all non Localised reactions are 
%                           present in all compartments, as generated by generateExtendedModel  
%    weightedRxns:          Indices of all transport or similar reactions which
%                           should be weighted during fastcore.
%    nonLocSets:            A double matrix with one row per reaction containing the 
%                           indices of that reaction in each compartment. each
%                           column corresponds to one compartment from the CompIDs
%                           cell array.
%    localisedReactions:    A Set of positions of reactions which are localised.
%    CompIDs:               A cell array of letters indicating compartment ids arrays
%    noncompmodel:          The uncompartmentalised model to obtain involvement of metabolites
%    nonLocToOrigMap:       A Cell array mapping the non localised reactions to the original reaction names in noncompmodel
%    ReacLocal:             The Compartmentalisation data, see the format
%                           description in fastComp
%    ComparisonModel:       a model for comparison in the preprocessing
%                           step (Essentially the original model with added
%                           transporters).
%
% OUTPUTS:
%    FinalComps:            A Cell array indicating localisation for each
%                           reaction from the nonLocSets.
%    FCComps:               A Cell Array for all reactions which were
%                           directly assigned by the first consistency
%                           based assignment run.
%    FCTime:                The Time it took to run the initial consistency
%                           based assignment step.
%
% NOTE:
%     FastComp will create a linear problem and solve it to obtain scores.
%     The structure of the LP will be as follows:
%     With 
%     m = numel(model.mets);
%     n = numel(model.rxns); 
%     k = numel(nonLocSets);
%     q = numel(localisedReactions)
%     t = numel(weightedRxns)
%     NP = unique(nonLocSets);
%     L = zeros(n,1); L(localisedReaction) = 1; L = diag(L); L(find(all(L==0,2))) = [];
%     TP = zeros(n,1); TP(weightedRxns) = 1; TP = diag(TP); TP(find(all(TP==0,2))) = [];
% 
%     A = 
%     [model.S  , sparse(m,q) , sparse(m,k) , sparse(m,t) ; % S*v = 0
%     L , speye(q,q), sparse(q,k), sparse(q,t) %  z <= v for all Localised Reactions
%     TP, sparse(t,q), -speye(t,t) ; % v - z <= 0
%     TP, sparse(t,q), speye(t,t)  , speye(k,k) ; % 0 <= v + z
%     ]
%     Order of variables:
%     R1 .. RN , u1..uq, z1 .. zk 
% 
% Author: Thomas Pfau 2017

                                                          

epsilon = 1;
ctime = clock;
orignonLocSets = nonLocSets;

% Perform one run of the FastCompLP to determine the used transporters
[Metabolites,AvailableInComp,TransDirections] = determineUsedTransportersCOBRA(model , weightedRxns, nonLocSets,...
localisedReactions, CompIDs, epsilon,...
noncompmodel);

% Perform an initial pre_localisation step to get an initial set of
% assigned reactions.
[new_localised] = compartment_localisation_step_COBRA( model,CompIDs, localisedReactions,...
                                         ComparisonModel, nonLocSets,ReacLocal,...
                                         epsilon, Metabolites, AvailableInComp,TransDirections);


t = numel(weightedRxns);
 
transpenalty = ones(t,1);
%get the sums of involved reactions for each original metabolite
%(otherwise this is too dependent on the number of compartments)
% these 
MetOcc = noncompmodel.S ~= 0;
MetSums = sum(MetOcc,2);
%Find the most common metabolites to exclude them from penalisation
relmets = find((MetSums > 0.05 * numel(noncompmodel.rxns)) & (MetSums > 10));
%Get the positions of these Metabolites in the Compartmentalised model
CompModMets = find(ismember(model.mets,noncompmodel.mets(relmets)));
%Get the positions of the reactions these metabolites are involved with
[~,InvolvedReacs] = find(model.S(CompModMets,:));
InvolvedReacs = unique(InvolvedReacs);
% These transports will not be punished.
InvolvedTrans = intersect(weightedRxns,InvolvedReacs);
transpenalty(ismember(weightedRxns,InvolvedReacs)) = 0;

%Extract all Reactions from nonLocSets which are found to be localised in the preprocessing step 
found = intersect(new_localised,nonLocSets);

%remove all nonLocSets which have entries in found
[rows,cols] = find(ismember(nonLocSets,found));
fixedrows= unique(rows);
predrows = 1:size(nonLocSets,1);
predrows = setdiff(predrows,fixedrows);
FinalComps = cell(size(nonLocSets,1),1);
FinalComps(:) = {{}};
for r = 1:numel(fixedrows)
    crow = fixedrows(r);
    ccomps = cols(find(rows == crow));
    FinalComps(crow) = {CompIDs(ccomps)};
end

%start timers
FCTime = etime(clock,ctime)
%This is the prediction based on the preprocessing step
FCComps = FinalComps;
%save('PreStepSol','FinalComps');
nonLocSets(fixedrows,:) = [];
%and add the now localised reactions
localisedReactions = union(found,localisedReactions);
%if everything is done, return the solution.
if numel(nonLocSets) == 0
    return
end

fixedreactions = setdiff(orignonLocSets,nonLocSets);
%These reactions are removed and should carry 0 flux.
removedreacs = setdiff(fixedreactions,found);

%Set up the linear program
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
%save('TestData.mat','transpenalty','epsilon','t')
obj = [zeros(n,1);ones(q,1);epsilon*-1/t*transpenalty];

%We only have continous variables
coltype = '';
coltype(1,1:n) = 'C';
coltype(1,(end+1):(end+q+t)) = 'C';

csense = repmat('E',size(A,1),1);
csense(rhs == inf) = 'G';
csense(lhs == -inf) = 'L';
b = zeros(size(A,1),1);
b(rhs == inf) = lhs(rhs == inf);
b(lhs == -inf) = rhs(lhs == -inf);

%Set up the cplex model.
LPproblem.A = A;
LPproblem.lb = lbs;
LPproblem.ub = ubs;
LPproblem.csense = csense;
LPproblem.b = b;
LPproblem.c = obj;
LPproblem.osense = -1;

%Select all reactions which have already been localised and turn off those,
%which are not the correct ones
% Those reactions are the reactions, which were present in the original
% NonLocSets, and are present in neither the localised reactions nor the
% remaining non Localised reactions.
% so these reactions had to be in sets that 
removeReacs = setdiff(setdiff(orignonLocSets(orignonLocSets~=0),localisedReactions),nonLocSets(nonLocSets~=0));
%Set the reactions which were determined to be off to 0.
if numel(removeReacs) > 0
    LPproblem.lb(removeReacs) = 0;
    LPproblem.ub(removeReacs) = 0;
    model.lb(removeReacs) = 0;
    model.ub(removeReacs) = 0;
end

%unloc is the current nonLocSet
unloc = nonLocSets;
%locrows will indicate the rows that need to be localised in unloc
locrows = 1:size(unloc,1);
%stepsize is the frequency of using the pre-localisation step
stepsize = numel(unloc)/10;
step = 0;
%save(['FCLP' indicator '.mat'])
allscores = struct;
R2posc = find(ismember(model.rxns,'R2(c)'));
R2posm = find(ismember(model.rxns,'R2(m)'));
R2posp = find(ismember(model.rxns,'R2(p)'));


while numel(unloc) > 0    
    %get a random row 
    row = randi(numel(locrows));
    curr = locrows(row);
    %Calculate the scores for this row and determine the compartment for
    %this row.
    
    scores = calcScores(model,LPproblem,nonLocSets(curr,:),epsilon);         
    [comp,pos] = GetCompartmentsFromScore(scores, CompIDs, lp.Param.simplex.tolerances.feasibility.Cur);
    creac = model.rxns{nonLocSets(curr,1)};    
    creac = creac(1:end-3);    
    allscores.(['R' creac]) = scores;
    %add those reactions to the set of localised reactions and remove the
    %non localised reactions.
    localisedReactions = [ localisedReactions ; nonLocSets(curr,pos)'];
    removedreacs = setdiff(nonLocSets(curr,:),nonLocSets(curr,pos));
    %Eliminate non existing reactions.
    removedreacs = removedreacs(removedreacs~=0);    
    %And adjust the upper and lower bounds of the associated reactions
    %which are not in the "target" compartments.
    LPproblem.lb(removedreacs) = 0;
    LPproblem.ub(removedreacs) = 0;
    model.lb(removedreacs) = 0;
    model.ub(removedreacs) = 0;
    %reduce the unloc set and the locrows set
    index = true(1, size(unloc, 1));
    index(row) = 0;
    unloc = unloc(index,:);
    locrows(row) = [];
    %Assign the FinalComps for this row. Note, that the GetCompMethod is
    %suitable for multiple rows, and thus the first element (i.e. the only
    %element of the cell array) needs to be used.
    FinalComps{predrows(curr)} = comp{1};
    %And update the ReactionCompartmentalisation Data
    ReacLocal(nonLocToOrigMap) = FinalComps;
    step = step + 1;    
    %If we exceed the step-size repeat the pre-localisation step.
    if step > stepsize && step > 10
        step = 0;
        %calculate a new round of pre_localisation step
        [Metabolites,AvailableInComp,TransDirections] = determineUsedTransportersCOBRA(model , weightedRxns, nonLocSets,...
        localisedReactions, CompIDs, epsilon,...
        noncompmodel);
        [new_localised] = compartment_localisation_step_COBRA( model,CompIDs,...
            localisedReactions, ComparisonModel, unloc,...
            ReacLocal, epsilon,Metabolites,AvailableInComp,TransDirections);
        % Get all newly localised reactions
        nloc = intersect(new_localised,unloc);
        locdel = true(numel(locrows),1);
        %and update the respective fields
        for curreac =1:numel(nloc)
            [crow,ccol] = find(nonLocSets == nloc(curreac));
            FinalComps{predrows(crow)} = [FinalComps{predrows(crow)} , CompIDs{ccol}];
            [crow,ccol] = find(unloc == nloc(curreac));
            locdel(crow) = 0;
        end        
        localisedReactions = union(localisedReactions,new_localised);     
        unloc = unloc(locdel,:);       
        locrows = locrows(locdel);
        %And now update the linear Problem 
        %i.e. set all reactions which are remains of localised reactions to 0 (thus disallowing them).
        removeReacs = setdiff(nonLocSets,union(unloc,localisedReactions));
        removedreacs = removedreacs(removedreacs~=0);
        
        LPproblem.lb(removedreacs) = 0;
        LPproblem.ub(removedreacs) = 0;          
        model.lb(removedreacs) = 0;
        model.ub(removedreacs) = 0;
    end
end
%save('FCModpred','allscores');
end

function scores = calcScores(model, LPproblem,nonLocReacs,epsilon)


%fprintf('CalcScores1: The Bounds of R2(c), R2(m) and R2(p) are %f/%f , %f/%f , %f/%f\n',MILP.Model.lb(R2posc),MILP.Model.ub(R2posc),MILP.Model.lb(R2posm),MILP.Model.ub(R2posm),MILP.Model.lb(R2posp),MILP.Model.ub(R2posp))
scores = -inf * ones(size(nonLocReacs));
%turn off the alternate reactions.
oldlbs = zeros(size(nonLocReacs));
oldubs = zeros(size(nonLocReacs));
oldlbs(nonLocReacs~=0) = LPproblem.lb(nonLocReacs(nonLocReacs~=0));
oldubs(nonLocReacs~=0) = LPproblem.ub(nonLocReacs(nonLocReacs~=0));
LPproblem.lb(nonLocReacs(nonLocReacs~=0)) = 0;
LPproblem.ub(nonLocReacs(nonLocReacs~=0)) = 0;

for i=1:numel(nonLocReacs)
    if nonLocReacs(i) == 0
        %This reaction is not included so ignore it.
        continue
    end
    %force positive flux through the reaction, this has always to be
    %possible    
    score = -inf;
    oldlb = LPproblem.lb(nonLocReacs(i));           
    LPproblem.ub(nonLocReacs(i)) = oldubs(i);
    LPproblem.lb(nonLocReacs(i)) = epsilon;    
    sol = solveCobraLP(LPproblem);   
    try
        score = sol.obj;
    catch
 	if oldlbs(i) == 0
	        fprintf('FASTCOMP: No Solution found for Reaction %s\n',model.rxns{nonLocReacs(i)});
	end
        %printRxnFormula(model,model.rxns{nonLocReacs(i)});
        %save(model.rxns{nonLocReacs(i)},'MILP')
    end
    LPproblem.lb(nonLocReacs(i)) = oldlb;    
    LPproblem.ub(nonLocReacs(i)) = 0;    
    if model.lb(nonLocReacs(i)) < 0
        oldub = LPproblem.ub(nonLocReacs(i));            
        LPproblem.ub(nonLocReacs(i)) = -epsilon;        
        LPproblem.lb(nonLocReacs(i)) = oldlbs(i);        
        sol = solveCobraLP(LPproblem);        
        try            
            if sol.status == 1
                score = max(score,sol.obj);
            end
        catch
		if score == -inf
	            fprintf('FASTCOMP: No Solution found for positive Reaction %s\n',model.rxns{nonLocReacs(i)});
		end
        end
        LPproblem.ub(nonLocReacs(i)) = oldub;
        LPproblem.lb(nonLocReacs(i)) = 0;
    end   
    scores(i) = score;
end
LPproblem.lb(nonLocReacs(nonLocReacs~=0)) = oldlbs(nonLocReacs~=0);
LPproblem.ub(nonLocReacs(nonLocReacs~=0)) = oldubs(nonLocReacs~=0);

end



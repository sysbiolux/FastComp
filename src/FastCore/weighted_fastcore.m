function  [A,score] = weighted_fastcore( C, model, epsilon, weightedRxns, weightratio, listreacs ) 
%
% A = fastcore( C, model, epsilon )
%
% The FASTCORE algorithm for context-specific metabolic network reconstruction
% Input C is the core set, and output A is the reconstruction

% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

if nargin < 6
    listreacs = [];
end

model_org = model;

N = 1:numel(model.rxns);
I = find(model.rev==0);

A = [];
flipped = false;
singleton = false;  

% start with I
J = intersect( C, I );% fprintf('|J|=%d  ', length(J));
P = setdiff( N, C);
Supp = weighted_findSparseMode( J, P, singleton, model, epsilon, weightedRxns, weightratio );
if ~isempty( setdiff( J, Supp ) ) 
  fprintf ('Error: Inconsistent irreversible core reactions.\n');
  %We need an error value here, otherwise we have no clue that there was an
  %error. 
  A=-1;
  return;
end
A = Supp; % fprintf('|A|=%d\n', length(A));
%printRxnFormula(model,model.rxns(intersect(A,listreacs)));
J = setdiff( C, A );% fprintf('|J|=%d  ', length(J));
%disp('---------------------------Starting Loop -------------------------\n')
% main loop     
while ~isempty( J )
    P = setdiff( P, A);
    Supp = weighted_findSparseMode( J, P, singleton, model, epsilon, weightedRxns, weightratio );
    A = union( A, Supp );  % fprintf('|A|=%d\n', length(A)); 
    if ~isempty( intersect( J, A ))
        J = setdiff( J, A );    % fprintf('|J|=%d  ', length(J));
        flipped = false;
    else
        if singleton
            JiRev = setdiff(J(1),I);
        else
            JiRev = setdiff(J,I);
        end
        if flipped || isempty( JiRev )
            if singleton
                fprintf('\nError: Global network is not consistent.\n');
                A = -1;
                return
            else
              flipped = false;
              singleton = true;
            end
        else
            model.S(:,JiRev) = -model.S(:,JiRev);
            tmp = model.ub(JiRev);
            model.ub(JiRev) = -model.lb(JiRev);
            model.lb(JiRev) = -tmp;
            flipped = true; % fprintf('(flip)  ');
        end
    end
end
%printRxnFormula(model,model.rxns(intersect(A,listreacs)));


score = 1* numel(setdiff(A,weightedRxns)) + weightratio * numel(intersect(A,weightedRxns));
%fprintf('|A|=%d\n', length(A));




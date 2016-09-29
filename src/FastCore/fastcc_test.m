function [A,model,incI,V] = fastcc_test( model, epsilon ) 
%
% A = fastcc( model, epsilon )
%
% Amodified version of the fastcc method by Vlassis,Pires Pacheco and
% Sauter testing whether there are reactions which cannot carry flux in the
% forward direction simultaneously and flipping their directionality.
% 
% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg
% Modification by Thomas Pfau, 2015, LSRU University of Luxembourg


tic

N = (1:numel(model.rxns));

A = [];

% start with all reactions and see whether they can all carry positive flux at the
% same time.
J = N ; % intersect( N, I ); fprintf('|J|=%d  ', numel(J));
V = LP7( J, model, epsilon ); 
Supp = find( V >= 0.99*epsilon );
[model,flips] = flipRxns(model,find(V < -0.99 *epsilon));
%flips
Supp = union(Supp,flips);
A = Supp;  %fprintf('|A|=%d\n', numel(A));
incI = setdiff( J, A );    
if ~isempty( incI )
    [model,flips] = flipRxns(model,incI);
end

end

function [model,flipped] = flipRxns(model,rxns)
            model.S(:,rxns) = -model.S(:,rxns);
            tmp = model.ub(rxns);
            model.ub(rxns) = -model.lb(rxns);
            model.lb(rxns) = -tmp;
            %fprintf('(flip)  %i\n', numel(rxns));            
            flipped =rxns;
end
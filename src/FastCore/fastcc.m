function [A] = fastcc( model, epsilon, ignore ) 
%
% A = fastcc( model, epsilon )
%
% The FASTCC algorithm for testing the consistency of an input model
% Output A is the consistent part of the model

% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

%tic
%N = (1:numel(model.rxns));
inact = model.lb~=0 | model.ub~=0;
N = find(inact);

if nargin > 2
    N = setdiff(N,ignore);
end


I = intersect(N,find(model.rev==0));

A = [];

% start with I
J = intersect( N, I ); %fprintf('|J|=%d  ', numel(J));
V = LP7( J, model, epsilon ); 
Supp = find( abs(V) >= 0.99*epsilon );  
A = Supp;  %fprintf('|A|=%d\n', numel(A));
incI = setdiff( J, A );    
if ~isempty( incI )
    %fprintf('\n(inconsistent subset of I detected)\n');
end
J = setdiff( setdiff( N, A ), incI);  %fprintf('|J|=%d  ', numel(J));

% reversible reactions
flipped = false;
singleton = false;     
firstsingle = true;
while ~isempty( J )
    if singleton
        if firstsingle
            firstsingle = false;
        end
        Ji = J(1);
        V = LP3( Ji, model ) ; 
    else
        Ji = J;
        V = LP7( Ji, model, epsilon ) ; 
    end    
    Supp = find( abs(V) >= 0.99*epsilon );  
    A = union( A, Supp); % fprintf('|A|=%d\n', numel(A)); 
    if ~isempty( intersect( J, A ))
        J = setdiff( J, A );%     fprintf('|J|=%d  ', numel(J));
        flipped = false;
    else
        JiRev = setdiff( Ji, I );
        if flipped || isempty( JiRev )
            flipped = false;
            if singleton
                J = setdiff( J, Ji );  
                %fprintf('\n(inconsistent reversible reaction detected)\n');
            else
                singleton = true;
 %               toc
 %               tic
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

if numel(A) == numel(N)
    %fprintf('\nThe input model is consistent.\n'); 
end

%toc

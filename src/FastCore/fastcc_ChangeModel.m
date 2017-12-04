function [A,model,incI] = fastcc_ChangeModel( model, epsilon ) 
%
% A = fastcc( model, epsilon )
%
% This version of fastcc will, if possible return a model where all 
% consistent reactions can simultaneously carry flux in the forward
% direction.
%
% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg
% Modifications by Thomas Pfau

N = (1:numel(model.rxns));
I = find(model.lb < 0);

A = [];
Vs = {};
% start with I
J = intersect( N, I ); %fprintf('|J|=%d  ', numel(J));
V = LP7( J, model, epsilon ); 
Supp = find( V >= 0.99*epsilon );
[model,flips] = flipRxns(model,find(V < -0.99 *epsilon));
Supp = union(Supp,flips);
A = Supp;  %fprintf('|A|=%d\n', numel(A));
incI = setdiff( J, A );    
if ~isempty( incI )
    %fprintf('\n(inconsistent subset of I detected)\n');
end
J = setdiff( setdiff( N, A ), incI);  %fprintf('|J|=%d  ', numel(J));

flipedRs = flips;

% reversible reactions
flipped = false;
singleton = false;     
firstsingle = true;
while ~isempty( J )
    if singleton
        Ji = J(1);
        V = LP3( Ji, model ) ; 
    else
        Ji = J;
        V = LP7( Ji, model, epsilon ) ; 
    end    
    %Lets try to flip the reactions which are used 
    Supp = intersect(J,find( V >= 0.99*epsilon )); 
    
    A = union( A, Supp); % fprintf('|A|=%d\n', numel(A)); 
    if ~isempty( intersect( J, A ))
        J = setdiff( J, A );     %fprintf('|J|=%d  ', numel(J));
        flipped = false;
    else
        JiRev = setdiff( Ji, I );
        if flipped || isempty( JiRev )
            flipped = false;
            if singleton
                J = setdiff( J, Ji );  
                %fprintf('\n(inconsistent reversible reaction detected)\n');
            else
               % fprintf('Singelton Step\n');
                singleton = true;
            end
        else
            model.S(:,JiRev) = -model.S(:,JiRev);
            tmp = model.ub(JiRev);
            model.ub(JiRev) = -model.lb(JiRev);
            model.lb(JiRev) = -tmp;
            flipped = true;% fprintf('(flip)  %i\n', numel(JiRev));            
            
        end
    end
end

if numel(A) == numel(N)
    %fprintf('\nThe input model is consistent.\n'); 
end

[A,model,incI] = fastcc_test(model,epsilon);
count = 0;
while numel(incI) > 0 
    %fprintf('Adjusting %i reactions in model\n',numel(incI));
    if(mod(count,50) == 0)
        fprintf('There are %i reactions remaining\n',numel(incI));
    end
    model = flipRxns(model,incI(randperm(numel(incI),randi(numel(incI)))));    
    [A,model,incI] = fastcc_test(model,epsilon);            
    count = count +1;
    
%    if count > 20
%        model2 = removeRxns(model,model.rxns(incI));
%        A = fastcc_test(model2,epsilon);
%        if isempty(setdiff(A,1:numel(model2.rxns)))
%            fprintf('The following reactions are likely inconsistent:\n')
%            model.rxns(incI)
            %model = model2;
%            return
%        end
%    end
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
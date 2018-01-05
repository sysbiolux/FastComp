function sol = getSolForAllP(model,P)
    model.c(:) = 0;
    model.lb(P) = 1;
    try
        sol = optimizeCbModel(model,'max','one');
    catch
        disp('ups')
    end
    if ~isfield(sol,'v') && ~isempty(sol.full)
        sol.v = sol.full(1:numel(model.rxns));
    end
end
function [comps,locs] = GetCompartmentsFromScore(scores, Compartments, threshold)
comps = {};
maxscores = max(scores,[],2);
for i = 1:numel(maxscores)
    locs = find(abs(scores(i,:)-maxscores(i)) < threshold);        
    if numel(locs) > 1
        locs(locs==1) = [];
    end
    if numel(locs) <=3
        comps{i} = Compartments(locs);        
    else
        comps{i} = {};
    end
    if numel(locs) == 0
        comps{i} = {};
    end
end
comps = comps';
end
function Exchangers = createExchangerTableFromModelAndType(compModel,decompModel)
%Assume model is decompartmentalised (everything in c)

Exchangers = cell(0,5);
CompEx = find(findExcRxns(compModel));
DecompEx = findExcRxns(decompModel);
%first, merge all exchangers.
NewInfo = cell(numel(CompEx),3);
notpresent = false(size(CompEx));
for i = 1:numel(CompEx)
    cmet = compModel.mets{find(compModel.S(:,CompEx(i)))};
    metName = regexprep(cmet,'\[[a-z]\]$','');
    metComp = regexprep(cmet,'.*\[([a-z])\]$','$1');
    if ~any(ismember(decompModel.mets,[metName,'[c]']))
        notpresent(i) = true;
        continue;
    end
    lb = compModel.lb(CompEx(i));
    ub = compModel.ub(CompEx(i));
    metcoef = sum(compModel.S(:,CompEx(i)));
    if lb < 0 && ub > 0
        %this is easy, its a free exchanger. directionality is
        %irrelevant.
        NewInfo(i,:) = {regexprep(cmet,'\[[a-z]\]$',''),regexprep(cmet,'.*\[([a-z])\]$','$1'),0};
    else %So, either lb >=0  or ub <=0
        if ub <=0
           if metcoef > 0
               %this is an exporter
               NewInfo(i,:) = {regexprep(cmet,'\[[a-z]\]$',''),regexprep(cmet,'.*\[([a-z])\]$','$1'),1};
           else %this is an importer
               NewInfo(i,:) = {regexprep(cmet,'\[[a-z]\]$',''),regexprep(cmet,'.*\[([a-z])\]$','$1'),-1};
           end
        else %lb >=0
           if metcoef > 0 %this is an importer
               NewInfo(i,:) = {regexprep(cmet,'\[[a-z]\]$',''),regexprep(cmet,'.*\[([a-z])\]$','$1'),-1};
           else
               NewInfo(i,:) = {regexprep(cmet,'\[[a-z]\]$',''),regexprep(cmet,'.*\[([a-z])\]$','$1'),1};
           end
        end
    end
end
NewInfo(notpresent,:) = [];
ModelEx = sum(decompModel.S~=0) == 1;
for i = 1:size(NewInfo)
    cmet = [NewInfo{i,1} '[c]'];
    metPos = ismember(decompModel.mets,cmet);
    ExchangerPos = (decompModel.S(metPos,:)~=0) & ModelEx;
    metProd = decompModel.S(metPos,ExchangerPos) > 0;
    if metProd
        metCoeff = 1;
    else
        metCoeff = -1;
    end
    %Now, we can have importers, exporters or both directionalities.
    if NewInfo{i,3} == 0
        %This is simple, we just add the information from the reaction, as
        %it can carry flux
        Exchangers(i,:) = {NewInfo{i,1:2},metCoeff,decompModel.lb(ExchangerPos),decompModel.ub(ExchangerPos)};
    elseif NewInfo{i,3} == 1
        %This is an exporter. So lets see
            if metCoeff == 1
                %now, this is an issue. But it indicates, 
                %that there are other reactions which are actually
                %importers for the metabolite. So in the end, we just add a
                %inverse reaction, and we will have a working flux through
                %transporters.
                Exchangers(i,:) = {NewInfo{i,1:2},-metCoeff,...
                    0,max(decompModel.ub)};
            else
                Exchangers(i,:) = {NewInfo{i,1:2},metCoeff,...
                    0,max(decompModel.ub)};
            end
    else
        %This is an importer
        
        if metCoeff == -1
            %As above this is an issue. But it indicates,
            %that there are other reactions which are actually
            %exporters for the metabolite. So in the end, we just add a
            %inverse reaction, and we will have a working flux through
            %transporters.
            Exchangers(i,:) = {NewInfo{i,1:2},-metCoeff,...
                0,max(decompModel.ub)};
        else
            Exchangers(i,:) = {NewInfo{i,1:2},metCoeff,...
                0,max(decompModel.ub)};
        end
    end
end

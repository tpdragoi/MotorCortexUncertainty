function trialNums = findTrials(obj, conditions)

% % older data objects have obj.bp.autowater.nums
% % newer data objects do not (from summer2021 pipeline)
% if ~isfield(obj.bp.autowater, 'nums')
%     tmp = obj.bp.autowater;
%     obj.bp = rmfield(obj.bp, 'autowater');
%     obj.bp.autowater.nums = tmp + (tmp-1)*-2;
% end

% if strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
% 
%     obj.bp.Ntrials = 298;
% 
% end


varnames = getStructVarNames(obj);
Ntrials = obj.bp.Ntrials;
for i = 1:numel(varnames)
    eval([varnames{i} ' = obj.bp.' varnames{i} ';']);
    
    if eval(['numel(' varnames{i} ')==obj.bp.Ntrials && isrow(' varnames{i} ')'])
        eval([varnames{i} '=' varnames{i} ''';']);
    end
    
    eval([varnames{i} '=' varnames{i} '(~isnan(' varnames{i} '))' ''';'])
end



if strcmp(obj.pth.dt,'2024-07-07') && strcmp(obj.pth.anm,'TD9si')

else

if isfield (obj.bp, 'fidx')
    obj.bp.Ntrials = sum(obj.bp.fidx==2);
    addToTrialNums = sum(obj.bp.fidx==1);
end
end


if strcmp(obj.pth.dt,'2024-08-25') && strcmp(obj.pth.anm,'TD7d')

    obj.bp.Ntrials = 340;

end

if strcmp(obj.pth.dt,'2024-12-03') && strcmp(obj.pth.anm,'TD16d')

    obj.bp.Ntrials = 288;

end

if strcmp(obj.pth.dt,'2025-07-14') && strcmp(obj.pth.anm,'TD25d')

    obj.bp.Ntrials = 171;
    addToTrialNums = 0;

end


mask = zeros(obj.bp.Ntrials, numel(conditions));

for i = 1:numel(conditions)
    try
        mask(:,i) = eval(conditions{i});
    catch
        mask(:,i) = eval(conditions{i}{1});
    end
    if exist('addToTrialNums','var') == 1
        trialNums{i} = find(mask(:,i)) + addToTrialNums;
    else
        trialNums{i} = find(mask(:,i));
    end
end


end % findTrials
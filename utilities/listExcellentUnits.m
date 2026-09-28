%% listExcellentUnits.m
%  Lists every unit labeled 'excellent' in <dataRoot>/<task>/*_obj.mat, and the
%  number of good / excellent / other units in each session.
%  Reads objFL only (clu lives there). Unmapped sessions are included, marked
%  as such. Behaviour-only sessions are skipped.

clear; clc

repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
addpath(repoRoot);
cfgPaths = setPaths();
dataRoot = cfgPaths.dataRoot;

files = dir(fullfile(dataRoot, '*', '*_obj.mat'));   % task folders one level down
exc  = cell(0, 8);   % task, anm, date, region, probe, cluid, cluidOriginal, label
sess = cell(0, 7);   % task, anm, date, regions, nGood, nExcellent, nOther
tAll = tic;

for i = 1:numel(files)
    f = fullfile(files(i).folder, files(i).name);
    [~, task] = fileparts(files(i).folder);
    w = whos('-file', f);  vnames = {w.name};
    if ismember('objFL', vnames)
        L = load(f, 'objFL');  O = L.objFL;
    else
        L = load(f, 'objFL300');  O = L.objFL300;
    end
    clear L
    regStr = char(O.regions);
    if strcmp(regStr, 'behaviour only') || isempty(O.probeOfUnit)
        continue
    end

    % quality of every unit, in unit-axis order (cluid values index into clu{k})
    if iscell(O.cluid), cells = O.cluid(:)'; else, cells = {O.cluid}; end
    if iscell(O.clu), clu = O.clu; else, clu = {O.clu}; end
    if isfield(O, 'cluidOriginal')
        if iscell(O.cluidOriginal), cellsOrig = O.cluidOriginal(:)'; else, cellsOrig = {O.cluidOriginal}; end
    else
        cellsOrig = cells;
    end
    qRaw = cell(0,1);  cid = zeros(0,1);  cidOrig = zeros(0,1);
    for k = 1:numel(cells)
        if numel(clu) == numel(cells), cK = clu{k}; else, cK = clu{1}; end
        labs = {cK.quality}';
        for j = 1:numel(labs)
            if ~ischar(labs{j}) && ~isstring(labs{j}), labs{j} = ''; end
        end
        qRaw    = [qRaw; labs(cells{k}(:))];            %#ok<AGROW>
        cid     = [cid; cells{k}(:)];                   %#ok<AGROW>
        cidOrig = [cidOrig; cellsOrig{k}(:)];           %#ok<AGROW>
    end
    q = lower(strtrim(qRaw));
    probeOf = O.probeOfUnit(:);

    % region of each unit: nonzero groupCodes, in order, match the '+'-separated regions
    regs = strsplit(regStr, '+');
    gc = O.groupCodes(:)';  probes = gc(gc > 0);
    unitReg = repmat({regStr}, numel(probeOf), 1);      % 'unmapped' sessions keep that label
    if numel(regs) == numel(probes)
        for r = 1:numel(regs)
            unitReg(probeOf == probes(r)) = regs(r);
        end
    end

    isG = strcmp(q, 'good');  isE = strcmp(q, 'excellent');
    anm = char(O.anm);  dt = char(O.date);
    sess(end+1, :) = {task, anm, dt, regStr, sum(isG), sum(isE), sum(~isG & ~isE)}; %#ok<SAGROW>
    for u = find(isE)'
        exc(end+1, :) = {task, anm, dt, unitReg{u}, probeOf(u), cid(u), cidOrig(u), qRaw{u}}; %#ok<SAGROW>
    end
    fprintf('%-11s %-26s %-10s good %4d  excellent %3d  other %3d\n', task, files(i).name, regStr, ...
        sum(isG), sum(isE), sum(~isG & ~isE));
    clear O
end

S = cell2table(sess, 'VariableNames', {'task','anm','date','regions','nGood','nExcellent','nOther'});
E = cell2table(exc,  'VariableNames', {'task','anm','date','region','probe','cluid','cluidOriginal','label'});

fprintf('\n==== %d excellent units in %d of %d sessions ====\n', height(E), sum(S.nExcellent > 0), height(S));
if height(E) > 0, disp(E); end
fprintf('==== Sessions with excellent units ====\n');
disp(S(S.nExcellent > 0, :))
fprintf('==== Excellent units per task x region ====\n');
if height(E) > 0, disp(groupsummary(E, {'task','region'})); end
fprintf('%.1f min\n', toc(tAll)/60);

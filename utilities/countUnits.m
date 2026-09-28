%% countUnits.m
%  Counts neurons per task and region, two ways:
%    good     quality label 'good' only
%    goodExc  'good' + 'excellent' (what every neural analysis uses)
%  Each is reported with every stored unit and after the analyses' lowFR
%  criterion: mean rate > 0.01 Hz over all trials, goCue-aligned -2..4 s
%  (as in shared\slimToLegacy).
%
%  Reads only <dataRoot>/<task>/*_obj.mat (dataRoot is set in setPaths.m).
%  Region per unit comes from obj.regions + obj.groupCodes (the probe each
%  region used). Sessions with regions 'unmapped' or 'behaviour only' are
%  skipped. The analyses do not use them either.
%
%  Output: printed tables, plus unitCounts_perSession.csv and
%  unitCounts_summary.csv in the output folder set in setPaths.m.
%
%  lowFR is computed on the 100 Hz data. The rate over a 6 s window equals
%  spikes / (6 s x trials) at any bin size, so it matches the 200/300 Hz
%  scripts except for spikes that fall exactly on a window edge.

clear; clc

repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
addpath(repoRoot);
cfgPaths = setPaths();
dataRoot = cfgPaths.dataRoot;
outDir   = cfgPaths.outputRoot;

tasks     = {'R1','R14','R16','Learning','VTA'};
taskNames = {'Simple Reward','Delayed Reward','Double Reward','Learning','VTA'};
expectedN = [22 35 26 20 18];            % sessions per task, from sessionCounts.m
lowFR     = 0.01;                        % Hz
frWin     = [-2 4];                      % s, goCue-aligned
regionMap = {'M1','tjM1'; 'S1','tjS1'; 'ALM','ALM'};   % file label -> paper label

rows = {};   % task, anm, date, region, nAll, nGood, nGoodExc, nGoodFR, nGoodExcFR
nSess = zeros(1, numel(tasks));
tAll = tic;

for t = 1:numel(tasks)
    files = dir(fullfile(dataRoot, tasks{t}, '*_obj.mat'));
    for i = 1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        w = whos('-file', f);  vnames = {w.name};

        % ---- firstLick struct: units, quality, regions --------------------
        if ismember('objFL', vnames)
            L = load(f, 'objFL');  O = L.objFL;
        else                                            % one-struct-per-rate layout
            L = load(f, 'objFL300');  O = L.objFL300;
        end
        clear L
        regStr = char(O.regions);
        if strcmp(regStr, 'unmapped') || strcmp(regStr, 'behaviour only') || isempty(O.probeOfUnit)
            fprintf('%-9s %-28s skipped (%s)\n', tasks{t}, files(i).name, regStr);
            clear O
            continue
        end
        nSess(t) = nSess(t) + 1;

        % unit axis = cluid cells concatenated; cluid values index into clu{k}
        if iscell(O.cluid), cells = O.cluid(:)'; else, cells = {O.cluid}; end
        if iscell(O.clu), clu = O.clu; else, clu = {O.clu}; end
        q = cell(0, 1);
        for k = 1:numel(cells)
            if numel(clu) == numel(cells), cK = clu{k};
            elseif numel(clu) == 1,        cK = clu{1};
            else, error('%s: clu has %d entries for %d cluid cells.', files(i).name, numel(clu), numel(cells));
            end
            labs = {cK.quality}';
            for j = 1:numel(labs)
                if ~ischar(labs{j}) && ~isstring(labs{j}), labs{j} = ''; end
            end
            labs = lower(strtrim(labs));
            q = [q; labs(cells{k}(:))];                 %#ok<AGROW>
        end
        probeOf = O.probeOfUnit(:);
        assert(numel(q) == numel(probeOf), '%s: %d quality labels for %d units.', files(i).name, numel(q), numel(probeOf));
        isGood    = strcmp(q, 'good');
        isGoodExc = isGood | strcmp(q, 'excellent');
        other = unique(q(~isGoodExc));
        if ~isempty(other)
            fprintf('   %s: %d units with other labels: %s\n', files(i).name, sum(~isGoodExc), strjoin(other, ', '));
        end

        % region -> probe: nonzero groupCodes, in order, match the '+'-separated regions
        regs   = strsplit(regStr, '+');
        gc     = O.groupCodes(:)';
        probes = gc(gc > 0);
        assert(numel(regs) == numel(probes), '%s: regions ''%s'' but groupCodes %s.', files(i).name, regStr, mat2str(gc));
        if isfield(O, 'groupLabels')                    % cross-check when present
            assert(isequal(regs(:), O.groupLabels(gc > 0)'), '%s: regions and groupLabels disagree.', files(i).name);
        end
        anm = char(O.anm);  dt = char(O.date);
        clear O

        % ---- goCue struct: lowFR criterion --------------------------------
        if ismember('objGC', vnames)
            L = load(f, 'objGC');  tm = L.objGC.time100(:);  X = L.objGC.trialdat100;
        else
            L = load(f, 'objGC100');  tm = L.objGC100.time(:);  X = L.objGC100.trialdat;
        end
        clear L
        assert(size(X, 2) == numel(probeOf), '%s: goCue export has %d units, firstLick %d.', files(i).name, size(X, 2), numel(probeOf));
        rIdx = tm > frWin(1) & tm < frWin(2);           % bin centers inside the window
        fr   = squeeze(mean(mean(double(X(rIdx, :, :)), 3), 1));
        fr   = fr(:);
        clear X
        passFR = fr > lowFR;

        % ---- one row per session x region --------------------------------
        for r = 1:numel(regs)
            reg = regs{r};
            if strcmp(tasks{t}, 'Learning'), reg = 'tjM1'; end      % 'Day N' groups, all tjM1
            hit = strcmp(regionMap(:, 1), reg);
            if any(hit), reg = regionMap{hit, 2}; end
            m = probeOf == probes(r);
            rows(end+1, :) = {taskNames{t}, anm, dt, reg, sum(m), sum(m & isGood), sum(m & isGoodExc), ...
                sum(m & isGood & passFR), sum(m & isGoodExc & passFR)}; %#ok<SAGROW>
        end
        unassigned = ~ismember(probeOf, probes);
        if any(unassigned)
            fprintf('   %s: %d units on probes not mapped to a region (not counted).\n', files(i).name, sum(unassigned));
        end
        fprintf('%-9s %-28s %-10s %4d units (%d good, %d good+exc)\n', tasks{t}, files(i).name, regStr, ...
            numel(probeOf), sum(isGood), sum(isGoodExc));
    end
end

T = cell2table(rows, 'VariableNames', {'task','anm','date','region','nAll','good','goodExc','good_lowFR','goodExc_lowFR'});
T.task   = categorical(T.task, taskNames);
T.region = categorical(T.region, unique([{'tjM1','ALM','tjS1'}, T.region'], 'stable'));
vars = {'good','goodExc','good_lowFR','goodExc_lowFR'};

% ---- session-count check against sessionCounts.m --------------------------
fprintf('\nSessions counted vs sessionCounts.m:\n');
for t = 1:numel(tasks)
    flag = '';
    if nSess(t) ~= expectedN(t), flag = '   <-- MISMATCH'; end
    fprintf('  %-15s %3d  (expected %d)%s\n', taskNames{t}, nSess(t), expectedN(t), flag);
end

% ---- summaries ------------------------------------------------------------
byTaskRegion = groupsummary(T, {'task','region'}, 'sum', vars);
byRegion     = groupsummary(T, 'region', 'sum', vars);
byTask       = groupsummary(T, 'task', 'sum', vars);
byTaskRegion.Properties.VariableNames{'GroupCount'} = 'nSessions';
byRegion.Properties.VariableNames{'GroupCount'}     = 'nSessionRegionPairs';
byTask.Properties.VariableNames{'GroupCount'}       = 'nSessionRegionPairs';
total = varfun(@sum, T, 'InputVariables', vars);

fprintf('\n==== Units per task per region ====\n');       disp(byTaskRegion)
fprintf('==== Units per region, all tasks ====\n');       disp(byRegion)
fprintf('==== Units per task, all regions ====\n');       disp(byTask)
fprintf('==== All tasks, all regions ====\n');            disp(total)

% ---- save -----------------------------------------------------------------
if ~exist(outDir, 'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir, 'unitCounts_perSession.csv'));
S1 = byTaskRegion;  S1.level = repmat({'task x region'}, height(S1), 1);
S2 = byRegion;      S2.task  = categorical(repmat({'<all>'}, height(S2), 1));  S2.level = repmat({'region'}, height(S2), 1);
S3 = byTask;        S3.region = categorical(repmat({'<all>'}, height(S3), 1)); S3.level = repmat({'task'}, height(S3), 1);
S2.Properties.VariableNames{'nSessionRegionPairs'} = 'nSessions';
S3.Properties.VariableNames{'nSessionRegionPairs'} = 'nSessions';
S4 = total;  S4.task = categorical({'<all>'});  S4.region = categorical({'<all>'});  S4.nSessions = height(T);  S4.level = {'total'};
cols = [{'level','task','region','nSessions'}, strcat('sum_', vars)];
S4.Properties.VariableNames(1:numel(vars)) = strcat('sum_', vars);
summary = [S1(:, cols); S2(:, cols); S3(:, cols); S4(:, cols)];
summary.task = cellstr(summary.task);  summary.region = cellstr(summary.region);
writetable(summary, fullfile(outDir, 'unitCounts_summary.csv'));
fprintf('Saved unitCounts_perSession.csv and unitCounts_summary.csv to %s  (%.1f min)\n', outDir, toc(tAll)/60);

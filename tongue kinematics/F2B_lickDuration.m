%% F2B_lickDuration.m
%  Duration of successive tongue protrusions, Delayed Reward Task.
%  Protrusions are detected from tongue length, summarized per session, then
%  averaged across sessions.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\loadBehavSession; Data = dataRoot in setPaths.m
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/200
%    params.smooth      1
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2.5 to 4 s
%  Run the whole file.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


% ---- data folder from setPaths.m ----
% Data\R14 holds the exported <ANM>_<DATE>_obj.mat / _kin.mat files.
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
assert(exist(fullfile(repoRoot, 'setPaths.m'), 'file') == 2, ...
    'Cannot find setPaths.m. Run this script from its file, or cd to the repository root first.');
addpath(repoRoot);
cfgPaths = setPaths();   % data locations are set once, in setPaths.m
dataDir = fullfile(cfgPaths.dataRoot, 'R14');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s (set dataRoot in setPaths.m)', dataDir);

%% PARAMETERS
params.alignEvent = 'goCue';
params.behav_only = 1;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1'};

params.tmin   = -2.5;
params.tmax   = 4;
params.dt     = 1/200;
params.smooth = 1;
params.quality = {'good'};   % good units only (findClusters trims blanks and ignores case)

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'}, ...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};
params.feat_varToExplain = 80;
params.N_varToExplain    = 80;
params.advance_movement  = 0;
params.fcut   = 10;
params.cond   = 5;
params.method = 'xcorr';
params.fa     = false;
params.bctype = 'reflect';

%% CONFIG
% BOUT DETECTION. Every threshold is a TIME; sample counts are derived from
% params.dt, so changing the sampling rate cannot silently change what counts
% as a lick.
cfg.nanBridge_sec   = 0.020;   % interpolate tracking gaps up to this long
cfg.minBout_sec     = 0.030;   % a contact run shorter than this is not a lick
cfg.maxBout_sec     = 0.400;   % runs longer than this after the contact-driven
% split are artifacts, not licks
cfg.contactTol_sec  = 0.020;   % a run must fall within this of a lick event
cfg.maxLickIdx      = 35;   % lick-indexed columns kept per trial
cfg.nLicksSummary   = 15;   % licks carried into the summary matrices

% ANALYSIS / FIGURE
cfg.nLicksPlot  = 8;   % licks shown and tested
cfg.alpha       = 0.05;
cfg.plotChoices = {'mean','sd'};   % both measures are plotted, one figure each
cfg.tail        = 'right';   % unused in this script (no test)
cfg.durScale    = 1000;   % seconds -> ms for display
cfg.durUnit     = 'ms';

% TRIAL CAPS: {animal, date, last usable trial}
cfg.exclusions = {
    'TD13d', '2024-11-11', 278;
    'TD8d',  '2024-09-07', 313;
    'TD8d',  '2024-09-09', 298;
};

%% SESSIONS TO LOAD
%% Expected for this panel: n = 36 sessions, 7 animals

% one empty placeholder per session slot; the ones a task uses are
% filled in below and the rest drop out of the all_meta concatenation
[meta, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11, ...
    meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19, meta20, meta21, ...
    meta22, meta23, meta24, meta25, meta26, meta27, meta28, meta29, meta30, meta31, ...
    meta32, meta33, meta34, meta35, meta36, meta37, meta38] = deal([]);

% TD1d
date = '2023-02-21';
meta1 = struct('anm','TD1d','date',date);
date = '2023-02-22';
meta2 = struct('anm','TD1d','date',date);
date = '2023-02-23';
meta3 = struct('anm','TD1d','date',date);
date = '2023-02-24';
meta4 = struct('anm','TD1d','date',date);

% TD4d
date = '2023-02-21';
meta5 = struct('anm','TD4d','date',date);
date = '2023-02-24';
meta6 = struct('anm','TD4d','date',date);
date = '2023-02-25';
meta7 = struct('anm','TD4d','date',date);
date = '2023-03-19';
meta8 = struct('anm','TD4d','date',date);

% TD13d
date = '2024-11-11';
meta9 = struct('anm','TD13d','date',date);
date = '2024-11-12';
meta11 = struct('anm','TD13d','date',date);
date = '2024-11-13';
meta12 = struct('anm','TD13d','date',date);
date = '2024-11-21';
meta13 = struct('anm','TD13d','date',date);
date = '2024-11-22';
meta14 = struct('anm','TD13d','date',date);
date = '2024-11-24';
meta15 = struct('anm','TD13d','date',date);
date = '2024-11-25';
meta16 = struct('anm','TD13d','date',date);

% TD15d
date = '2024-11-24';
meta17 = struct('anm','TD15d','date',date);
date = '2024-11-25';
meta18 = struct('anm','TD15d','date',date);
date = '2024-11-26';
meta19 = struct('anm','TD15d','date',date);
date = '2024-11-27';
meta20 = struct('anm','TD15d','date',date);

% TD8d
date = '2024-09-06';
meta21 = struct('anm','TD8d','date',date);
date = '2024-09-07';
meta22 = struct('anm','TD8d','date',date);
date = '2024-09-08';
meta23 = struct('anm','TD8d','date',date);
date = '2024-09-09';
meta24 = struct('anm','TD8d','date',date);
date = '2024-09-10';
meta25 = struct('anm','TD8d','date',date);
date = '2024-09-22';
meta26 = struct('anm','TD8d','date',date);

date = '2025-06-17';
meta27 = struct('anm','TD22d','date',date);
date = '2025-06-18';
meta28 = struct('anm','TD22d','date',date);
date = '2025-06-19';
meta29 = struct('anm','TD22d','date',date);
date = '2025-06-20';
meta30 = struct('anm','TD22d','date',date);
date = '2025-06-21';
meta31 = struct('anm','TD22d','date',date);

date = '2025-06-17';
meta32 = struct('anm','TD23d','date',date);
date = '2025-06-18';
meta33 = struct('anm','TD23d','date',date);
date = '2025-06-19';
meta34 = struct('anm','TD23d','date',date);
date = '2025-06-20';
meta35 = struct('anm','TD23d','date',date);
date = '2025-06-21';
meta36 = struct('anm','TD23d','date',date);

y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12;meta13 ...
    ;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24;meta25;meta26;meta27;meta28;meta29;meta30; ...
    meta31;meta32;meta33;meta34;meta35;meta36];

%% MAIN LOOP
nSess    = numel(all_meta);
allDur   = cell(1, nSess);   allDurSD = cell(1, nSess);
allLen   = cell(1, nSess);   allILI   = cell(1, nSess);
allAnm   = cell(1, nSess);   allDate  = cell(1, nSess);
allDur4  = cell(1, nSess);   allDurSD4 = cell(1, nSess);
allLen4  = cell(1, nSess);   allILI4   = cell(1, nSess);

for sessnum = 1:nSess

    clear obj kin me Length

    meta         = all_meta(sessnum, 1);
    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);
    sessix = 1;

% ---- trial cap (cfg.exclusions) ----
    capN     = trialCapFor(obj.pth.anm, obj.pth.dt, cfg.exclusions);
    condtrix = (1:obj.bp.Ntrials)';
    if ~isnan(capN)
        condtrix(condtrix > capN) = [];
        logf('  [cap] %s %s: keeping trials 1-%d of %d\n', ...
            obj.pth.anm, obj.pth.dt, capN, obj.bp.Ntrials);
    end

    kinix = find(strcmp(kin(sessix).featLeg, 'tongue_length'));
    assert(~isempty(kinix), 'tongue_length not found in kin.featLeg');
    Length = kin(sessix).dat(:, condtrix, kinix);

% ---- trial sets ----
% Length is indexed by condtrix, so its columns are 1..numel(condtrix).
% Trial lists must be expressed in the SAME indexing (ismember/find below).
    allTr = (1:obj.bp.Ntrials)';
    hitTr = allTr(obj.bp.hit == 1);
    P1 = intersect(hitTr, allTr(obj.bp.rewardedLick == 1));   % reward on lick 1
    P4 = intersect(hitTr, allTr(obj.bp.rewardedLick == 4));   % reward on lick 4
    P  = {P1, P4};

    for k = 1:numel(P)
        P{k} = find(ismember(condtrix, P{k}));   % -> columns of Length
    end

% ---- per-condition extraction ----
    for k = 1:numel(P)
        trialList = P{k}(:)';
        dur_all = [];  len_all = [];  ili_all = [];

        for tr = trialList
            trAbs = condtrix(tr);   % absolute trial number, for obj.bp.ev
            [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial( ...
                Length(:,tr), obj.time, obj.bp.ev.lickL{trAbs,1}, obj.bp.ev.goCue(trAbs), ...
                params.dt, cfg);
            if why > 0
                if why < 5, continue; end   % why == 5 is a collision warning only
            end
            dur_all = [dur_all; durRow];   %#ok<AGROW>
            len_all = [len_all; lenRow];   %#ok<AGROW>
            ili_all = [ili_all; iliRow];   %#ok<AGROW>
        end

        nL = cfg.nLicksSummary;
        if isempty(dur_all)
            warning('session %d (%s %s) condition %d produced no usable trials.', ...
                sessnum, obj.pth.anm, obj.pth.dt, k);
            dur_all = nan(1,nL); len_all = nan(1,nL); ili_all = nan(1,nL);
        end
        dur_all = padTo(dur_all, nL);
        len_all = padTo(len_all, nL);
        ili_all = padTo(ili_all, nL);

% ONE VALUE PER SESSION PER LICK: mean (and SD) across trials. Sessions, not
% trials, are the unit of analysis.
        valDur   = mean(dur_all(:,1:nL), 1, 'omitnan')';
        valDurSD = std( dur_all(:,1:nL), 0, 1, 'omitnan')';
        valLen   = mean(len_all(:,1:nL), 1, 'omitnan')';
        valILI   = mean(ili_all(:,1:nL), 1, 'omitnan')';

        if k == 1
            allDur{sessnum}   = valDur;
            allDurSD{sessnum} = valDurSD;
            allLen{sessnum}   = valLen;
            allILI{sessnum}   = valILI;
        else
            allDur4{sessnum}   = valDur;
            allDurSD4{sessnum} = valDurSD;
            allLen4{sessnum}   = valLen;
            allILI4{sessnum}   = valILI;
        end
    end

    allAnm{sessnum}  = obj.pth.anm;
    allDate{sessnum} = obj.pth.dt;
    fprintf('Session %d\n', sessnum);
end

anmList = unique(allAnm(~cellfun(@isempty, allAnm)));

%% FIGURE
% Fig. 2B: reward-on-lick-1 vs reward-on-lick-4 trials at each lick; no test.
for pc = 1:numel(cfg.plotChoices)
    plotChoice = cfg.plotChoices{pc};

    switch lower(plotChoice)
        case 'mean', m1 = cell2mat(allDur)   * cfg.durScale;  m4 = cell2mat(allDur4)   * cfg.durScale;
                     yLab = sprintf('Mean protrusion duration (%s)', cfg.durUnit);
        case 'sd',   m1 = cell2mat(allDurSD) * cfg.durScale;  m4 = cell2mat(allDurSD4) * cfg.durScale;
                     yLab = sprintf('SD of protrusion duration (%s)', cfg.durUnit);
        otherwise,   error('plotChoice must be ''mean'' or ''sd''.');
    end
    m1 = m1(1:cfg.nLicksPlot, :);
    m4 = m4(1:cfg.nLicksPlot, :);
    nSessUsed = size(m1, 2);

    n1  = sum(isfinite(m1),2);   n4 = sum(isfinite(m4),2);
    mu1 = mean(m1,2,'omitnan');  mu4 = mean(m4,2,'omitnan');
    c1  = std(m1,0,2,'omitnan')./sqrt(n1) .* tinv(1-cfg.alpha/2, max(n1-1,1));
    c4  = std(m4,0,2,'omitnan')./sqrt(n4) .* tinv(1-cfg.alpha/2, max(n4-1,1));

    figure('Color','w','Position',[75 + 380*(pc-1), 75, 360, 620]); hold on
    h1 = []; h4 = [];
    for i = 1:cfg.nLicksPlot
        fill([i-0.22 i+0.22 i+0.22 i-0.22], [mu1(i)-c1(i) mu1(i)-c1(i) mu1(i)+c1(i) mu1(i)+c1(i)], ...
            [0.05 0.20 0.60], 'FaceAlpha',0.30, 'EdgeColor','none');   % R1 dark blue
        h1 = plot(i, mu1(i), 'o', 'Color',[0.05 0.20 0.60], ...
            'MarkerFaceColor',[0.05 0.20 0.60], 'MarkerSize',8);
        fill([i-0.22 i+0.22 i+0.22 i-0.22], [mu4(i)-c4(i) mu4(i)-c4(i) mu4(i)+c4(i) mu4(i)+c4(i)], ...
            [0.35 0.70 0.90], 'FaceAlpha',0.30, 'EdgeColor','none');   % R4 light blue
        h4 = plot(i, mu4(i), 'o', 'Color',[0.35 0.70 0.90], ...
            'MarkerFaceColor',[0.35 0.70 0.90], 'MarkerSize',8);
    end
    xlabel('Lick number'); ylabel(yLab);
    xlim([0.5 cfg.nLicksPlot+0.5]); xticks(1:cfg.nLicksPlot);
    legend([h1 h4], {'reward on lick 1', 'reward on lick 4'}, 'Location','best', 'Box','off');
    title(sprintf('Fig. 2B  (n = %d sessions)', nSessUsed), 'FontSize',10, 'FontWeight','normal');
    set(gca,'TickDir','out'); box off

end

%% LOCAL FUNCTIONS

function [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial(lenTrace, tAxis, licksAbs, gcTime, dt, cfg)
% Per-trial protrusion features, indexed BY LICK NUMBER so that lick identity is
% preserved across trials: column 3 is always the third post-cue lick, in every
% trial and every session.
% why: 0 ok | 1 fewer than 3 licks | 2 tongue all NaN | 3 no post-cue licks |
%      4 no run matched a lick | 5 several runs matched one lick (closest kept)
% nDrop: [runs dropped for being too short, runs dropped for being too long]
    durRow = nan(1, cfg.maxLickIdx);
    lenRow = nan(1, cfg.maxLickIdx);
    iliRow = nan(1, cfg.maxLickIdx);
    why = 0;
    nDrop = [0 0];

    if numel(licksAbs) <= 2,  why = 1; return; end
    if all(isnan(lenTrace)),  why = 2; return; end

    licks = licksAbs(:) - gcTime;
    licks = licks(licks > 0);
    if isempty(licks),        why = 3; return; end

% ---- bridge short tracking gaps ----
    bridgeSamps = round(cfg.nanBridge_sec / dt);
    val = bridgeNaNGaps(lenTrace(:), bridgeSamps);

% ---- lick events -> sample indices on the same (go-cue-aligned) axis ----
    lickIdx = zeros(numel(licks), 1);
    for q = 1:numel(licks)
        [~, lickIdx(q)] = min(abs(tAxis(:) - licks(q)));
    end

    minSamps = round(cfg.minBout_sec / dt);
    maxSamps = round(cfg.maxBout_sec / dt);
    tolSamps = round(cfg.contactTol_sec / dt);

% ---- contact runs, split only where the PORT CONTACTS say to ----
    runs = contactRuns(find(~isnan(val)), lickIdx, tolSamps);

% ---- length filter, applied AFTER splitting ----
    keepLen = false(1, numel(runs));
    for j = 1:numel(runs)
        n = numel(runs{j});
        keepLen(j) = n >= minSamps && n <= maxSamps;
        if n < minSamps, nDrop(1) = nDrop(1) + 1; end
        if n > maxSamps, nDrop(2) = nDrop(2) + 1; end
    end
    runs = runs(keepLen);
    if isempty(runs),         why = 4; return; end

% ---- match each surviving run to the nearest lick event ----
    nR = numel(runs);
    matchLick = nan(1, nR);
    matchDist = inf(1, nR);
    for j = 1:nR
        dists         = abs(runs{j}(:) - lickIdx(:)');
        minPerLick    = min(dists, [], 1);
        [dBest, best] = min(minPerLick);
        if dBest <= tolSamps
            matchLick(j) = best;
            matchDist(j) = dBest;
        end
    end
    keep = ~isnan(matchLick);
    runs = runs(keep);  matchLick = matchLick(keep);  matchDist = matchDist(keep);
    if isempty(runs),         why = 4; return; end

% ---- one bout per lick: closest wins ----
    [uL, ~, grp] = unique(matchLick);
    if numel(uL) < numel(matchLick), why = 5; end
    keep2 = false(1, numel(runs));
    for g = 1:numel(uL)
        ix     = find(grp == g);
        [~, b] = min(matchDist(ix));
        keep2(ix(b)) = true;
    end
    runs = runs(keep2);  matchLick = matchLick(keep2);

% ---- write features into lick-indexed columns ----
    firstSample = nan(1, cfg.maxLickIdx);
    for j = 1:numel(runs)
        c = matchLick(j);
        if c <= cfg.maxLickIdx
            durRow(c)      = numel(runs{j}) * dt;
            lenRow(c)      = max(val(runs{j}));
            firstSample(c) = runs{j}(1);
        end
    end

    v = find(~isnan(firstSample));
    if numel(v) > 1
        iliRow(v(1:end-1)) = diff(firstSample(v)) * dt;
    end
end

function runs = contactRuns(idx, lickIdx, tolSamps)
% Runs of consecutive tongue-visible samples, split only where the port-contact
% times show more than one lick inside a run. If two or more contacts fall
% inside a run (within tolSamps), it is split at the midpoints between
% consecutive contacts. With one contact the run is kept whole, and the length
% filter in the caller removes runs that are still too long.
    runs = {};
    if isempty(idx), return; end
    idx = idx(:)';
    brk = find(diff(idx) > 1);
    s   = idx([1, brk+1]);
    e   = idx([brk, numel(idx)]);

    lickIdx = sort(lickIdx(:))';
    for r = 1:numel(s)
        inside = lickIdx(lickIdx >= s(r) - tolSamps & lickIdx <= e(r) + tolSamps);
        if numel(inside) <= 1
            runs{end+1} = s(r):e(r);   %#ok<AGROW>
            continue
        end
        cuts = floor((inside(1:end-1) + inside(2:end)) / 2);
        prev = s(r);
        for c = 1:numel(cuts)
            hi = min(cuts(c), e(r));
            if hi >= prev, runs{end+1} = prev:hi; end   %#ok<AGROW>
            prev = hi + 1;
        end
        if e(r) >= prev, runs{end+1} = prev:e(r); end   %#ok<AGROW>
    end
end

function val = bridgeNaNGaps(val, maxGap)
% Linearly interpolate NaN gaps of at most maxGap samples. Gaps at the very
% start or end of the trace are left alone -- there is nothing to interpolate
% between.
    nanIdx = find(isnan(val));
    if isempty(nanIdx), return; end
    brk = find(diff(nanIdx) > 1);
    gs  = [nanIdx(1); nanIdx(brk+1)];
    ge  = [nanIdx(brk); nanIdx(end)];
    for g = 1:numel(gs)
        n  = ge(g) - gs(g) + 1;
        i0 = gs(g) - 1;
        i1 = ge(g) + 1;
        if n <= maxGap && i0 >= 1 && i1 <= numel(val)
            val(gs(g):ge(g)) = linspace(val(i0), val(i1), n);
        end
    end
end

function M = padTo(M, n)
    if size(M,2) < n, M(:, end+1:n) = NaN; end
end

function capN = trialCapFor(anm, dte, caps)
    capN = NaN;
    for r = 1:size(caps,1)
        if strcmpi(caps{r,1}, anm) && strcmpi(caps{r,2}, dte)
            capN = caps{r,3};  return
        end
    end
end

function [p, pAdj, sig] = bhCorrect(p, alpha)
% Benjamini-Hochberg adjusted p-values; m = numel(p), the number of licks in
% the family.
    m = numel(p);
    [ps, ix] = sort(p);
    pAdj = nan(1, m);
    for i = 1:m
        pAdj(ix(i)) = min(ps(i:end) .* m ./ (i:m));
    end
    pAdj = min(pAdj, 1);
    sig  = pAdj < alpha;
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

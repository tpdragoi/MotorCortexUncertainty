%% F4C_lickDuration.m
%  Duration of successive tongue protrusions, VTA Reward Task.
%  Protrusions are detected from tongue length, summarized per session, then
%  averaged across sessions.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/200
%    params.smooth      1
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2.5 to 4 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


% ---- SELF-CONTAINED: data and functions come only from this folder ----
% Data\VTA holds the exported <ANM>_<DATE>_obj.mat / _kin.mat; shared\ holds the
% pipeline functions it needs (shared\pipelineCopies). uninstructedMovements_v2-main
% and the raw data tree are not on the path.
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'loadBehavSession.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));
dataDir = fullfile(v2Root, 'Data', 'VTA');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s', dataDir);

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
% BOUT DETECTION -- shared with the R14 / R16 scripts.
cfg.nanBridge_sec   = 0.020;
cfg.minBout_sec     = 0.030;
cfg.maxBout_sec     = 0.400;
cfg.contactTol_sec  = 0.020;
cfg.maxLickIdx      = 35;
cfg.nLicksSummary   = 15;

% TRIAL SELECTION -- specific to this panel.
% A trial is used only if at least cfg.minContacts port contacts fall within
% cfg.contactWin_s AFTER the go cue. Contacts before the cue are ignored; the
% previous version counted them, which is the bug in note A.
cfg.minContacts  = 4;
cfg.contactWin_s = 1.5;
cfg.hitOnly      = false;   % true restricts to hit == 1

% ANALYSIS / FIGURE
cfg.nLicksPlot  = 6;
cfg.alpha       = 0.05;
cfg.plotChoices = {'mean','sd'};   % both measures are plotted, one figure each
cfg.tail        = 'both';   % the old script ran a two-tailed test here
cfg.durScale    = 1000;
cfg.durUnit     = 'ms';

cfg.exclusions = {
    'TDv1', '2025-02-15', 211;
};

%% SESSIONS TO LOAD
% loader | date | animal label | use
sessionTable = {
    'loadTDv1_many', '2025-02-15', 'TDv1', true
    'loadTDv1_many', '2025-02-17', 'TDv1', true
    'loadTDv1_many', '2025-02-18', 'TDv1', true
    'loadTDv1_many', '2025-02-19', 'TDv1', true

    'loadTDv4_many', '2025-02-25', 'TDv4', true
    'loadTDv4_many', '2025-02-26', 'TDv4', true
    'loadTDv4_many', '2025-02-27', 'TDv4', true
    'loadTDv4_many', '2025-02-28', 'TDv4', true

    'loadTDv6_many', '2025-08-21', 'TDv6', true
    'loadTDv6_many', '2025-08-22', 'TDv6', true
    'loadTDv6_many', '2025-08-23', 'TDv6', true
    'loadTDv6_many', '2025-08-24', 'TDv6', true
    'loadTDv6_many', '2025-08-25', 'TDv6', true
    'loadTDv6_many', '2025-08-26', 'TDv6', true

    'loadTDv5_many', '2025-08-25', 'TDv5', true
    'loadTDv5_many', '2025-08-26', 'TDv5', true
    'loadTDv5_many', '2025-08-28', 'TDv5', false
    'loadTDv5_many', '2025-08-30', 'TDv5', true
    'loadTDv5_many', '2025-08-31', 'TDv5', true
};

% ---- build all_meta from the table; the metadata travels with the row ----
all_meta = [];
sessAnm  = {};
sessDate = {};
sessTag  = [];
for r = 1:size(sessionTable,1)
    if ~sessionTable{r,4}, continue; end
    m        = struct('anm', sessionTable{r,3}, 'date', sessionTable{r,2});   % was feval(loader)
    all_meta = [all_meta; m];   %#ok<AGROW>
    sessAnm{end+1}  = sessionTable{r,3};   %#ok<SAGROW>
    sessDate{end+1} = sessionTable{r,2};   %#ok<SAGROW>
    sessTag(end+1)  = 0;   %#ok<SAGROW>
end
nSess = size(all_meta,1);
assert(nSess == numel(sessAnm), 'session labels and all_meta are out of step');

%% MAIN LOOP
allDur   = cell(1, nSess);   allDurSD = cell(1, nSess);
allLen   = cell(1, nSess);   allILI   = cell(1, nSess);
allAnm   = cell(1, nSess);   allDate  = cell(1, nSess);


for sessnum = 1:nSess

    clear obj kin me Length

    meta         = all_meta(sessnum, 1);
    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);
    sessix = 1;

% ---- trial cap ----
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

% ---- TRIAL SELECTION: at least cfg.minContacts contacts within
%      cfg.contactWin_s AFTER the go cue (see note A) ----
    allTr  = (1:obj.bp.Ntrials)';
    nEarly = zeros(obj.bp.Ntrials, 1);
    for t = 1:obj.bp.Ntrials
        lks = obj.bp.ev.lickL{t};
        if isempty(lks), continue; end
        rel = lks(:) - obj.bp.ev.goCue(t);
        nEarly(t) = sum(rel > 0 & rel <= cfg.contactWin_s);
    end
    P1 = allTr(nEarly >= cfg.minContacts);
    if cfg.hitOnly
        P1 = intersect(P1, allTr(obj.bp.hit == 1));
    end
    P = {P1};

% Length is indexed by condtrix, so its columns are 1..numel(condtrix).
% Trial lists must be expressed in the SAME indexing.
    for k = 1:numel(P)
        P{k} = find(ismember(condtrix, P{k}));
    end

    for k = 1:numel(P)
        trialList = P{k}(:)';
        dur_all = [];  len_all = [];  ili_all = [];

        for tr = trialList
            trAbs = condtrix(tr);
            [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial( ...
                Length(:,tr), obj.time, obj.bp.ev.lickL{trAbs,1}, obj.bp.ev.goCue(trAbs), ...
                params.dt, cfg);
            if why > 0 && why < 5, continue; end
            dur_all = [dur_all; durRow];   %#ok<AGROW>
            len_all = [len_all; lenRow];   %#ok<AGROW>
            ili_all = [ili_all; iliRow];   %#ok<AGROW>
        end

        nL = cfg.nLicksSummary;
        if isempty(dur_all)
            warning('session %d (%s %s) produced no usable trials.', ...
                sessnum, obj.pth.anm, obj.pth.dt);
            dur_all = nan(1,nL); len_all = nan(1,nL); ili_all = nan(1,nL);
        end
        dur_all = padTo(dur_all, nL);
        len_all = padTo(len_all, nL);
        ili_all = padTo(ili_all, nL);

% ONE VALUE PER SESSION PER LICK. The tests below are over sessions, so
% the session is the unit of analysis and individual trials must not
% enter the test as independent samples.
        allDur{sessnum}   = mean(dur_all(:,1:nL), 1, 'omitnan')';
        allDurSD{sessnum} = std( dur_all(:,1:nL), 0, 1, 'omitnan')';
        allLen{sessnum}   = mean(len_all(:,1:nL), 1, 'omitnan')';
        allILI{sessnum}   = mean(ili_all(:,1:nL), 1, 'omitnan')';
    end

    allAnm{sessnum}  = obj.pth.anm;
    allDate{sessnum} = obj.pth.dt;
    fprintf('Session %d\n', sessnum);
end

anmList = unique(allAnm(~cellfun(@isempty, allAnm)));

%% FIGURE + STATISTICS
% Each protrusion compared to the FIRST, Wilcoxon signed-rank paired within
% set from the hypothesis BEFORE looking at the data.
for pc = 1:numel(cfg.plotChoices)
    plotChoice = cfg.plotChoices{pc};

    switch lower(plotChoice)
        case 'mean', mat = cell2mat(allDur)   * cfg.durScale;
                     yLab = sprintf('Mean protrusion duration (%s)', cfg.durUnit);
        case 'sd',   mat = cell2mat(allDurSD) * cfg.durScale;
                     yLab = sprintf('SD of protrusion duration (%s)', cfg.durUnit);
        otherwise,   error('plotChoice must be ''mean'' or ''sd''.');
    end
    mat = mat(1:cfg.nLicksPlot, :);
    nSessUsed = size(mat, 2);

    nOK = sum(isfinite(mat), 2);
    mu  = mean(mat, 2, 'omitnan');
    sem = std(mat, 0, 2, 'omitnan') ./ sqrt(nOK);
    ci  = sem .* tinv(1 - cfg.alpha/2, max(nOK - 1, 1));

    switch lower(cfg.tail)
        case 'both',  tailA = {};                 tailWord = 'two-tailed';
        case 'left',  tailA = {'tail','left'};    tailWord = 'one-tailed (lick i < lick 1)';
        case 'right', tailA = {'tail','right'};   tailWord = 'one-tailed (lick i > lick 1)';
        otherwise,    error('cfg.tail must be ''both'', ''left'' or ''right''.');
    end

    licksTested = 2:cfg.nLicksPlot;
    p     = nan(1, numel(licksTested));
    p2    = nan(1, numel(licksTested));
    nPair = zeros(1, numel(licksTested));
    for ii = 1:numel(licksTested)
        a  = mat(licksTested(ii), :);
        b  = mat(1, :);
        ok = isfinite(a) & isfinite(b);
        nPair(ii) = sum(ok);
        if nPair(ii) < 2, continue; end
        p(ii)  = signrank(a(ok), b(ok), tailA{:});
        p2(ii) = signrank(a(ok), b(ok));
    end
    [p, pAdj, sig] = bhCorrect(p, cfg.alpha);

    fprintf('\n%s\n', repmat('=',1,78));
    fprintf('VTA | Wilcoxon signed-rank, each protrusion vs the first (%s)\n', plotChoice);
    fprintf('%s | m = %d comparisons | n = %d sessions\n', tailWord, numel(licksTested), nSessUsed);
    fprintf('trials: >= %d contacts within %.2f s of the go cue', cfg.minContacts, cfg.contactWin_s);
    if cfg.hitOnly, fprintf(', hit trials only\n'); else, fprintf(', all outcomes\n'); end
    fprintf('%s\n', repmat('-',1,78));
    fprintf('%-6s %-8s %-10s %-11s %-11s %-6s %s\n', ...
        'Lick', 'n sess', 'mean', 'p (raw)', 'p (BH)', 'sig', 'p two-tailed');
    for ii = 1:numel(licksTested)
        fprintf('%-6d %-8d %-10.2f %-11.4f %-11.4f %-6s %.4f\n', licksTested(ii), nPair(ii), ...
            mu(licksTested(ii)), p(ii), pAdj(ii), string(sig(ii)), p2(ii));
    end
    fprintf('%s\n', repmat('=',1,78));

    figure('Color','w','Position',[75 + 360*(pc-1), 75, 340, 620]); hold on
    for i = 1:cfg.nLicksPlot
        fill([i-0.22 i+0.22 i+0.22 i-0.22], ...
            [mu(i)-ci(i) mu(i)-ci(i) mu(i)+ci(i) mu(i)+ci(i)], ...
            [0.10 0.60 0.25], 'FaceAlpha',0.30, 'EdgeColor','none');   % VTA green
        plot(i, mu(i), 'o', 'Color',[0.10 0.60 0.25], ...
            'MarkerFaceColor',[0.10 0.60 0.25], 'MarkerSize',8);
    end
    yl = ylim;
    for ii = 1:numel(licksTested)
        if sig(ii)
            text(licksTested(ii), yl(2) - 0.03*range(yl), '*', ...
                'HorizontalAlignment','center', 'FontSize',16, 'FontWeight','bold');
        end
    end
    xlabel('Lick number'); ylabel(yLab);
    xlim([0.5 cfg.nLicksPlot+0.5]); xticks(1:cfg.nLicksPlot);
    title(sprintf('VTA  (n = %d sessions, %d animals)', nSessUsed, numel(anmList)), ...
        'FontSize',10, 'FontWeight','normal');
    set(gca,'TickDir','out'); box off

end

%% LOCAL FUNCTIONS

function [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial(lenTrace, tAxis, licksAbs, gcTime, dt, cfg)
% Per-trial protrusion features, indexed BY LICK NUMBER so that lick identity is
% preserved across trials: column 3 is always the third post-cue lick, in every
% trial and every session.
% why: 0 ok | 1 fewer than 3 licks | 2 tongue all NaN | 3 no post-cue licks
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
% The old code let a later bout overwrite an earlier one silently, so which
% protrusion a lick ended up with depended on loop order.
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
% Runs of consecutive tongue-visible samples, split ONLY where the port-contact
% times say there is more than one lick inside a run.
% THIS IS THE FIX THAT MATTERS. The previous findConsecutiveSets did:
% became a 200 ms lick plus a 50 ms lick, and every lick after it in that trial
% and it pulled mean duration at low lick numbers DOWNWARD.
% Simply DISCARDING long runs is no better: it throws away the genuinely long
% first protrusions and biases lick 1 downward as well.
% What the old code never used is the information that decides the question.
% A run is too long for one of two reasons: two protrusions ran together
% because tracking never dropped out between them, or it is an artefact. The
% port-contact times distinguish them. If two or more contacts fall inside a
% run, it is a merged pair and it is split at the MIDPOINT between consecutive
% contacts. If only one contact falls inside it, it is one lick and it is kept
% whole however long it is -- the length filter in the caller then removes what
% is left as an artefact.
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
% Benjamini-Hochberg adjusted p-values. m is numel(p) -- the number of licks in
% this family -- which is the number the Methods should quote.
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

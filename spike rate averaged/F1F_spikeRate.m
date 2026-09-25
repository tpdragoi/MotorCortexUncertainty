%% F1F_spikeRate.m
%  Session-normalized mean spike rate across units, Simple Reward Task.
%  Rates are averaged within a session and normalized, then averaged across
%  sessions and plotted against mean tongue length.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'firstLick'
%    params.dt          1/200
%    params.smooth      30
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2.5 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc;

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.

rewardLickB = 4;
%% PATHS

% ---- SELF-CONTAINED: data and functions come only from this folder ----
% Sessions are read from the exported obj/kin files in <MATLAB Codes _ v2>\Data
% through slimMeta / slimToLegacy / loadSlimSession (in <MATLAB Codes _ v2>\shared).
% slimToLegacy rebuilds obj/params/kin for THIS script's params (dt, window,
% alignment, smooth via the pipeline's mySmooth, quality via findClusters, lowFR,
% conditions via findTrials) from the exported spike times, as the pipeline did.
% Nothing here reads uninstructedMovements_v2-main or the raw data tree: the
% pipeline functions still needed are byte-identical copies in shared\pipelineCopies.
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'slimToLegacy.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));
spec.dataDir = '';   % raw data folder not used (was: spec.dataDir = fullfile(projectRoot,'data');)
%% SETTINGS (identical in all five scripts)
% ONE parameter set across R1, R14, R16, VTA and Learning, matching the
% tongue/jaw DECODING scripts, so a firing-rate figure and a decoding figure
% from the same session describe the same trials, units and clock.
params.alignEvent = 'firstLick';   % t = 0 is the FIRST LICK CONTACT
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 8;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz
params.quality    = {'good'};   % good units only (findClusters trims blanks and ignores case)
params.tmin   = -2.5;
params.tmax   = 5;
params.dt     = 1/200;   % 5 ms bins
params.smooth = 30;   % causal gaussian, in bins (neural AND tongue)
params.traj_features = { ...
    {'tongue','left_tongue','right_tongue','jaw','trident','nose'}, ...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'} };
params.feat_varToExplain = 80;
params.N_varToExplain    = 80;
params.advance_movement  = 0;
params.fcut   = 10;
params.cond   = 5;
params.method = 'xcorr';
params.fa     = false;
params.bctype = 'reflect';

params.condition(1)     = {'hit==1 | hit==0'};   % 1
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 1'};   % 2
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 1'};   % 3
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 1'};   % 4
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 1 & rewardedLick == %d', rewardLickB)};   % 5
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 2 & rewardedLick == %d', rewardLickB)};   % 6
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 3 & rewardedLick == %d', rewardLickB)};   % 7
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};   % 8
params.condition(end+1) = {sprintf('hit==1 & rewardedLick == %d', rewardLickB)};   % 9
params.condition(end+1) = {'hit==1'};   % 10

%% THE TRACE
% One trace per (session, probe, condition set): mean firing rate over UNITS,
% then over TRIALS, in Hz. One number per session enters every average, so a
% 200-unit session and an 8-unit session count equally -- deliberate, matching
% the old scripts, and the reason unit counts are printed per region.
cfg.minTrials = 5;   % a (session, probe, condition) with fewer trials is skipped

%% NORMALISATION: TRIM, SCALE, THEN AVERAGE
% This is the recipe, in the order it runs:
%   1. TRIM the session's mean trace to [go cue, contact cfg.normLastContact].
%      Both edges are medians across that condition's trials, on the trace clock
%      where t = 0 is contact 1, so the go cue edge is negative.
%   2. MIN-MAX that window to [0, 1]. The window's minimum is the quiet gap
%      between the go cue and contact 1; its maximum is the early-lick peak. So
%      0 means "back to pre-lick rest" and 1 means "peak of the bout", both
%      measured on this animal, in this session.
%   3. AVERAGE the normalised trace inside each epoch.
% WHY THE TRIM MATTERS. Scaling over the whole -2.5 to 5 s trace, or over the
% epoch values themselves, both give a much narrower axis: a session whose firing
% falls from 18 to 12 Hz across the bout reads 1.00 / 0.87 / 0.71 when each epoch
% is divided by epoch 1, but 0.73 / 0.59 / 0.42 under this recipe. The trimmed
% window is what puts a real floor under the scale.
cfg.normLastContact = 8;   % the trim ends at this contact
cfg.normMode        = 'window';   % 'window' (above) | 'trace' (whole trace) | 'none' (Hz)

% Divide the epoch vector by its FIRST epoch afterwards. true pins every curve to
% exactly 1.0 at epoch 1 and shows the fraction of the early-bout response
% remaining -- the published panel. false leaves the 0-1 window scale intact, so
% the height of epoch 1 is itself readable and comparable across groups.
cfg.epochDivideByFirst = true;

% The time-course figure (FIGURE 1) is scaled separately, over the same trimmed
% window, so the two figures agree about what 0 and 1 mean.
cfg.traceNorm        = 'window';   % 'window' | 'trace' | 'none'
cfg.traceNorm_tongue = true;   % rescale the tongue trace to [0,1] for display

%% DISPLAY SMOOTHING
% The smoothing of main_spikeRate_manySessions_R1_forScience.m (its
% first figure). After the per-session traces are averaged ACROSS SESSIONS, the
% mean spike-rate trace and the mean tongue trace are each passed through a
% moving average whose width shrinks as more sessions go in:
% (17 = the session count those widths were set on). The 95% CI band is NOT
% smoothed, as in that script. The widths are 5 ms samples; at another
% params.dt 1/200
% Per-session steps of that script were already identical here: causal
% gaussian params.smooth = 20, tongue NaN -> 0, tongue min-max per session,
% t-based 95% CI across sessions. Only the per-session RATE scaling differs
% (that script min-maxed the whole trace; this one uses cfg.traceNorm) and is
cfg.dispSmooth           = true;
cfg.dispSmoothBaseN      = 17;
cfg.dispSmoothRateBase   = 6;
cfg.dispSmoothTongueBase = 22;
cfg.dispSmoothRefDt      = 0.005;

%% EPOCHS: THE 'PORT CONTACTS' FIGURE
% Binned on REAL CONTACT TIMES (medians across that condition's trials), not on
% hardcoded sample indices. The old scripts used fixed windows under an axis
% labelled by lick number; lick rate changes across animals and across learning
% days, so a window covering contacts 1-3 on day 5 can cover only 1-2 on day 1.
cfg.epochContacts = [1 3; 3 5; 5 8];   % one [firstContact lastContact] row per epoch
cfg.epochLabels   = {};   % {} = built automatically ('1-3', '3-5', ...)

%% STATISTICS
% WILCOXON SIGNED RANK throughout, paired, and that is what the figure says.
%   within group : every epoch against EPOCH 1, paired across sessions
%   between      : spec.compareGroups at every epoch, paired BY ANIMAL
cfg.alpha  = 0.05;
cfg.useFDR = true;   % Benjamini-Hochberg across the epochs of one comparison

%% FIGURE / RUN CONTROL
cfg.showTrace   = true;
cfg.traceXlim_s = [-0.5 2.5];
cfg.rngSeed     = [];
cfg.sessionsToRun = [];   % [] = all; otherwise a list of session indices
cfg.useCache    = true;
cfg.cacheDir    = fullfile(tempdir, 'spikeRate_cache_clean');
%% STUDY: R1 (reward on lick 1), M1 + ALM
% Session list, loaders, group maps, trial caps and single-probe sessions are
% taken VERBATIM from tongue_r1.m / jaw_r1.m, so this figure and the decoding
% figures describe the same 26 sessions and the same probes.
spec.name = 'R1';
spec.sessionDates = { ...
    '2024-07-09','2024-07-10','2024-07-11','2024-07-13','2024-07-14', ...   % TD10si
    '2024-07-05','2024-07-06','2024-07-07','2024-07-08','2024-07-09','2024-07-10', ...   % TD9si
    '2025-07-25','2025-07-26','2025-07-27','2025-07-28','2025-07-29','2025-07-30','2025-07-31', ...   % TD27d
    '2025-07-30','2025-07-31','2025-08-01','2025-08-02','2025-08-03','2025-08-05','2025-08-06','2025-08-07' };   % TD26d
spec.sessionLoaders = { ...
    @loadTD10s_neur, @loadTD10s_neur, @loadTD10s_neur, @loadTD10s_neur713, @loadTD10s_neur, ...
    @loadTD9s_neur,  @loadTD9s_neur,  @loadTD9s_neur,  @loadTD9s_neur,  @loadTD9s_neur709, @loadTD9s_neur, ...
    @loadTD27_neur,  @loadTD27_neur,  @loadTD27_neur,  @loadTD27_neur,  @loadTD27_neur, @loadTD27_neur, @loadTD27_neur, ...
    @loadTD26_neur,  @loadTD26_neur,  @loadTD26_neur,  @loadTD26_neur222, @loadTD26_neur222, ...
    @loadTD26_neur,  @loadTD26_neur,  @loadTD26_neur };
spec.groupMaps = { ...
    [0 1 0 2 2  2 0 2 2 1 0  2 1 2 0 2 0 1  1 0 2 0 1 2 2 1] , ...   % M1
    [2 2 2 0 1  1 0 0 0 0 1  0 0 1 1 1 0 0  2 0 0 0 0 1 1 0] };   % ALM
spec.groupLabels  = {'M1','ALM'};
spec.groupColours = [0.85 0.10 0.10 ; 0.10 0.35 0.75];
spec.trialCaps           = { 'TD26d','2025-08-07', 187 };
spec.singleProbeSessions = { 'TD10si','2024-07-13' ; 'TD9si','2024-07-09' ; ...
                             'TD26d','2025-08-02' ; 'TD26d','2025-08-03' };
spec.condSets = struct('name', {'R1'}, 'condIdx', {8}, ...
                       'colour', {[0.85 0.10 0.10]});
spec.showEpochFigure = false;   % the 3-point panel is Learning-only
spec.compareGroups   = [];
%% MAIN LOOP (shared)
nSess = numel(spec.sessionDates);
if ~isempty(cfg.rngSeed), rng(cfg.rngSeed); end

if isempty(cfg.sessionsToRun)
    runList = 1:nSess;
else
    runList = unique(cfg.sessionsToRun(cfg.sessionsToRun >= 1 & cfg.sessionsToRun <= nSess));
end
assert(~isempty(runList), 'cfg.sessionsToRun left no valid sessions.');

% ONE group definition per script, asserted against the session count. In the old
% main_spikeRate_* files m1/alm were re-declared in several cells with different
% contents and sometimes the wrong length -- R16 carried the VTA script's
% 19-element vectors for its 33 sessions. This is what makes that impossible.
for d = 1:numel(spec.groupMaps)
    assert(numel(spec.groupMaps{d}) == nSess, ...
        'group "%s" map has %d entries but there are %d sessions.', ...
        spec.groupLabels{d}, numel(spec.groupMaps{d}), nSess);
end
runMask = false(1, nSess);  runMask(runList) = true;
for d = 1:numel(spec.groupMaps)
    spec.groupMaps{d}(~runMask) = 0;
end
probesOf = cell(1, nSess);
for s = runList
    pr = cellfun(@(m) m(s), spec.groupMaps);
    probesOf{s} = unique(pr(pr > 0));
end

nG   = numel(spec.groupMaps);
nCS  = numel(spec.condSets);
epRows   = cfg.epochContacts;
nEp      = size(epRows,1);
epLabels = epochLabels(cfg, epRows);

fprintf('\n=== %s | %d of %d sessions | %d (session,probe) pairs | %d condition set(s) ===\n', ...
    spec.name, numel(runList), nSess, sum(cellfun(@numel, probesOf)), nCS);
for c = 1:nCS
end

res = struct('session',{},'probe',{},'condSet',{},'anm',{},'date',{},'region',{}, ...
             'nUnits',{},'nTrials',{},'time',{},'rateHz',{},'rateNorm',{}, ...
             'tongue',{},'epMean',{},'epN',{},'winEdges',{});

% Trials per condition, per session. Kept so that a condition selecting nothing
% reports WHICH condition and what the alternatives hold, instead of failing at
% the end with an empty result and no way to tell why.
census    = nan(nSess, numel(params.condition));
censusHdr = false;

sessCount = 0;
for sessionIdx = runList
    sessCount = sessCount + 1;
    fprintf('Session %d\n', sessCount);
if isempty(probesOf{sessionIdx}), continue; end

S = srLoadSession(spec.sessionLoaders{sessionIdx}, spec.sessionDates{sessionIdx}, spec, params, cfg);
anm = S.anm;  dte = S.date;
traceTime = S.time;
nT = numel(traceTime);
dtS = median(diff(traceTime));

assert(size(S.trialdat,1) == nT, 'sess %d: trialdat %d samples, time %d.', sessionIdx, size(S.trialdat,1), nT);
assert(size(S.tongueRaw,1) == nT, 'sess %d: kinematics %d samples, time %d.', sessionIdx, size(S.tongueRaw,1), nT);

maxUsable = min([size(S.tongueRaw,2), size(S.trialdat,3), S.Ntrials]);
capRow = find(strcmp(anm, spec.trialCaps(:,1)) & strcmp(dte, spec.trialCaps(:,2)), 1);
if ~isempty(capRow), maxUsable = min(maxUsable, spec.trialCaps{capRow,3}); end

% ---------------- CONTACT AND GO-CUE TIMES ON THE TRACE CLOCK ------------
% alignEvent is 'firstLick', so t = 0 is the first post-go-cue contact: contact k
% sits at post(k) - post(1) and the go cue at goCue - post(1), which is negative.
contactRel = cell(maxUsable,1);
goCueRel   = nan(maxUsable,1);
for tr = 1:maxUsable
    lk = S.lickL{tr};
    if isempty(lk), continue; end
    post = sort(lk(lk > S.goCue(tr)));
    if isempty(post), continue; end
    contactRel{tr} = post(:) - post(1);
    goCueRel(tr)   = S.goCue(tr) - post(1);
end
hasContacts = ~cellfun(@isempty, contactRel);

% ---------------- CONDITION CENSUS ----------------
for ci = 1:numel(S.trialid)
    tr = S.trialid{ci};
    tr = tr(tr >= 1 & tr <= maxUsable);
    census(sessionIdx, ci) = sum(hasContacts(tr));
end
if ~censusHdr
    parts = cell(1, numel(S.trialid));
    for ci = 1:numel(S.trialid)
        parts{ci} = sprintf('c%d=%d', ci, census(sessionIdx, ci));
    end
    censusHdr = true;
end

% ---------------- UNITS PER PROBE (once per session) ---------------------
isSingle = any(strcmp(anm, spec.singleProbeSessions(:,1)) & strcmp(dte, spec.singleProbeSessions(:,2)));

for probeIdx = probesOf{sessionIdx}(:)'
regionStr = groupOfProbe(spec, sessionIdx, probeIdx);
if isSingle
    unitIds = 1:size(S.cluid,1);
elseif probeIdx == 1
    unitIds = 1:size(S.cluid{1,1},1);
else
    unitIds = size(S.cluid{1,1},1)+1 : S.nUnitsTotal;
end
if isempty(unitIds)
    logf('  [skip] sess %2d probe %d (%s): no units\n', sessionIdx, probeIdx, regionStr);
    continue
end

for cs = 1:nCS
    trials = unique(cell2mat(S.trialid(spec.condSets(cs).condIdx)'));
    trials = trials(trials <= maxUsable);
    trials = trials(hasContacts(trials));
    if numel(trials) < cfg.minTrials
        logf('  [skip] sess %2d probe %d %-4s: %d trials\n', ...
            sessionIdx, probeIdx, spec.condSets(cs).name, numel(trials));
        continue
    end

% ---- population trace, in Hz -----------------------------------------
% Kept in SINGLE through both reductions. The old code cast the whole
% nT x nUnits x nTrials block to double first, which for a 1500-sample,
% 120-unit, 260-trial session allocated ~370 MB per probe per condition and
    rateHz = double(mean(mean(S.trialdat(:, unitIds, trials), 2), 3));

% ---- tongue length ----------------------------------------------------
% NaN means the tongue is in the mouth, i.e. length 0 -- a real observation,
% zero-filled exactly as the tongue decoding scripts do.
    tg = S.tongueRaw(:, trials);
    tg(isnan(tg)) = 0;
    tg = mySmooth(tg, params.smooth, params.bctype);   % same causal kernel as the spikes, per trial
    tongue = mean(tg, 2);

% ---- the trim window and the epoch boundaries -------------------------
% Medians across this condition's trials, because the trace being scaled is
% itself a trial average.
    cMed  = medianContacts(contactRel(trials), max(cfg.normLastContact, max(epRows(:))));
    gcMed = median(goCueRel(trials), 'omitnan');
    if ~isfinite(gcMed) || ~isfinite(cMed(cfg.normLastContact))
        logf('  [skip] sess %2d probe %d %-4s: no contact %d\n', ...
            sessionIdx, probeIdx, spec.condSets(cs).name, cfg.normLastContact);
        continue
    end
    winEdges = [gcMed, cMed(cfg.normLastContact)];

    rateNorm = scaleTrace(rateHz, traceTime, winEdges, cfg.normMode);
    [epMean, epN] = epochFromTrace(rateNorm, traceTime, dtS, cMed, epRows);
    if cfg.epochDivideByFirst && isfinite(epMean(1)) && epMean(1) ~= 0
        epMean = epMean ./ epMean(1);
    end

    k = numel(res) + 1;
    res(k).session = sessionIdx;   res(k).probe   = probeIdx;   res(k).condSet = cs;
    res(k).anm     = anm;          res(k).date    = dte;        res(k).region  = regionStr;
    res(k).nUnits  = numel(unitIds);
    res(k).nTrials = numel(trials);
    res(k).time     = traceTime;
    res(k).rateHz   = rateHz;
    res(k).rateNorm = scaleTrace(rateHz, traceTime, winEdges, cfg.traceNorm);
    res(k).tongue   = tongue;
    res(k).epMean   = epMean;
    res(k).epN      = epN;
    res(k).winEdges = winEdges;

end
end   % probe
clear S
end   % session

if isempty(res)
    fprintf('\n%s\n  NO SESSIONS PRODUCED DATA -- trials per condition, per session\n%s\n', ...
        repmat('=',1,72), repmat('-',1,72));
    fprintf('  %-6s', 'sess');
    for k = 1:size(census,2), fprintf(' %6s', sprintf('c%d', k)); end
    fprintf('\n');
    for s = runList
        if all(isnan(census(s,:))), continue; end
        fprintf('  %-6d', s);  fprintf(' %6d', census(s,:));  fprintf('\n');
    end
    used = unique([spec.condSets.condIdx]);
    error(['Every condition set selected fewer than cfg.minTrials (%d) trials.\n' ...
           'This script asks for condition %s: "%s".\n' ...
           'If that column is zero above, spec.condSets is pointing at the wrong ' ...
           'condition -- it must match the matching decoding script''s ' ...
           'spec.poolCondIdx (8 for R1, [8 9] for R14/R16, 1 for VTA and Learning).'], ...
           cfg.minTrials, mat2str(used), params.condition{used(1)});
end

%% NEURON COUNTS PER REGION
% Totals over the (session, probe) pairs that actually entered the analysis --
% each pair counted ONCE, not once per condition set.
fprintf('\n%s\n  NEURONS ENTERING THE ANALYSIS\n%s\n', repmat('=',1,72), repmat('-',1,72));
fprintf('  %-10s %9s %9s %11s %11s %11s\n', 'region','sessions','neurons','per session','min','max');
grandPairs = 0;  grandUnits = 0;
for g = 1:nG
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    if isempty(rows)
        fprintf('  %-10s %9d %9s %11s %11s %11s\n', spec.groupLabels{g}, 0, '-', '-', '-', '-');
        continue
    end
% one row per (session, probe): collapse the condition sets
    key = arrayfun(@(r) r.session*100 + r.probe, res(rows));
    [~, ia] = unique(key, 'stable');
    u = [res(rows(ia)).nUnits];
    fprintf('  %-10s %9d %9d %11.1f %11d %11d\n', ...
        spec.groupLabels{g}, numel(ia), sum(u), mean(u), min(u), max(u));
    grandPairs = grandPairs + numel(ia);
    grandUnits = grandUnits + sum(u);
end
fprintf('%s\n  %-10s %9d %9d\n%s\n', repmat('-',1,72), 'TOTAL', grandPairs, grandUnits, repmat('=',1,72));
%% FIGURE 1: RATE vs TONGUE LENGTH
% One figure per GROUP. Where a task has two reward conditions (R1 and R4, R1
% and R6) both are drawn ON THE SAME AXES rather than in separate figures, so
% the comparison is made by eye in one place.
if cfg.showTrace
for g = 1:nG
    rowsG = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    if isempty(rowsG), continue; end

    fh = figure('Name', sprintf('%s - %s - population rate', spec.name, spec.groupLabels{g}), 'Color','w');
    ax = axes('Parent', fh);  hold(ax,'on')
    h = gobjects(0);  lbl = {};  nUsed = 0;

    for cs = 1:nCS
        rows = rowsG(arrayfun(@(r) r.condSet == cs, res(rowsG)));
        if isempty(rows), continue; end
        [R, tCommon] = stackCols(res(rows), 'rateNorm');
        col = spec.condSets(cs).colour;

        yyaxis(ax,'left')
        [m, ci] = meanCI(R);
        [wR, wT] = dispSmoothWidths(size(R,2), median(diff(tCommon)), cfg);
        m = movmean(m, wR);   % CI band stays unsmoothed
        gd = isfinite(m) & isfinite(ci);
        if any(gd)
            fill(ax, [tCommon(gd); flipud(tCommon(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                col, 'FaceAlpha',0.22, 'EdgeColor','none', 'HandleVisibility','off');
        end
        h(end+1) = plot(ax, tCommon, m, '-', 'Color', col, 'LineWidth', 3);   %#ok<SAGROW>
        lbl{end+1} = sprintf('%s rate (n=%d)', spec.condSets(cs).name, size(R,2));   %#ok<SAGROW>
        nUsed = max(nUsed, size(R,2));

        yyaxis(ax,'right')
        T = stackCols(res(rows), 'tongue');
        if cfg.traceNorm_tongue, T = normaliseCols(T, 'minmax'); end
% tongue traces are greyscale, never the condition colour, and
% solid rather than dashed: the first trial type is black and any further
% ones step towards grey, so two types read as black and grey.
        if nCS > 1
            tcol = repmat(0.55*(cs-1)/(nCS-1), 1, 3);
        else
            tcol = [0 0 0];
        end
        h(end+1) = plot(ax, tCommon, movmean(mean(T,2,'omitnan'), wT), '-', ...   % tongue smoothed
            'Color', tcol, 'LineWidth', 2);   %#ok<SAGROW>
        lbl{end+1} = sprintf('%s tongue', spec.condSets(cs).name);   %#ok<SAGROW>
    end
    if isempty(h), close(fh); continue; end

    yyaxis(ax,'left');   ylabel(ax, normLabel(cfg.traceNorm));
    ax.YAxis(1).Color = spec.condSets(1).colour;
    yyaxis(ax,'right');  ylabel(ax, 'Mean tongue length');
    ax.YAxis(2).Color = [0 0 0];

    xline(ax, 0, 'k--', 'LineWidth', 1);
    xlim(ax, cfg.traceXlim_s);
    xlabel(ax, 'Time from first lick contact (s)');
    aa = animalsInGroup(res, spec, g);
    title(ax, sprintf('%s | %d sessions, %d animals', ...
        spec.groupLabels{g}, nUsed, aa), 'FontSize', 10, 'FontWeight','normal');
    legend(ax, h, lbl, 'Location','northeast', 'Box','off', 'FontSize', 8);
    box(ax,'off');  set(ax, 'FontSize', 11);
    set(fh, 'Position', [60 60 480 420]);
end
end

%% FIGURE 2: THE 'PORT CONTACTS' FIGURE
% Drawn only where spec.showEpochFigure is true -- the Learning script. The
% other four tasks have a single time course and no epoch panel.
if spec.showEpochFigure
E = cell(1, nG);  anmOf = cell(1, nG);
for g = 1:nG
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe && r.condSet == 1, res));
    if isempty(rows), E{g} = [];  anmOf{g} = {};  continue; end
    E{g}     = cat(2, res(rows).epMean);   % nEp x nSess
    anmOf{g} = {res(rows).anm};
end

fh = figure('Name', sprintf('%s - mean spike rate by port contact', spec.name), 'Color','w');
ax = axes('Parent', fh);  hold(ax,'on')
x = 1:nEp;  h = gobjects(0);  lbl = {};  allY = [];

for g = 1:nG
    if isempty(E{g}), continue; end
    M = E{g}(:, any(isfinite(E{g}), 1));
    if isempty(M), continue; end
    m   = mean(M, 2, 'omitnan');
    sem = std(M, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(M), 2));
    col = spec.groupColours(g,:);

    gd = (isfinite(m) & isfinite(sem))';
    if any(gd)
        xg = x(gd);  hi = m(gd)' + sem(gd)';  lo = m(gd)' - sem(gd)';
        fill(ax, [xg, fliplr(xg)], [hi, fliplr(lo)], col, ...
            'FaceAlpha',0.20, 'EdgeColor','none', 'HandleVisibility','off');
        allY = [allY, hi, lo];   %#ok<AGROW>
    end
    h(end+1) = plot(ax, x, m, '-o', 'Color', col, 'MarkerFaceColor', col, ...
        'LineWidth', 2, 'MarkerSize', 6);   %#ok<SAGROW>
    lbl{end+1} = sprintf('%s (n=%d)', spec.groupLabels{g}, size(M,2));   %#ok<SAGROW>
end

set(ax, 'XTick', x, 'XTickLabel', epLabels, 'FontSize', 12);
xlim(ax, [0.7 nEp+0.3]);
xlabel(ax, 'Port contacts', 'FontSize', 13);
ylabel(ax, epochLabelY(cfg), 'FontSize', 13);
legend(ax, h, lbl, 'Location','northeast', 'Box','off', 'FontSize', 9);
box(ax,'off');
set(fh, 'Position', [560 60 360 430]);

%% STATISTICS: WILCOXON SIGNED RANK, PAIRED, THROUGHOUT
fprintf('\n%s\n  WILCOXON SIGNED RANK (paired)%s\n%s\n', repmat('=',1,72), ...
    ternary(cfg.useFDR, ', Benjamini-Hochberg across epochs', ', uncorrected'), repmat('-',1,72));
fprintf('  within group: epoch L vs epoch 1, paired across sessions\n');
for g = 1:nG
    if isempty(E{g}) || size(E{g},2) < 3, continue; end
    p = signrankVsFirst(E{g});
    if cfg.useFDR, p = bhFDR(p); end
    fprintf('    %-8s', spec.groupLabels{g});
    for e = 2:nEp
        fprintf('  %s vs %s: p=%s%s |', epLabels{e}, epLabels{1}, pstr(p(e)), ...
            ternary(isfinite(p(e)) && p(e) < cfg.alpha, ' *', ''));
    end
    fprintf('\n');
end

if ~isempty(spec.compareGroups)
    ga = spec.compareGroups(1);  gb = spec.compareGroups(2);
    fprintf('%s\n  between groups: %s vs %s at each epoch, paired BY ANIMAL\n', ...
        repmat('-',1,72), spec.groupLabels{ga}, spec.groupLabels{gb});
    [p, nPair] = signrankGroupsByAnimal(E{ga}, anmOf{ga}, E{gb}, anmOf{gb});
    if cfg.useFDR, p = bhFDR(p); end
    for e = 1:nEp
        fprintf('    %-8s n=%d pairs  p=%s%s\n', epLabels{e}, nPair(e), pstr(p(e)), ...
            ternary(isfinite(p(e)) && p(e) < cfg.alpha, '  *', ''));
    end
    drawEpochStars(ax, x, p, cfg.alpha, allY);
    drawGroupBracket(ax, nEp, E{ga}, E{gb});
end
fprintf('%s\n', repmat('=',1,72));
end
%% LOCAL FUNCTIONS (shared)

function L = epochLabels(cfg, rows)
% Axis labels: '1-3', '3-5', '5-8'.
    if ~isempty(cfg.epochLabels)
        L = cfg.epochLabels;
        assert(numel(L) == size(rows,1), 'cfg.epochLabels has %d entries for %d epochs.', ...
            numel(L), size(rows,1));
        return
    end
    L = cell(1, size(rows,1));
    for e = 1:size(rows,1)
        L{e} = sprintf('%d-%d', rows(e,1), rows(e,2));
    end
end

function s = normLabel(mode)
% The y label names THIS figure's normalisation, so a panel can never imply an
% amplitude claim the scaling already removed.
    switch lower(mode)
        case 'window', s = 'Mean spike rate (norm, GC-to-contact window)';
        case 'trace',  s = 'Mean spike rate (norm, whole trace)';
        otherwise,     s = 'Mean spike rate (Hz)';
    end
end

function s = epochLabelY(cfg)
    if cfg.epochDivideByFirst
        s = 'Mean spike rate (norm)';
    else
        s = normLabel(cfg.normMode);
    end
end

function c = medianContacts(contactCell, nWant)
% Median time of contact k across this condition's trials, k = 1..nWant. Trials
% that never reached contact k contribute nothing to it, so a late contact is a
% median over the trials that got there -- and NaN if none did.
    c = nan(nWant,1);
    n = numel(contactCell);
    M = nan(nWant, n);
    for i = 1:n
        v = contactCell{i};
        k = min(numel(v), nWant);
        if k > 0, M(1:k, i) = v(1:k); end
    end
    for k = 1:nWant
        v = M(k, isfinite(M(k,:)));
        if ~isempty(v), c(k) = median(v); end
    end
end

function y = scaleTrace(x, traceTime, winEdges, mode)
% Trim, then min-max, then hand the WHOLE trace back on the same scale.
% The min and max come from the trimmed window only -- [go cue, contact N] --
% but they are applied to every sample, so the trace outside the window is still
% drawn, just off the 0-1 range. That is what makes FIGURE 1 and the epoch
    switch lower(mode)
        case 'none'
            y = x;  return
        case 'window'
            m = traceTime >= winEdges(1) & traceTime <= winEdges(2);
        case 'trace'
            m = true(size(traceTime));
        otherwise
            error('normalisation mode must be ''window'', ''trace'' or ''none''.');
    end
    v = x(m & isfinite(x));
    if isempty(v), y = nan(size(x));  return; end
    lo = min(v);  hi = max(v);
    if hi <= lo, y = nan(size(x));  return; end
    y = (x - lo) ./ (hi - lo);
end

function [epMean, epN] = epochFromTrace(y, traceTime, dt, cMed, rows)
% Mean of the normalised trace inside each epoch. Epoch [a b] runs from the
% median time of contact a to the median time of contact b.
% Index arithmetic rather than a logical mask per epoch: the clock is uniform, so
% the sample range is a division. On a 1500-sample trace this is the difference
% between three comparisons over the whole vector and three subtractions.
    nEp    = size(rows,1);
    epMean = nan(nEp,1);
    epN    = zeros(nEp,1);
    t0 = traceTime(1);  nT = numel(traceTime);
    for e = 1:nEp
        a = cMed(rows(e,1));  b = cMed(rows(e,2));
        if ~isfinite(a) || ~isfinite(b) || b <= a, continue; end
        i0 = max(1,  floor((a - t0)/dt) + 1);
        i1 = min(nT, ceil( (b - t0)/dt) + 1);
        if i1 <= i0, continue; end
        v = y(i0:i1);
        v = v(isfinite(v));
        if isempty(v), continue; end
        epMean(e) = mean(v);
        epN(e)    = numel(v);
    end
end

function M = normaliseCols(M, mode)
% Column-wise rescale to [0,1]; a degenerate column becomes NaN rather than a
% divide by zero.
    if isempty(M) || strcmpi(mode,'none'), return; end
    mn = min(M, [], 1, 'omitnan');
    mx = max(M, [], 1, 'omitnan');
    rg = mx - mn;
    rg(rg == 0 | ~isfinite(rg)) = NaN;
    M  = (M - mn) ./ rg;
end

function [R, tCommon] = stackCols(rows, field)
% One column per (session, probe, condition), trimmed to the shortest common
% length. With one shared params block every session has the same clock, so the
% trim is a guard rather than a routine operation.
    R = [];  tCommon = [];
    if isempty(rows), return; end
    n = min(arrayfun(@(r) numel(r.(field)), rows));
    R = zeros(n, numel(rows));
    for i = 1:numel(rows)
        v = rows(i).(field);
        R(:,i) = v(1:n);
    end
    tCommon = rows(1).time(1:n);
end

function u = unitsInGroup(res, spec, g)
% Total neurons in one group, counting each (session, probe) ONCE however many
% condition sets it contributed.
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    if isempty(rows), u = 0;  return; end
    key = arrayfun(@(r) r.session*100 + r.probe, res(rows));
    [~, ia] = unique(key, 'stable');
    u = sum([res(rows(ia)).nUnits]);
end

function a = animalsInGroup(res, spec, g)
% Distinct animals contributing to one group, counted from the same rows
% unitsInGroup uses: a session's probe must be the one this group maps to.
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    if isempty(rows), a = 0;  return; end
    a = numel(unique({res(rows).anm}));
end


function [m, ci] = meanCI(M)
% Mean and 95% CI half-width across columns.
    m = mean(M, 2, 'omitnan');
    n = sum(isfinite(M), 2);
    s = std(M, 0, 2, 'omitnan');
    if exist('tinv','file'), tc = tinv(0.975, max(n-1,1)); else, tc = 1.96*ones(size(n)); end
    ci = tc .* s ./ sqrt(max(n,1));
    ci(n < 2) = NaN;
end

function p = signrankVsFirst(E)
% Epoch e against epoch 1, paired across sessions. Wilcoxon signed rank -- the
% test the figure and the console both name.
    nEp = size(E,1);
    p = nan(nEp,1);
    for e = 2:nEp
        a = E(e,:);  b = E(1,:);
        g = isfinite(a) & isfinite(b);
        if sum(g) < 3, continue; end
        if exist('signrank','file'), p(e) = signrank(a(g), b(g)); end
    end
end

function [p, nPair] = signrankGroupsByAnimal(Ea, anmA, Eb, anmB)
% Two groups at every epoch, PAIRED BY ANIMAL. Each animal contributes one
% session per learning day, so day 1 vs day 5 is within-animal and pairing is
% the correct test. An animal in only one group is dropped from the comparison
% rather than matched to someone else's session.
    nEp   = size(Ea,1);
    p     = nan(nEp,1);
    nPair = zeros(nEp,1);
    if isempty(Ea) || isempty(Eb), return; end
    shared = intersect(anmA, anmB);
    if isempty(shared), return; end
    A = nan(nEp, numel(shared));  B = nan(nEp, numel(shared));
    for i = 1:numel(shared)
        A(:,i) = Ea(:, find(strcmp(anmA, shared{i}), 1));
        B(:,i) = Eb(:, find(strcmp(anmB, shared{i}), 1));
    end
    for e = 1:nEp
        g = isfinite(A(e,:)) & isfinite(B(e,:));
        nPair(e) = sum(g);
        if nPair(e) < 3, continue; end
        if exist('signrank','file'), p(e) = signrank(A(e,g), B(e,g)); end
    end
end

function q = bhFDR(p)
% Benjamini-Hochberg over the finite entries only.
    q = p;
    ok = find(isfinite(p));
    if isempty(ok), return; end
    [ps, ord] = sort(p(ok));
    m = numel(ps);
    adj = ps .* m ./ (1:m)';
    for i = m-1:-1:1, adj(i) = min(adj(i), adj(i+1)); end
    out = nan(m,1);  out(ord) = min(adj, 1);
    q(ok) = out;
end

function drawEpochStars(ax, x, p, alpha, allY)
    allY = allY(isfinite(allY));
    if isempty(allY), return; end
    yTop = max(allY);  yBot = min(allY);
    pad  = 0.06 * max(yTop - yBot, eps);
    for e = 1:numel(x)
        if isfinite(p(e)) && p(e) < alpha
            text(ax, x(e), yTop + pad, '*', 'HorizontalAlignment','center', ...
                'FontSize', 16, 'FontWeight','bold', 'Color', [0 0 0]);
        end
    end
    ylim(ax, [yBot - 2*pad, yTop + 4*pad]);
end

function drawGroupBracket(ax, nEp, Ea, Eb)
% The vertical bracket at the last epoch, spanning the two compared groups.
    if isempty(Ea) || isempty(Eb), return; end
    ma = mean(Ea(nEp,:), 'omitnan');
    mb = mean(Eb(nEp,:), 'omitnan');
    if ~isfinite(ma) || ~isfinite(mb), return; end
    xb = nEp + 0.14;  tick = 0.05;
    plot(ax, [xb xb], [min(ma,mb) max(ma,mb)], 'k-', 'LineWidth',1.2, 'HandleVisibility','off');
    plot(ax, [xb-tick xb], [ma ma], 'k-', 'LineWidth',1.2, 'HandleVisibility','off');
    plot(ax, [xb-tick xb], [mb mb], 'k-', 'LineWidth',1.2, 'HandleVisibility','off');
    xlim(ax, [0.7 nEp + 0.35]);
end

function s = pstr(p)
    if ~isfinite(p), s = '  n/a ';
    elseif p < 1e-4, s = '<1e-4';
    else, s = sprintf('%.4f', p);
    end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function s = groupOfProbe(spec, sessionIdx, probeIdx)
% Which group this probe carries in THIS session -- the region for the region
% scripts, the learning day for the learning script. Read off the maps rather
% than assumed, because probe 2 is ALM in one session and M1 in the next.
    hit = cellfun(@(m) m(sessionIdx) == probeIdx, spec.groupMaps);
    if ~any(hit)
        s = 'no group';
    else
        s = strjoin(spec.groupLabels(hit), '+');
    end
end

function S = srLoadSession(loader, dateStr, spec, params, cfg)
% Load one session, with an on-disk cache keyed on every params field that can
% epoch definition or a normalisation re-runs in seconds off the warm cache.
% The cache is SEPARATE from the decoding scripts' cache (different directory)
% because params.lowFR differs -- 0.001 here against 0.01 there -- and mixing
% them would silently serve the wrong unit set.
    key = cacheKey(loader, dateStr, params);

    useCache = cfg.useCache && ~isempty(cfg.cacheDir);
    if useCache
        if ~exist(cfg.cacheDir, 'dir'), mkdir(cfg.cacheDir); end
        fname = fullfile(cfg.cacheDir, sprintf('%s_%s.mat', spec.name, matlab.lang.makeValidName(dateStr)));
        if exist(fname, 'file')
            C = load(fname, 'S', 'key');
            if isfield(C,'key') && strcmp(C.key, key)
                S = C.S;
                logf('  [cache] hit  %s\n', dateStr);
                return
            end
            logf('  [cache] stale %s (params changed) -- reloading\n', dateStr);
        end
    end

    tLoad = tic;
    meta  = slimMeta(loader, dateStr);
    params.probe = {meta.probe};
    params.cluid = {};
    [obj, params, kin] = slimToLegacy(meta, params);
    obj = obj(1);  params = params(1);

    featCol = find(strcmp(kin.featLeg, 'tongue_length'), 1);
    assert(~isempty(featCol), 'tongue_length not found in kin.featLeg for %s', dateStr);

    S.time        = obj.time(:);
    S.trialdat    = single(obj.trialdat);
    S.nUnitsTotal = size(obj.psth, 2);
    S.tongueRaw   = kin.dat(:,:,featCol);
    S.cluid       = params.cluid;
    S.trialid     = params.trialid;
    S.anm         = obj.pth.anm;
    S.date        = obj.pth.dt;
    S.Ntrials     = obj.bp.Ntrials;
    S.lickL       = obj.bp.ev.lickL;
    S.goCue       = obj.bp.ev.goCue;
    S.trialTypes  = obj.bp.trialTypes;
    S.hit         = obj.bp.hit;
    logf('  [load]  %s in %.1f s\n', dateStr, toc(tLoad));

    if useCache
        save(fname, 'S', 'key', '-v7.3');
    end
end

function key = cacheKey(loader, dateStr, params)
% Everything that can change what gets loaded, and nothing that cannot.
    k.loader  = func2str(loader);
    k.date    = dateStr;
    k.align   = params.alignEvent;
    k.behav   = params.behav_only;
    k.warp    = params.timeWarp;
    k.nLicks  = params.nLicks;
    k.lowFR   = params.lowFR;
    k.quality = params.quality;
    k.cond    = params.condition;
    k.tmin    = params.tmin;
    k.tmax    = params.tmax;
    k.dt      = params.dt;
    k.smooth  = params.smooth;
    k.traj    = params.traj_features;
    k.featVar = params.feat_varToExplain;
    k.NVar    = params.N_varToExplain;
    k.adv     = params.advance_movement;
    k.fcut    = params.fcut;
    k.condN   = params.cond;
    k.method  = params.method;
    k.fa      = params.fa;
    k.bctype  = params.bctype;
    k.source  = 'slimExport';   % never reuse a cache entry built by the pipeline
    key = jsonencode(k);
end

function [wRate, wTongue] = dispSmoothWidths(nSess, dt, cfg)
% Moving-average widths in samples, from
% main_spikeRate_manySessions_R1_forScience.m:
% rescaled from 5 ms samples to params.dt so the width in SECONDS is kept.
    if ~cfg.dispSmooth || nSess < 1 || ~isfinite(dt) || dt <= 0
        wRate = 1;  wTongue = 1;  return
    end
    s = cfg.dispSmoothRefDt / dt;
    wRate   = max(1, round(round(cfg.dispSmoothRateBase   * cfg.dispSmoothBaseN / nSess) * s));
    wTongue = max(1, round(round(cfg.dispSmoothTongueBase * cfg.dispSmoothBaseN / nSess) * s));
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

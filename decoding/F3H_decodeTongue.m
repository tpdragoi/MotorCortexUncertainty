%% F3H_decodeTongue.m
%  Ridge decoder of tongue length from population spike rates, Double Reward Task.
%  Produces Fig. 3H (test sets R1 and R6).
%  Spike rates are z-scored and lagged and fit per session by ridge regression,
%  with lambda chosen by k-fold cross-validation over training trials. The
%  decoding index per lick contact is 1 - RMSE(decoded) / RMSE(zero baseline).
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    via shared\slimMeta and shared\slimToLegacy; Data = dataRoot in setPaths.m
%  ANALYSIS SETTINGS
%    params.alignEvent 'firstLick' | dt 1/300 s | smooth 10 bins | window -2.5 to 5 s
%    params.quality {'good','excellent'} | lowFR 0.01 Hz
%  Run the whole file.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.

%% RUN SETTINGS (edit these two)
% Number of full passes through the analysis (1 = a single run). Each pass
% re-draws the train/test split and the cross-validation folds and gets its own
% figure windows; passes are not averaged or pooled.
RUN.nRepeats = 2;

% Cross-validation folds for choosing lambda, within the training set.
RUN.cvFolds  = 5;

%% PATHS
% Sessions are read from the exported obj/kin files in the data folder set in
% setPaths.m, through slimMeta / slimToLegacy (in shared/). slimToLegacy rebuilds
% obj/params/kin for this script's params from the exported spike times; the
% pipeline functions it needs are in shared\pipelineCopies.
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
assert(exist(fullfile(repoRoot, 'setPaths.m'), 'file') == 2, ...
    'Cannot find setPaths.m. Run this script from its file, or cd to the repository root first.');
addpath(repoRoot);
cfgPaths = setPaths();   % data locations are set once, in setPaths.m
spec.dataDir = '';   % raw data folder not used

%% SETTINGS
rewardLickB = 6;   % rewardedLick value of the second reward-lick condition

%% TIMING: the three windows, all in seconds
% All times are relative to t = 0, the first lickport contact after the go cue
% (params.alignEvent = 'firstLick').

% FITTING WINDOW (per trial, variable length): from preGoCue_s before the go cue
% to postTongueRetract_s after the tongue retracts from contact 1. Contact 1
% occurs part way through the protrusion, so the window closes at the end of the
% run of visible (non-zero) tongue length that contains t = 0, plus the pad.
cfg.preGoCue_s          = 0.20;   % window opens this long BEFORE the go cue
cfg.postTongueRetract_s = 0.005;   % window closes this long after the retraction
cfg.minWindowSamples    = 10;   % trials with a shorter window are dropped

%% FITTING WINDOW MODE
% 'contact1'        go cue - preGoCue_s  ->  retraction after contact 1
%                   + postTongueRetract_s
% 'bout'            contact 2 + winBoutStartPad_s  ->  bout-ending contact
%                   + winBoutEndPad_s. The bout ends at the first contact not
%                   followed by another within winBoutGap_s; contact 1 is
%                   outside the window.
% 'lick2ToBoutEnd'  retraction of lick 2 + winBoutStartPad_s  ->  same end as
%                   'bout'
% Under 'contact1' the model is trained on the first protrusion, so contacts 2-8
% are extrapolation; under 'bout' contact 1 is the extrapolation. Everything
% else (z-scoring, lags, lambda rule, decoding index) is the same in all modes.
cfg.fitWindowMode     = 'contact1';   % 'contact1' | 'bout' | 'lick2ToBoutEnd'
cfg.winBoutStartPad_s = 0.02;   % window opens this long AFTER contact 2
cfg.winBoutEndPad_s   = 0.20;   % window closes this long AFTER the bout's last contact
cfg.winBoutGap_s      = 0.75;   % a gap longer than this is what ENDS the bout

% Runs of visible tongue separated by this many zero samples or fewer are merged
% into one protrusion, so a brief tracking dropout does not split one lick into
% two. 2 samples = 6.7 ms at dt = 1/300 s; 0 disables the merge.
cfg.mergeGapSamples = 3;

% NEURAL Z-SCORING WINDOW
cfg.zscoreWin_s = [-2.0 2.0];   % s, around the first port contact
cfg.minBaselineStd = 1e-6;   % units with a smaller z-score SD are dropped

% NEURAL LAG CONTEXT (lag set is -preBins : +postBins, inclusive)
cfg.lagPre_s  = 0.06;
cfg.lagPost_s = 0.04;

% ANALYSIS WINDOW for the per-contact decoding index
cfg.analysisWin_s = [-0.09 3.0];

%% TARGET: tongue length
% Tongue length is normalized between two percentiles of the visible samples.
% Frames where the tongue is not visible are set to length zero (an observation,
% not missing data).
cfg.tongueZeroMode = 'p_low';   % 'p_low' | 'anatomical'
cfg.normPctLow  = 2;   % lower normalization percentile
cfg.normPctHigh = 99;   % upper normalization percentile

%% RIDGE
cfg.ridgeGrid       = logspace(-2, 6, 80);
cfg.standardizeCols = true;   % penalize standardized coefficients (required by ridgeImpl = 'matlab')

cfg.ridgeImpl = 'eig';   % 'eig' | 'matlab'
% 'eig' is the closed-form solution from one eigendecomposition of the Gram
% matrix; 'matlab' calls ridge(), which solves once per lambda. Same estimator
% (see ridgePath); 'eig' is much faster when there are many predictors.
cfg.ridgeCheckEquivalence = false;   % compare 'eig' and ridge() on the first probe fit

% LAMBDA: 'cv5' = k-fold cross-validation inside the training set
% (k = cfg.cvFolds), with folds over whole trials (see selectLambdaCV5).
% Test trials are never used for fitting or for choosing lambda.
cfg.lambdaRule  = 'cv5';   % 'cv5' (k-fold CV, k = cfg.cvFolds)
cfg.cvFolds     = RUN.cvFolds;   % set in the RUN block at the top
cfg.cvSeed      = [];   % [] = folds in trial order; a number seeds the fold shuffle
RUN.cvSeed0     = cfg.cvSeed;   % base seed; pass k adds k - 1
cfg.lambdaFixed = [];   % [] = choose lambda by cfg.lambdaRule; a number pins it (nearest grid point)

%% CONTACT DETECTION + DECODING INDEX
% decodingIndex(L) = 1 - RMSE_decoded(L) / RMSE_zero(L), where RMSE_zero is the
% error of a constant-zero prediction over contact L's samples. RMSE_zero scales
% with the excursion, so an amplitude decline alone lowers the index.
cfg.nLicksAnalyze     = 8;
cfg.minContactSamples = 6;   % samples; 20 ms at dt = 1/300 s

%% BOUT REQUIREMENT (one toggle, on or off)
% Keep a trial only if it contains a bout: at least boutMinLicks consecutive port
% contacts, no more than boutMaxILI_s apart, all within boutWin_s of the go cue.
% Applied to the whole trial pool before the train/test split, so the model is
% trained on the same kind of trial it is scored on.
cfg.requireBout   = true;   % true = apply the bout requirement
cfg.boutMinLicks  = 6;
cfg.boutWin_s     = 1.50;
cfg.boutMaxILI_s  = 0.25;

%% STATISTICS
cfg.alpha      = 0.05;

%% SESSION CACHE
% Loaded sessions are cached on disk. The cache key holds every params field that
% affects loading, so changing one forces a reload. Budget roughly
% nTime x nUnits x nTrials x 4 bytes per session (trialdat is stored as single).
% Delete the folder to force a reload.
cfg.useCache = true;
cfg.cacheDir = fullfile(cfgPaths.cacheRoot, 'tlDecodeCache');

%% HOUSEKEEPING
cfg.rngSeed          = [];
cfg.sessionsToRun    = [];   % [] = all

% DIAGNOSTIC FIGURES. cfg.diagSessions lists session numbers explicitly; [] picks
% them by spec.diagRule:
%   'random'  cfg.nDiagSessions sessions drawn at random (seeded by cfg.rngSeed)
%   'groups'  every session in spec.diagGroups
% Each diagnostic session is a tab in two figures: the six-panel fit diagnostics
% and the warped lick train.
cfg.diagSessions     = [];
cfg.nDiagSessions    = 0;

% Sessions on the warped lick-train figure:
%   true   every fitted session
%   false  the same sessions as the six-panel diagnostics
cfg.warpAllSessions  = false;
cfg.warnPredictorCount = 8000;

% Warped lick-train layout: first slot at sample 10, slots every 50 samples (a
% 20-sample NaN gap between licks); up to 10 licks extracted, the first 8 drawn.
cfg.warpedSegmentLength = 30;   % samples each excursion is resampled onto
cfg.warpedSlotFirst     = 10;   % where lick 1's slot starts
cfg.warpedSlotSpacing   = 50;   % slot pitch; spacing - segment = the NaN gap
cfg.maxLicksToExtract   = 10;
cfg.maxLickShow         = 8;

%   warpedInterp      'pchip' or 'linear' ('linear' cannot overshoot at segment ends)
%   warpedTickCenter  false = contact label at the left edge of each slot,
%                     true  = centered under the excursion
cfg.warpedInterp     = 'pchip';   % 'pchip' | 'linear'
cfg.warpedTickCenter = false;

% With two test sets, draw each on its own stacked axes (actual black, predicted
% in that test set's color, labeled in-plot). Ignored with one test set.
cfg.warpSplitTestSets = true;

%% PARAMS (field names read by the loaders)
params.alignEvent = 'firstLick';
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 8;
params.lowFR  = 0.01;   % minimum mean firing rate, Hz
params.quality    = {'good','excellent',' good','good '};   % good + excellent units (the exported unit list)
params.tmin   = -2.5;
params.tmax   = 5;
params.dt     = 1/300;   % 300 Hz analysis grid (video acquired at 400 Hz)
params.smooth = 10;   % bins; 10 x (1/300 s) = 33 ms
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

% Condition 8 is the R1 pool; condition 9 is the second reward-lick pool
% (rewardedLick == rewardLickB).
b = rewardLickB;
params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 1'};
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 1 & rewardedLick == %d', b)};
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 2 & rewardedLick == %d', b)};
params.condition(end+1) = {sprintf('hit==1 & trialTypes == 3 & rewardedLick == %d', b)};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};   % 8
params.condition(end+1) = {sprintf('hit==1 & rewardedLick == %d', b)};   % 9
params.condition(end+1) = {'hit==1'};   % 10

%% STUDY: R1 vs R6
spec.name = 'R1R6';

spec.sessionDates = { ...
    '2023-06-03','2023-06-04','2023-06-05','2023-06-06', ...   % YH2
    '2023-05-08','2023-05-10','2023-05-11','2023-05-12','2023-06-03','2023-06-04', ...   % YH1
    '2023-12-06','2023-12-07','2023-12-09','2023-12-13','2023-12-14','2023-12-16', ...   % TD4f
    '2024-01-15','2024-01-16','2024-01-17','2024-01-18','2024-01-22','2024-01-23','2024-01-24', ...   % TD7f
    '2025-06-25','2025-06-26','2025-06-27','2025-06-28','2025-06-29', ...   % TD24
    '2025-07-14','2025-07-15','2025-07-16','2025-07-17','2025-07-18' };   % TD25

spec.sessionLoaders = [ ...
    repmat({@loadYH2_neur},   1,4), repmat({@loadYH1_neur},   1,6), ...
    repmat({@loadTD4f_neur},  1,6), repmat({@loadTD7f_neur},  1,7), ...
    repmat({@loadTD24_neural},1,5), repmat({@loadTD25_neural},1,5) ];

spec.groupMaps = { ...
    [1 1 1 1  0 1 1 1 1 1  1 1 1 0 0 1  1 0 1 1 2 0 2  0 0 0 0 0  0 0 0 0 0] , ...   % M1
    [2 0 2 2  0 0 0 0 2 2  0 2 2 1 0 0  0 0 0 0 0 0 0  2 0 2 2 0  1 1 1 1 0] };   % ALM
spec.groupLabels = {'M1','ALM'};

% TRAIN pool = conditions 8 and 9 together (R1 + R6), type-stratified quota.
spec.poolCondIdx = [8 9];
spec.splitRule   = 'typeStratified';
spec.trainFrac   = 0.60;   % fraction of the pool used for training
spec.testSets    = struct('name', {'R1','R6'}, 'condIdx', {8, 9});

spec.trialCaps = { 'TD4f', '2023-12-16', 203 ; ...
                   'TD24d','2025-06-30', 200 ; ...
                   'TD25d','2025-07-17', 199 };
spec.singleProbeSessions = { 'TD24d','2025-06-30' ; 'TD24d','2025-07-01' ; 'TD25d','2025-07-14' };
spec.lickCountWin_s      = 1.25;
spec.warpShowLicks    = 8;   % contacts LABELED on the warped figure

%% bout filter + diagnostics
spec.applyBoutFilter = true;   % cfg.requireBout is the on/off switch
spec.diagRule        = 'random';   % 'random' = cfg.nDiagSessions sessions at random
spec.diagGroups      = [];   % used only by diagRule 'groups'

%% colors
% testSets are {'R1','R6'} in that order: R1 PINK, R6 PURPLE.
spec.tsColors  = [0.95 0.45 0.70 ;   % R1  pink
                   0.50 0.15 0.65];   % R6  purple
spec.grpColors = [];
spec.predColor = [0.50 0.15 0.65];   % prediction on the diagnostics

%% figure behavior
spec.overlayGroups = false;   % false = never draw two groups on one axes
spec.showFigCI     = true ;   % Figure 1: decoding index with error bars
spec.showFigDelta  = false;   % Figure 3: delta from contact 1
spec.deltaGroups   = [];   % groups on the delta figures ([] = all)
spec.figRef       = 'Fig. 3H';   % figure panel this file produces
spec.vsC1Tail     = 'both';   % 'both' | 'left' | 'right'
spec.vsC1Test     = false;   % per-contact test vs contact 1
spec.betweenTest  = false;   % per-contact test BETWEEN the two trial types
spec.groupTest    = 'crossTask';   % 'none' | 'ranksum' | 'crossTask'
spec.groupSummaryLicks = 3:8;   % contacts averaged into the single summary value
spec.summaryFile  = fullfile(cfgPaths.cacheRoot, 'decodeSummary_tongue_r16.mat');
spec.crossTaskFile = fullfile(cfgPaths.cacheRoot, 'decodeSummary_tongue_r1.mat');
spec.compareGroups = [];   % [a b] = per-contact paired test between two groups

%% MAIN LOOP
nSess = numel(spec.sessionDates);
if ~isempty(cfg.rngSeed), rng(cfg.rngSeed); end

if isempty(cfg.sessionsToRun)
    runList = 1:nSess;
else
    runList = unique(cfg.sessionsToRun(cfg.sessionsToRun >= 1 & cfg.sessionsToRun <= nSess));
end
assert(~isempty(runList), 'cfg.sessionsToRun left no valid sessions.');

% Restrict the group maps to the run set, then list the probes each session needs.
runMask = false(1, nSess);  runMask(runList) = true;
for d = 1:numel(spec.groupMaps)
    assert(numel(spec.groupMaps{d}) == nSess, '%s map has %d entries, %d sessions.', ...
        spec.groupLabels{d}, numel(spec.groupMaps{d}), nSess);
    spec.groupMaps{d}(~runMask) = 0;
end

probesOf = cell(1, nSess);
for s = runList
    pr = cellfun(@(m) m(s), spec.groupMaps);
    probesOf{s} = unique(pr(pr > 0));
end

nTS  = numel(spec.testSets);
nLA  = cfg.nLicksAnalyze;
nK   = numel(cfg.ridgeGrid);
logK = log10(cfg.ridgeGrid(:));

% ---------------- WHICH SESSIONS GET A DIAGNOSTIC FIGURE ----------------
% Chosen up front, not as the loop goes, so the random draw does not depend on
% which sessions happen to survive the trial-count check.
if ~isempty(cfg.diagSessions)
    diagList = cfg.diagSessions(:)';
elseif strcmpi(spec.diagRule, 'groups')
    diagList = [];
    for g = spec.diagGroups(:)'
        diagList = [diagList, find(spec.groupMaps{g} > 0)];   %#ok<AGROW>
    end
    diagList = unique(diagList);
else
    cand  = runList(~cellfun(@isempty, probesOf(runList)));
    nPick = min(cfg.nDiagSessions, numel(cand));
    diagList = sort(cand(randperm(numel(cand), nPick)));
end

applyBout = cfg.requireBout && spec.applyBoutFilter;

fprintf('\n=== %s | %d of %d sessions | %d (session,probe) pairs ===\n', ...
    spec.name, numel(runList), nSess, sum(cellfun(@numel, probesOf)));
switch lower(cfg.fitWindowMode)
    case 'contact1'
        logf('  fitting window [contact1]: go cue -%.0f ms  ->  tongue retracts after contact 1, +%.0f ms (per trial)\n', ...
            1000*cfg.preGoCue_s, 1000*cfg.postTongueRetract_s);
    case 'bout'
        logf('  fitting window [bout]: contact 2 +%.0f ms  ->  bout-end contact +%.0f ms (bout ends on a gap > %.2f s)\n', ...
            1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s, cfg.winBoutGap_s);
        fprintf('    CONTACT 1 IS OUTSIDE THE TRAINING WINDOW in this mode -- it is the extrapolation, not contacts 2-8.\n');
    case 'lick2toboutend'
        logf('  fitting window [lick2ToBoutEnd]: lick 2 RETRACTS +%.0f ms  ->  bout-end contact +%.0f ms (bout ends on a gap > %.2f s)\n', ...
            1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s, cfg.winBoutGap_s);
        fprintf('    CONTACT 1 AND LICK 2 ARE OUTSIDE THE TRAINING WINDOW -- contact 1 is the extrapolation, not contacts 3-8.\n');
    otherwise
        error('cfg.fitWindowMode must be ''contact1'', ''bout'' or ''lick2ToBoutEnd''.');
end
logf('  z-score window: %.1f to %.1f s | lag context %.0f to +%.0f ms | baseline: ZERO only\n', ...
    cfg.zscoreWin_s(1), cfg.zscoreWin_s(2), -1000*cfg.lagPre_s, 1000*cfg.lagPost_s);
if applyBout
    logf('  BOUT FILTER ON: >= %d contacts in a row, ILI <= %.0f ms, all within %.2f s of the go cue\n', ...
        cfg.boutMinLicks, 1000*cfg.boutMaxILI_s, cfg.boutWin_s);
elseif cfg.requireBout
    logf('  bout filter: not applied in this script (spec.applyBoutFilter = false)\n');
else
    logf('  bout filter: OFF (cfg.requireBout = false)\n');
end
logf('  lambda: chosen on the training trials (cfg.lambdaRule)\n');
logf('  diagnostics: sessions %s\n', mat2str(diagList));

% Sessions on the warped lick-train figure (see cfg.warpAllSessions).
if cfg.warpAllSessions
    warpList = runList(~cellfun(@isempty, probesOf(runList)));
else
    warpList = diagList;
end
logf('  warped figure: sessions %s\n', mat2str(warpList));

% Colors come from the study block.
tsCols  = spec.tsColors;
predCol = spec.predColor;

%% REPEAT PASSES
for repIdx = 1:RUN.nRepeats
logf('\n\n%s\n', repmat('#',1,78));
logf('#  PASS %d of %d   (%s)\n', repIdx, RUN.nRepeats, spec.name);
logf('%s\n\n', repmat('#',1,78));

% Re-draw BOTH sources of randomness for this pass. The trial split uses the
% global stream; the cross-validation folds use cfg.cvSeed (selectLambdaCV5
% builds its own RandStream from it). Offsetting both by repIdx makes each pass
% an independent draw while keeping the whole run reproducible: pass k is
% always the same pass k.
if isempty(cfg.rngSeed)
    rng(90210 + 1000*repIdx);
else
    rng(cfg.rngSeed + repIdx - 1);
end
cfg.cvSeed = RUN.cvSeed0 + repIdx - 1;
logf('[pass %d] split seed %d | cv fold seed %s | %d-fold\n\n', ...
    repIdx, 90210 + 1000*repIdx, mat2str(cfg.cvSeed), cfg.cvFolds);
if RUN.nRepeats > 1
    passTag = sprintf('pass %d/%d - ', repIdx, RUN.nRepeats);
else
    passTag = '';
end

% The two diagnostic figures are created ONCE and get a tab per session.
diagFig = [];  diagTG = [];
warpFig = [];  warpTG = [];

res = struct([]);  nFitted = 0;

sessCount = 0;
for sessionIdx = runList
    sessCount = sessCount + 1;
    fprintf('Session %d\n', sessCount);
    if isempty(probesOf{sessionIdx}), continue; end

    % ---------------- LOAD (cached; see tlLoadSession) ----------------
    S = tlLoadSession(spec.sessionLoaders{sessionIdx}, spec.sessionDates{sessionIdx}, spec, params, cfg);
    anm = S.anm;  dte = S.date;
    traceTime = S.time;
    nT = numel(traceTime);

    % X and y are paired by LINEAR INDEX into (nT x nTrials) layouts. If the two
    % arrays disagree about nT that pairing silently offsets the target against the
    % predictors. Cheap to check, fatal if wrong.
    assert(size(S.trialdat,1) == nT, 'sess %d: trialdat %d time samples, time %d.', ...
        sessionIdx, size(S.trialdat,1), nT);
    assert(size(S.tongueRaw,1) == nT, 'sess %d: kinematics %d time samples, time %d.', ...
        sessionIdx, size(S.tongueRaw,1), nT);

    % ---------------- TARGET: TONGUE LENGTH ----------------
    tlRaw = S.tongueRaw;
    visMask = ~isnan(tlRaw);
    visVals = tlRaw(visMask);
    pLo = prctile(visVals, cfg.normPctLow);
    pHi = prctile(visVals, cfg.normPctHigh);
    switch cfg.tongueZeroMode
        case 'p_low',      zeroRef = pLo;
        case 'anatomical', zeroRef = 0;
        otherwise, error('cfg.tongueZeroMode must be ''p_low'' or ''anatomical''.');
    end
    tongueLen = min(max((tlRaw - zeroRef) ./ (pHi - zeroRef), 0), 1);
    tongueLen(~visMask) = 0;   % tongue not visible -> length 0

    % ---------------- PER-TRIAL FITTING WINDOWS ----------------
    % Computed once for every trial in the session, then indexed by whichever set
    % needs it. goCueRel is negative: the cue precedes contact 1.
    maxUsable = min([size(S.tongueRaw,2), size(S.trialdat,3), S.Ntrials]);
    capRow = find(strcmp(anm, spec.trialCaps(:,1)) & strcmp(dte, spec.trialCaps(:,2)), 1);
    if ~isempty(capRow), maxUsable = min(maxUsable, spec.trialCaps{capRow,3}); end

    goCueRel  = nan(maxUsable,1);
    winCell   = cell(maxUsable,1);
    lickCount = zeros(maxUsable,1);
    boutRun   = zeros(maxUsable,1);   % longest consecutive-contact run in the window

    for tr = 1:maxUsable
        lk = S.lickL{tr};
        if isempty(lk), continue; end
        post = sort(lk(lk > S.goCue(tr)));
        if isempty(post), continue; end   % no contact -> cannot align
        goCueRel(tr) = S.goCue(tr) - post(1);   % t = 0 is post(1)
        lickCount(tr) = sum(post - S.goCue(tr) <= spec.lickCountWin_s);
        boutRun(tr)   = longestLickRun(lk, S.goCue(tr), cfg);
        winCell{tr} = fitWindowForTrial(tongueLen(:,tr), traceTime, goCueRel(tr), post - post(1), cfg);
    end
    hasWin = ~cellfun(@isempty, winCell);

    % ---------------- TRIAL POOL ----------------
    poolTrials = unique(cell2mat(S.trialid(spec.poolCondIdx)'));
    poolTrials = poolTrials(poolTrials <= maxUsable);
    poolTrials = poolTrials(hasWin(poolTrials));

    % ---------------- BOUT FILTER (cfg.requireBout) ----------------
    % Applied to the WHOLE pool, before the train/test split, so it removes training
    % trials too and the model is trained on the same kind of trial it is scored on.
    if applyBout
        nBefore = numel(poolTrials);
        poolTrials = poolTrials(boutRun(poolTrials) >= cfg.boutMinLicks);
        logf('  [bout]  %d of %d pool trials have >= %d contacts in a row within %.2f s (ILI <= %.0f ms); %d dropped\n', ...
            numel(poolTrials), nBefore, cfg.boutMinLicks, cfg.boutWin_s, ...
            1000*cfg.boutMaxILI_s, nBefore - numel(poolTrials));
    end

    % ---------------- MEAN TONGUE LENGTH INSIDE THE FITTING WINDOW ----------------
    % Every trial's own window, laid back on the real time axis and averaged across
    % trials at each sample. Windows are variable length (the go-cue latency and the
    % contact-1 duration both vary), so a sample is averaged over however many trials
    % actually cover it -- tlWinN records that count.
    Mtl = nan(nT, numel(poolTrials));
    for tt = 1:numel(poolTrials)
        ki = winCell{poolTrials(tt)};
        Mtl(ki,tt) = tongueLen(ki, poolTrials(tt));
    end
    tlWinMean = mean(Mtl, 2, 'omitnan');
    tlWinN    = sum(isfinite(Mtl), 2);
    clear Mtl

    % ---------------- TRIAL SPLIT ----------------
    switch spec.splitRule
        case 'typeStratified'
            [trainTrials, testTrials] = splitTypeStratified(poolTrials, S.trialTypes, spec.trainFrac);
        case 'trialType'
            typePool = poolTrials(S.trialTypes(poolTrials) == spec.trainTrialType);
            nTake    = min(numel(typePool), floor(numel(poolTrials)*spec.trainFrac));
            if nTake < 8, trainTrials = []; testTrials = []; else
                trainTrials = randsample(typePool, nTake, false);
                testTrials  = poolTrials(~ismember(poolTrials, trainTrials));
            end
        case 'lickCount'
            testTrials  = poolTrials(lickCount(poolTrials) >= spec.lickCountMin);
            trainTrials = poolTrials(lickCount(poolTrials) <  spec.lickCountMin);
        otherwise
            error('Unknown spec.splitRule ''%s''.', spec.splitRule);
    end

    nTrain = numel(trainTrials);  nTest = numel(testTrials);
    if nTrain < 10 || nTest < 5
        logf('\n[skip] sess %2d %s %s: %d train / %d test trials\n', sessionIdx, anm, dte, nTrain, nTest);
        continue
    end

    for probeIdx = probesOf{sessionIdx}(:)'

        regionStr = groupOfProbe(spec, sessionIdx, probeIdx);
        logf('\n===== sess %2d | %s %s | probe %d | %s =====\n', sessionIdx, anm, dte, probeIdx, regionStr);
        logf('  [split] pool %d (bound %d) | train %d | test %d\n', numel(poolTrials), maxUsable, nTrain, nTest);

        % ---------------- UNITS ON THIS PROBE ----------------
        isSingle = any(strcmp(anm, spec.singleProbeSessions(:,1)) & strcmp(dte, spec.singleProbeSessions(:,2)));
        if isSingle
            unitIds = 1:size(S.cluid,1);
        elseif probeIdx == 1
            unitIds = 1:size(S.cluid{1,1},1);
        else
            unitIds = size(S.cluid{1,1},1)+1 : S.nUnitsTotal;
        end

        % ---------------- PREDICTORS: z-score, then lag ----------------
        zIdx = find(traceTime >= cfg.zscoreWin_s(1) & traceTime <= cfg.zscoreWin_s(2));
        assert(~isempty(zIdx), 'cfg.zscoreWin_s selects no samples.');

        spikes = S.trialdat(:, unitIds, :);
        B  = reshape(spikes(zIdx,:,:), [], numel(unitIds));
        mu = mean(B, 1, 'omitnan');
        sd = std(B, 0, 1, 'omitnan');
        clear B

        dead = ~isfinite(sd) | sd < cfg.minBaselineStd;
        if any(dead)
            logf('  [units] dropping %d/%d unit(s) with z-score sd < %g\n', sum(dead), numel(unitIds), cfg.minBaselineStd);
            unitIds(dead) = [];  spikes = spikes(:,~dead,:);  mu = mu(~dead);  sd = sd(~dead);
        end
        nUnits = numel(unitIds);
        spikes = (spikes - reshape(mu,1,[],1)) ./ reshape(sd,1,[],1);

        lags = -round(cfg.lagPre_s/params.dt) : round(cfg.lagPost_s/params.dt);
        nPred = nUnits * numel(lags);
        if nPred > cfg.warnPredictorCount
            warning('p = %d predictors: the Gram matrix is %.1f GB.', nPred, 8*nPred^2/1e9);
        end

        winTrain = winCell(trainTrials);
        winTest  = winCell(testTrials);

        [Xtrain, rowTrain] = laggedDesign(spikes(:,:,trainTrials), lags, winTrain);
        Xtest              = laggedDesign(spikes(:,:,testTrials),  lags, winTest);
        yTrain = stackWindows(tongueLen, winTrain, trainTrials);
        yTest  = stackWindows(tongueLen, winTest,  testTrials);

        okTr = all(isfinite(Xtrain),2) & isfinite(yTrain);
        okTe = all(isfinite(Xtest), 2) & isfinite(yTest);
        winLen = cellfun(@numel, winTrain);

        logf('  [design] %d units x %d taps = %d predictors | window %d-%d samples (%.0f-%.0f ms) | %d train rows\n', ...
            nUnits, numel(lags), nPred, min(winLen), max(winLen), ...
            1000*params.dt*min(winLen), 1000*params.dt*max(winLen), sum(okTr));
        logf('  [target] zero on %.1f%% of window samples\n', 100*mean(yTrain(okTr) == 0));

        % ---------------- FIT: RIDGE PATH ON ALL TRAINING ROWS ----------------
        % One ridge path over all training rows gives the coefficients at every
        % candidate lambda. Lambda is chosen by cfg.lambdaRule on the training trials and
        % the coefficients are read off this path at that lambda. The test set is not
        % used.
        if cfg.ridgeCheckEquivalence && ~exist('ridgeChecked','var')
            checkRidgeEquivalence(Xtrain(okTr,:), yTrain(okTr), cfg.ridgeGrid, cfg);
            ridgeChecked = true;
        end
        R  = ridgePath(Xtrain(okTr,:), yTrain(okTr), cfg.ridgeGrid, cfg);
        Ptr = R.b0 + Xtrain*R.beta;
        Pte = R.b0 + Xtest *R.beta;

        switch lower(cfg.lambdaRule)
            case 'cv5'
                % k-fold CV inside the training set (k = cfg.cvFolds), squared errors pooled
                % across folds. The coefficients at the selected lambda are already in R.
                [kSel, cvSSE] = selectLambdaCV5(Xtrain(okTr,:), yTrain(okTr), rowTrain(okTr), cfg);
                lamHow = sprintf('%d-fold trial-wise CV', cfg.cvFolds);
            otherwise
                error('cfg.lambdaRule must be ''cv5''; got ''%s''.', cfg.lambdaRule);
        end
        b0 = R.b0(kSel);  beta = R.beta(:,kSel);

        ssTr = sum((yTrain(okTr) - mean(yTrain(okTr))).^2);
        ssTe = sum((yTest(okTe)  - mean(yTest(okTe))).^2);
        r2Tr = nan(nK,1);  r2Te = nan(nK,1);
        for kk = 1:nK
            r2Tr(kk) = 1 - sum((yTrain(okTr) - Ptr(okTr,kk)).^2) / ssTr;
            r2Te(kk) = 1 - sum((yTest(okTe)  - Pte(okTe,kk)).^2) / ssTe;
        end
        logf('  [ridge] lambda %.3g (%s, edf %.1f of %d) | train R^2 %.4f | test R^2 %.4f\n', ...
            cfg.ridgeGrid(kSel), lamHow, R.edf(kSel), nPred+1, r2Tr(kSel), r2Te(kSel));

        % ---------------- FULL-TRIAL TEST PREDICTION ----------------
        % The lagged linear model is applied as a filter over each whole test trial
        % (see predictFullTrials), without building the full lagged design matrix.
        predFull = predictFullTrials(spikes(:,:,testTrials), lags, b0, beta, nUnits);

        % ---------------- DECODING INDEX vs the ZERO baseline, per test set --------
        aIdx = find(traceTime >= cfg.analysisWin_s(1) & traceTime <= cfg.analysisWin_s(2));
        [~, aZero] = min(abs(traceTime(aIdx)));

        decIdx  = nan(nLA, nTS);
        peakL   = nan(nLA, nTS);
        nTrialL = zeros(nLA, nTS);

        for ts = 1:nTS
            if isempty(spec.testSets(ts).condIdx)
                sel = true(nTest,1);
            else
                member = unique(cell2mat(S.trialid(spec.testSets(ts).condIdx)'));
                sel = ismember(testTrials, member);
            end
            if ~any(sel), continue; end
            [decIdx(:,ts), peakL(:,ts), nTrialL(:,ts)] = ...
                contactDecodingIndex(tongueLen(aIdx, testTrials(sel)), predFull(aIdx, sel), aZero, cfg);
            logf('  [decIdx %-6s] n=%3d trials |', spec.testSets(ts).name, sum(sel));
            logf(' %6.3f', decIdx(:,ts));  logf('\n');
        end

        % ---------------- STORE ----------------
        nFitted = nFitted + 1;
        k = numel(res) + 1;
        res(k).session = sessionIdx;   res(k).probe = probeIdx;
        res(k).anm = anm;              res(k).date  = dte;
        res(k).region = regionStr;
        res(k).nUnits = nUnits;        res(k).nPred = nPred;
        res(k).nTrain = nTrain;        res(k).nTest = nTest;
        res(k).lambda = cfg.ridgeGrid(kSel);
        res(k).edf    = R.edf(kSel);
        res(k).r2Tr = r2Tr;            res(k).r2Te = r2Te;
        res(k).r2TrSel = r2Tr(kSel);   res(k).r2TeSel = r2Te(kSel);
        res(k).decIdx = decIdx;        res(k).peakL = peakL;   res(k).nTrialL = nTrialL;
        res(k).winLenSamples = winLen;
        res(k).tlTime    = traceTime;   % for the tongue-length panel
        res(k).tlWinMean = tlWinMean;
        res(k).tlWinN    = tlWinN;

        % ---------------- DIAGNOSTIC FIGURES (tabbed, one tab per session) --------
        wantDiag = ismember(sessionIdx, diagList);
        wantWarp = ismember(sessionIdx, warpList);
        if ~wantDiag && ~wantWarp, clear spikes Xtrain Xtest Ptr Pte predFull; continue; end

        % Every title carries the REGION this probe belongs to in this session, because
        % a probe index alone is ambiguous -- probe 2 is ALM in one session and M1 in the
        % next. groupOfProbe reads it off the group maps.
        tabTitle = sprintf('%s | s%d %s %s p%d', regionStr, sessionIdx, anm, dte, probeIdx);
        hdr      = sprintf('%s | %s | sess %d | %s %s | probe %d', spec.name, regionStr, ...
            sessionIdx, anm, dte, probeIdx);

        if wantDiag
            % =========================== FIGURE A: FIT DIAGNOSTICS =====================
            % One window, one TAB PER SESSION -- click along the tab strip to go through
            % them. subplot() only targets a figure, so gridAxes() places each panel inside
            % the tab at the position subplot(2,3,k) would have used.
            if isempty(diagFig) || ~ishandle(diagFig)
                diagFig = figure('Name', sprintf('%s%s - fit diagnostics', passTag, spec.name), ...
                    'Position', [60 70 1400 800]);
                diagTG  = uitabgroup(diagFig);
            end
            tb = uitab(diagTG, 'Title', tabTitle);

            % ---- 1: the lambda path and where the selection rule landed ----
            ax = gridAxes(tb, 2, 3, 1);  hold(ax,'on')
            plot(ax, logK, r2Tr, '--', 'Color', predCol, 'LineWidth',1.5);
            plot(ax, logK, r2Te, '-',  'Color', predCol, 'LineWidth',2);
            gn = cvSSE;  gn(~isfinite(gn)) = NaN;
            gn = (gn - min(gn)) ./ max(range(gn(isfinite(gn))), eps);   % rescaled for shape only
            plot(ax, logK, gn, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
            ylim(ax, [-0.2 1]);
            line(ax, [logK(kSel) logK(kSel)], [-0.2 1], 'Color','k', 'HandleVisibility','off');
            line(ax, [logK(1) logK(end)], [0 0], 'Color','k', 'LineStyle',':', 'HandleVisibility','off');
            xlabel(ax, 'log_{10} \lambda'); ylabel(ax, 'R^2   (CV error rescaled)');
            legend(ax, {'train','test','CV error (scaled)'}, 'Location','southwest','FontSize',7);
            if isempty(cfg.lambdaFixed), lamSrc = sprintf('%d-fold CV on TRAIN', cfg.cvFolds); else, lamSrc = 'FIXED'; end
            title(ax, sprintf('\\lambda by %s = %.3g | edf %.0f of %d', lamSrc, cfg.ridgeGrid(kSel), ...
                R.edf(kSel), nPred+1), 'FontSize',9,'FontWeight','normal');
            box(ax,'off')

            % ---- 2: how long the per-trial windows came out ----
            ax = gridAxes(tb, 2, 3, 2);
            histogram(ax, winLen * params.dt * 1000, 20, 'FaceColor',[0.35 0.35 0.35], 'EdgeColor','none');
            xlabel(ax, 'fitting window length (ms)'); ylabel(ax, 'train trials');
            if strcmpi(cfg.fitWindowMode,'bout')
                winTxt = sprintf('contact 2 +%.0f ms  ->  bout end +%.0f ms', 1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s);
            else
                winTxt = sprintf('go cue -%.0f ms  ->  tongue retracts +%.0f ms', 1000*cfg.preGoCue_s, 1000*cfg.postTongueRetract_s);
            end
            title(ax, ['per-trial window: ' winTxt], 'FontSize',9,'FontWeight','normal');
            box(ax,'off')

            % ---- 3: TRAIN, inside the fitting window, on the real time axis ----
            % The edges of this average rest on a handful of trials, so the x range is
            % trimmed to samples at least covFrac of the trials actually cover.
            ax = gridAxes(tb, 2, 3, 3);  hold(ax,'on')
            Mact = nan(nT, nTrain);  Mpred = nan(nT, nTrain);
            for tt = 1:nTrain
                ki = winTrain{tt};
                Mact(ki,tt)  = tongueLen(ki, trainTrials(tt));
                Mpred(ki,tt) = Ptr(rowTrain == tt, kSel);
            end
            covFrac = 0.2;
            covN    = sum(isfinite(Mact), 2);
            hA = ciBand(ax, traceTime, Mact,  [0 0 0]);
            hP = ciBand(ax, traceTime, Mpred, predCol);
            gd = find(covN >= covFrac*nTrain);
            if ~isempty(gd), xlim(ax, [traceTime(gd(1))-0.02, traceTime(gd(end))+0.02]); end
            ylim(ax, [-0.25 1.05]);  vline0(ax, 'r');
            xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'tongue length (norm)');
            legend(ax, [hA hP], {'actual','predicted'}, 'Location','northwest','FontSize',7);
            title(ax, {'TRAIN, in the fitting window only', ...
                sprintf('>= %.0f%% trial coverage | window mode: %s', 100*covFrac, cfg.fitWindowMode)}, ...
                'FontSize',9,'FontWeight','normal');
            box(ax,'off')

            % ---- 4: the whole trial, TEST ----
            ax = gridAxes(tb, 2, 3, 4);  hold(ax,'on')
            h1 = plot(ax, traceTime, mean(tongueLen(:,testTrials), 2, 'omitnan'), 'k', 'LineWidth', 2);
            h2 = plot(ax, traceTime, mean(predFull, 2, 'omitnan'), 'Color', predCol, 'LineWidth', 2);
            xlim(ax, [-0.5 1.5]);  vline0(ax, 'r');
            xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'tongue length (norm)');
            legend(ax, [h1 h2], {'actual','predicted'}, 'Location','northeast','FontSize',7);
            title(ax, {'full trial, TEST', 'contacts 2+ are extrapolation, not held-out prediction'}, ...
                'FontSize',9,'FontWeight','normal');
            box(ax,'off')

            % ---- 5 and 6: the FIRST and the LAST test trial, one panel each ----
            % Averages hide what a single trial looks like. If the fit is drifting across the
            % session it shows here and nowhere else on this tab.
            showTr = unique([1, numel(testTrials)], 'stable');
            for ii = 1:2
                ax = gridAxes(tb, 2, 3, 4+ii);  hold(ax,'on')
                if ii > numel(showTr), axis(ax,'off'); continue; end
                tt = showTr(ii);
                h1 = plot(ax, traceTime, tongueLen(:, testTrials(tt)), '-', 'Color',[0 0 0], 'LineWidth',1.25);
                h2 = plot(ax, traceTime, predFull(:,tt), '-', 'Color', predCol, 'LineWidth',1.75);
                xlim(ax, [-0.5 1.5]);  vline0(ax, 'r');
                xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'tongue length (norm)');
                legend(ax, [h1 h2], {'actual','predicted'}, 'Location','northeast','FontSize',7);
                if ii == 1, which1 = 'FIRST'; else, which1 = 'LAST'; end
                title(ax, sprintf('%s test trial of the session (trial %d)', which1, testTrials(tt)), ...
                    'FontSize',9,'FontWeight','normal');
                box(ax,'off')
            end

            annotation(tb, 'textbox', [0 0.955 1 0.04], 'String', hdr, 'EdgeColor','none', ...
                'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',11);
        end   % wantDiag

        % =========================== FIGURE B: WARPED LICK TRAIN ===================
        % Each protrusion is resampled to cfg.warpedSegmentLength points and placed in
        % its own slot; mean +/- 1.96 SEM, drawn segment by segment so nothing is joined
        % across the gaps.
        if wantWarp
            if isempty(warpFig) || ~ishandle(warpFig)
                warpFig = figure('Name', sprintf('%s%s - warped lick train', passTag, spec.name), ...
                    'Position', [110 110 1050 620]);
                warpTG  = uitabgroup(warpFig);
            end
            wt = uitab(warpTG, 'Title', tabTitle);

            [Wa, Wp, slotStart] = warpedLickMatrix(tongueLen(aIdx, testTrials), predFull(aIdx,:), aZero, cfg);
            nShow = min(spec.warpShowLicks, cfg.maxLicksToExtract);

            % ---- ONE PANEL PER CONDITION, STACKED ----
            % With two test sets (R1 vs R4, R1 vs R6) the tab holds two axes, R1 on top, in
            % the order spec.testSets declares them. Each panel is one condition: actual
            % black, predicted in that condition's color, labeled in-plot. Only the bottom
            % panel carries the x axis. Both panels are cut from the SAME warped matrix, so
            % the slot geometry is identical and the two are directly comparable.
            if nTS >= 2 && cfg.warpSplitTestSets, panels = 1:nTS; else, panels = 1; end
            nPan = numel(panels);
            nStr = '';
            for pp = 1:nPan
                ts = panels(pp);
                if nPan == 1 || isempty(spec.testSets(ts).condIdx)
                    selTS = true(numel(testTrials),1);
                else
                    memberTS = unique(cell2mat(S.trialid(spec.testSets(ts).condIdx)'));
                    selTS = ismember(testTrials, memberTS);
                end
                axP = gridAxes(wt, nPan, 1, pp);
                if ~any(selTS), axis(axP,'off'); continue; end
                if nPan == 1
                    colP = predCol;  nameStr = 'Model';
                else
                    colP = tsCols(min(ts, size(tsCols,1)), :);  nameStr = spec.testSets(ts).name;
                end
                drawWarpedSlotPanel(axP, Wa(:,selTS), Wp(:,selTS), slotStart, nShow, ...
                    colP, nameStr, pp == nPan, cfg);
                nStr = [nStr sprintf('%s n = %d   ', nameStr, sum(selTS))];   %#ok<AGROW>
            end

            annotation(wt, 'textbox', [0 0.955 1 0.04], 'String', ...
                sprintf('%s  |  warped licks 1-%d, TEST  |  mean \\pm 1.96 SEM  |  %s', hdr, nShow, strtrim(nStr)), ...
                'EdgeColor','none', 'HorizontalAlignment','center', ...
                'FontWeight','bold', 'FontSize',10, 'Interpreter','tex');
        end   % wantWarp

        clear spikes Xtrain Xtest Ptr Pte predFull Mact Mpred Wa Wp
    end   % probe
end   % session

%% REPORT
nLA  = cfg.nLicksAnalyze;
nG   = numel(spec.groupMaps);
nTS  = numel(spec.testSets);
lickVec = (1:nLA)';

% Colors come from the study block:
%   spec.tsColors   one row per test set
%   spec.grpColors  one row per group, or [] for the default gradient
tsCols = spec.tsColors;
if isempty(spec.grpColors)
    grpCols = [linspace(0.10,0.85,max(nG,2))', zeros(max(nG,2),1), linspace(0.85,0.10,max(nG,2))'];
else
    grpCols = spec.grpColors;
end
assert(size(tsCols,1) >= nTS, 'spec.tsColors has %d rows for %d test sets.', size(tsCols,1), nTS);

out.curves  = cell(nG, nTS);
out.peaks   = cell(nG, nTS);
out.keep    = cell(nG, nTS);
out.sessIdx = cell(1, nG);
out.anmName = cell(1, nG);

for g = 1:nG
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    out.sessIdx{g} = [res(rows).session];
    out.anmName{g} = {res(rows).anm};   % animal behind each column, for FIGURE 4
    for ts = 1:nTS
        C = nan(nLA, numel(rows));  P = nan(nLA, numel(rows));
        for cc = 1:numel(rows)
            C(:,cc) = res(rows(cc)).decIdx(:,ts);
            P(:,cc) = res(rows(cc)).peakL(:,ts);
        end
        out.curves{g,ts} = C;
        out.peaks{g,ts}  = P;
    end
    % Every session with a finite contact-1 value is kept.
    for ts = 1:nTS, out.keep{g,ts} = isfinite(out.curves{g,ts}(1,:)); end
end

%% PER-CONTACT TEST: contact L vs contact 1
% Contact L vs contact 1 at every contact from 2 up: Wilcoxon signed rank
% (tail = spec.vsC1Tail), BH-FDR across contacts 2..nLicksAnalyze.
out.pVsC1_sr = nan(nLA, nG, nTS);
out.qVsC1_sr = nan(nLA, nG, nTS);

fprintf('\n%s\n', repmat('=',1,84));
fprintf('CONTACT L vs CONTACT 1 | %s paired signed-rank | %s | %s\n', spec.vsC1Tail, spec.figRef, spec.name);
fprintf('%s\n', repmat('-',1,84));
for g = 1:nG
  if ~spec.vsC1Test, break; end
    for ts = 1:nTS
        C = out.curves{g,ts}(:, out.keep{g,ts});
        if size(C,2) < 2, continue; end
        for L = 2:nLA
            a = C(L,:);  b = C(1,:);
            ok = isfinite(a) & isfinite(b);
            if sum(ok) < 2, continue; end
            switch lower(spec.vsC1Tail)
                case 'both',  out.pVsC1_sr(L,g,ts) = signrank(a(ok), b(ok));
                case 'left',  out.pVsC1_sr(L,g,ts) = signrank(a(ok), b(ok), 'Tail','left');
                case 'right', out.pVsC1_sr(L,g,ts) = signrank(a(ok), b(ok), 'Tail','right');
                otherwise, error('spec.vsC1Tail must be ''both'', ''left'' or ''right''.');
            end
        end
        out.qVsC1_sr(:,g,ts) = bhFDR(out.pVsC1_sr(:,g,ts));
        fprintf('  %-7s %-9s : ', spec.groupLabels{g}, spec.testSets(ts).name);
        for L = 2:nLA
            mk = ' ';
            if isfinite(out.qVsC1_sr(L,g,ts)) && out.qVsC1_sr(L,g,ts) < cfg.alpha, mk = '*'; end
            fprintf('C%d p=%.3g (q=%.3g)%s  ', L, out.pVsC1_sr(L,g,ts), out.qVsC1_sr(L,g,ts), mk);
        end
        fprintf('| n = %d\n', size(C,2));
    end
end
logf('  (* = BH-FDR q < %.2f; the figure stars these same q values)\n', cfg.alpha);
fprintf('%s\n', repmat('=',1,84));

%% BETWEEN TEST SETS, per contact
% First vs second test set at each contact: paired signed rank, two-tailed,
% BH-FDR across contacts.
out.pBetween = nan(nLA, nG);
out.qBetween = nan(nLA, nG);
if nTS >= 2 && spec.betweenTest
    fprintf('\n%s\n', repmat('-',1,84));
    fprintf('%s vs %s at each contact | two-tailed paired signrank\n', ...
        spec.testSets(1).name, spec.testSets(2).name);
    for g = 1:nG
        A = out.curves{g,1};  B = out.curves{g,2};
        k = out.keep{g,1} & out.keep{g,2};
        if sum(k) < 2, continue; end
        A = A(:,k); B = B(:,k);
        for L = 1:nLA
            ok = isfinite(A(L,:)) & isfinite(B(L,:));
            if sum(ok) < 2, continue; end
            out.pBetween(L,g) = signrank(B(L,ok), A(L,ok));
        end
        out.qBetween(:,g) = bhFDR(out.pBetween(:,g));
        fprintf('  %-7s : ', spec.groupLabels{g});
        for L = 1:nLA
            mk = ' ';
            if isfinite(out.qBetween(L,g)) && out.qBetween(L,g) < cfg.alpha, mk = '*'; end
            fprintf('C%d p=%.3g (q=%.3g)%s  ', L, out.pBetween(L,g), out.qBetween(L,g), mk);
        end
        fprintf('| n = %d, m = %d\n', sum(k), nLA);
    end
    fprintf('%s\n', repmat('-',1,84));
end

out.pSummary = nan(1, max(nG, nTS));

%% SESSION-LEVEL SUMMARY, TEST AND SUMMARY FILE
% One value per session: the change in decoding index from contact 1, averaged
% over spec.groupSummaryLicks,
%   value(session) = mean over L of ( decIdx(L) - decIdx(1) ).
% Contact 1 is the within-session reference, so differences in baseline level
% between sessions do not enter.
out.summary = cell(nG, nTS);
for g = 1:nG
    for ts = 1:nTS
        C = out.curves{g,ts};
        if isempty(C), continue; end
        L = spec.groupSummaryLicks;
        L = L(L >= 1 & L <= nLA);
        d = C(L,:) - C(1,:);
        out.summary{g,ts} = mean(d, 1, 'omitnan');   % 1 x nSessions
    end
end

fprintf('\n%s\n', repmat('=',1,84));
fprintf('%s | single session-level value per session\n', spec.figRef);
fprintf('mean change in decoding index from contact 1, averaged over contacts %s\n', ...
    mat2str(spec.groupSummaryLicks));
fprintf('%s\n', repmat('-',1,84));
for g = 1:nG
    for ts = 1:nTS
        v = out.summary{g,ts};
        if isempty(v), continue; end
        v = v(isfinite(v));
        fprintf('  %-8s %-9s : n = %2d   mean %+0.4f\n', ...
            spec.groupLabels{g}, spec.testSets(ts).name, numel(v), mean(v));
    end
end

switch lower(spec.groupTest)

case 'ranksum'
    % No session-level test is run for this option.

case 'crosstask'
    % Fig. 3H / S4F: per-session summary values of the C1 test set in this task
    % against those of the Simple Reward Task, which the Simple Reward script (F1H)
    % writes to spec.summaryFile and this reads from spec.crossTaskFile, so run that
    % script first. One-tailed rank-sum, one test per region, uncorrected.
    if exist(spec.crossTaskFile, 'file')
        ref = load(spec.crossTaskFile);
        fprintf('\n  cross-task reference: %s (%s)\n', ref.figRef, ref.specName);
        for g = 1:nG
            b = out.summary{g,1};  b = b(isfinite(b));   % this task, C1 test set
            if g > numel(ref.summaryByGroup), continue; end
            a = ref.summaryByGroup{g};  a = a(isfinite(a));
            if numel(a) < 2 || numel(b) < 2
                fprintf('  %-8s: test not run, n = %d vs %d.\n', ...
                    spec.groupLabels{g}, numel(a), numel(b));
                continue
            end
            p1 = ranksum(a, b, 'tail', 'left');   % Simple < Double
            out.pSummary(g) = p1;
            fprintf('  %-8s: Simple n = %2d (%+0.4f) vs Double n = %2d (%+0.4f)\n', ...
                spec.groupLabels{g}, numel(a), mean(a), numel(b), mean(b));
            fprintf('            one-tailed rank-sum p = %.4f%s | floor %.4f\n', ...
                p1, repmat('  *', 1, p1 < cfg.alpha), ...
                1 / nchoosek(numel(a)+numel(b), numel(a)));
        end
    else
        fprintf(['\n  !! %s needs the Simple Reward task''s per-session summary and it is\n' ...
                 '     not on disk yet. Run the matching Simple Reward script first; it\n' ...
                 '     writes:\n       %s\n'], spec.figRef, spec.crossTaskFile);
    end

otherwise
    fprintf('\n  (%s reports a per-contact test; see the table above.)\n', spec.figRef);
end
fprintf('%s\n', repmat('=',1,84));

% ---- save this script's per-session summary (read by the cross-task test) ----
summaryByGroup = cell(1, nG);
for g = 1:nG, summaryByGroup{g} = out.summary{g,1}; end
figRef   = spec.figRef;   %#ok<NASGU>
specName = spec.name;   %#ok<NASGU>
sessionsByGroup = out.sessIdx;   %#ok<NASGU>
groupLabels  = spec.groupLabels;        %#ok<NASGU>   % lets a reader check the regions match
summaryLicks = spec.groupSummaryLicks;  %#ok<NASGU>   % contacts averaged into each value
save(spec.summaryFile, 'summaryByGroup', 'sessionsByGroup', 'groupLabels', ...
     'summaryLicks', 'figRef', 'specName');
logf('summary written to %s\n\n', spec.summaryFile);

%% FIGURE 1: DECODING INDEX, TABBED
if spec.showFigCI
    % One figure per group with two tabs: decoding index (mean +/- 95% CI) and mean
    % tongue length inside the fitting window.
    if spec.overlayGroups
        figs = {1:nG};  figNames = {'all groups'};
    else
        figs = num2cell(1:nG);  figNames = spec.groupLabels;
    end

    for fi = 1:numel(figs)
        gl = figs{fi};
        if all(cellfun(@(c) isempty(c), out.curves(gl,1))), continue; end

        fh = figure('Name', sprintf('%s%s - %s - decoding index (zero baseline)', passTag, spec.name, figNames{fi}));
        tg = uitabgroup(fh);

        % ---- TAB: mean with error bars ----
        ax = axes('Parent', uitab(tg, 'Title', 'mean +/- 95% CI'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = gl
            for ts = 1:nTS
                C = out.curves{g,ts};  k = out.keep{g,ts};
                if isempty(C), continue; end
                col = pickColor(spec, g, ts, grpCols, tsCols);
                [m, ci] = meanCI(C(:,k));
                gd = isfinite(m) & isfinite(ci);
                if any(gd)
                    fill(ax, [lickVec(gd); flipud(lickVec(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                        col, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off');
                end
                h(end+1) = plot(ax, lickVec, m, '-o', 'Color', col, 'MarkerFaceColor', col, ...
                    'LineWidth', 2.5);   %#ok<SAGROW>
                lbl{end+1} = seriesLabel(spec, g, ts, sum(k));   %#ok<SAGROW>
            end
        end
        finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{zero}', h, lbl, [0 1]);

        % stars: contact L vs contact 1, BH-corrected q (see drawStars)
        if numel(gl) == 1
            if spec.vsC1Test
                drawStars(ax, 2:nLA, [], out.qVsC1_sr(:,gl,1), cfg.alpha);
            elseif spec.betweenTest
                drawStars(ax, 1:nLA, [], out.qBetween(:,gl),   cfg.alpha);
            end
        end
        ttl = sprintf('%s | %s | zero baseline | mean \\pm 95%% CI', spec.name, figNames{fi});
        if spec.vsC1Test
            sub = sprintf('* = BH-FDR q < %.2f, %s signed-rank, contact L vs contact 1', cfg.alpha, spec.vsC1Tail);
        elseif spec.betweenTest
            sub = sprintf('* = BH-FDR q < %.2f, two-tailed signed-rank, %s vs %s', ...
                cfg.alpha, spec.testSets(1).name, spec.testSets(2).name);
        else
            sub = '';
        end
        title(ax, {ttl, sub}, 'FontSize',9, 'FontWeight','normal');

        % ---- TAB: mean tongue length inside the fitting window ----
        % Each trial contributes only the samples inside its own window, so each sample
        % is averaged over the trials that cover it; that count is on the right axis.
        ax = axes('Parent', uitab(tg, 'Title', 'tongue length in window'));  hold(ax,'on')
        rows = [];
        for g = gl
            rows = [rows, find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res))];   %#ok<AGROW>
        end
        if ~isempty(rows)
            tt = res(rows(1)).tlTime;
            M  = nan(numel(tt), numel(rows));
            N  = nan(numel(tt), numel(rows));
            for cc = 1:numel(rows)
                M(:,cc) = res(rows(cc)).tlWinMean;
                N(:,cc) = res(rows(cc)).tlWinN;
            end
            [m, ci] = meanCI(M);
            gd = isfinite(m) & isfinite(ci);
            yyaxis(ax,'left')
            if any(gd)
                fill(ax, [tt(gd); flipud(tt(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                    [0 0 0], 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
            end
            plot(ax, tt, m, '-', 'Color', [0 0 0], 'LineWidth', 2.5);
            ylabel(ax, 'tongue length (norm)', 'FontSize', 11);  ylim(ax, [0 1]);
            yyaxis(ax,'right')
            plot(ax, tt, mean(N, 2, 'omitnan'), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.25);
            ylabel(ax, 'trials covering the sample', 'FontSize', 10);
            yyaxis(ax,'left')
            good = find(mean(N,2,'omitnan') > 0);
            if ~isempty(good)
                xlim(ax, [tt(good(1))-0.05, tt(good(end))+0.05]);
            end
        end
        plot(ax, [0 0], [0 1], 'r--', 'HandleVisibility','off');
        xlabel(ax, 'time from contact 1 (s)', 'FontSize', 11);  box(ax,'off')
        title(ax, {sprintf('%s | %s | mean tongue length inside the fitting window', spec.name, figNames{fi}), ...
            'mean \pm 95% CI across sessions; red line = contact 1 (t = 0)'}, ...
            'FontSize',9,'FontWeight','normal');
    end
end

%% FIGURE 2: MEANS ONLY (no error bars)
% The means only, all groups on one axes (drawn when spec.overlayGroups is true).
if spec.overlayGroups
    fh = figure('Name', sprintf('%s%s - mean traces only', passTag, spec.name));
    ax = axes('Parent', fh); hold(ax,'on')
    [h, lbl] = deal(gobjects(0), {});
    for g = 1:nG
        for ts = 1:nTS
            C = out.curves{g,ts};  k = out.keep{g,ts};
            if isempty(C), continue; end
            col = pickColor(spec, g, ts, grpCols, tsCols);
            h(end+1) = plot(ax, lickVec, mean(C(:,k),2,'omitnan'), '-o', 'Color', col, ...
                'MarkerFaceColor', col, 'LineWidth', 2.5);   %#ok<SAGROW>
            lbl{end+1} = seriesLabel(spec, g, ts, sum(k));   %#ok<SAGROW>
        end
    end
    finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{zero}', h, lbl, [0 1]);
    title(ax, sprintf('%s | mean decoding index', spec.name), ...
        'FontSize',10,'FontWeight','normal');
end

%% FIGURE 3: DELTA FROM CONTACT 1, TABBED
if spec.showFigDelta
    % Every session is 0 at contact 1 by construction; this is the quantity the
    % contact-L-vs-contact-1 test is run on. One figure per group unless
    % spec.overlayGroups is true.
    if isempty(spec.deltaGroups), dgAll = 1:nG; else, dgAll = spec.deltaGroups; end
    dgAll = dgAll(dgAll >= 1 & dgAll <= nG);
    if spec.overlayGroups, dFigs = {dgAll}; dNames = {'compared groups'};
    else,                  dFigs = num2cell(dgAll); dNames = spec.groupLabels(dgAll);
    end

    for di = 1:numel(dFigs)
        dg = dFigs{di};
        if all(cellfun(@isempty, out.curves(dg,1))), continue; end

        fh = figure('Name', sprintf('%s%s - %s - delta from contact 1', passTag, spec.name, dNames{di}));
        tg = uitabgroup(fh);

        ax = axes('Parent', uitab(tg, 'Title', 'mean +/- 95% CI'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = dg
            for ts = 1:nTS
                C = out.curves{g,ts};  k = out.keep{g,ts};
                if isempty(C), continue; end
                D = C(:,k) - C(1,k);
                col = pickColor(spec, g, ts, grpCols, tsCols);
                [m, ci] = meanCI(D);
                gd = isfinite(m) & isfinite(ci);
                if any(gd)
                    fill(ax, [lickVec(gd); flipud(lickVec(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                        col, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off');
                end
                h(end+1) = plot(ax, lickVec, m, '-o', 'Color', col, 'MarkerFaceColor', col, ...
                    'LineWidth', 2.5);   %#ok<SAGROW>
                lbl{end+1} = seriesLabel(spec, g, ts, sum(k));   %#ok<SAGROW>
            end
        end
        finishAxes(ax, nLA, '\Delta decoding index (from contact 1)', h, lbl, []);
        if numel(dg) == 1
            drawStars(ax, 2:nLA, [], out.qVsC1_sr(:,dg,1), cfg.alpha);
        elseif ~isempty(spec.compareGroups) && all(ismember(spec.compareGroups, dg))
        end
        title(ax, {sprintf('%s | %s | delta from contact 1', spec.name, dNames{di}), ...
            'contact 1 is 0 for every session by construction'}, ...
            'FontSize',9,'FontWeight','normal');

        ax = axes('Parent', uitab(tg, 'Title', 'per session'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = dg
            C = out.curves{g,1};  k = out.keep{g,1};
            if isempty(C), continue; end
            D = C - C(1,:);
            col = pickColor(spec, g, 1, grpCols, tsCols);
            for cc = 1:size(D,2)
                if sum(isfinite(D(:,cc))) < 2 || ~k(cc), continue; end
                hs = plot(ax, lickVec, D(:,cc), '-', 'Color', [col 0.35], 'LineWidth', 1, ...
                    'HandleVisibility','off');
                hs.DisplayName = sprintf('sess %d', out.sessIdx{g}(cc));
            end
            h(end+1) = plot(ax, lickVec, mean(D(:,k),2,'omitnan'), '-o', 'Color', col, ...
                'MarkerFaceColor', col, 'LineWidth', 3);   %#ok<SAGROW>
            lbl{end+1} = sprintf('%s mean (n = %d)', spec.groupLabels{g}, sum(k));   %#ok<SAGROW>
        end
        finishAxes(ax, nLA, '\Delta decoding index (from contact 1)', h, lbl, []);
        title(ax, sprintf('%s | %s | delta, every session', spec.name, dNames{di}), ...
            'FontSize',9,'FontWeight','normal');
    end

%% FIGURE 4: PER ANIMAL
    % Sessions are averaged within each animal first, then the mean and 95% CI are
    % taken across animals, so the animal is the unit of replication. One thin line
    % per animal per test set; the bold line is the mean across animals.
    % ---- controls for this figure only ----
    fig4.showCI    = true;   % shaded 95% CI across animals behind the bold mean.
    fig4.lineAlpha = 0.35;   % transparency of the individual animal traces
    fig4.lineWidth = 1.25;   % width of the individual animal traces
    fig4.yLim      = [0 1];   % [] lets MATLAB choose

    out.anmNames  = cell(nG, nTS);
    out.anmCurves = cell(nG, nTS);
    out.anmNSess  = cell(nG, nTS);

    for g = 1:nG
        nmAll = out.anmName{g};
        if isempty(nmAll), continue; end
        for ts = 1:nTS
            C = out.curves{g,ts};  k = out.keep{g,ts};
            if isempty(C) || ~any(k), continue; end
            [an, M, ns] = collapseByAnimal(nmAll(k), C(:,k));
            out.anmNames{g,ts}  = an;
            out.anmCurves{g,ts} = M;
            out.anmNSess{g,ts}  = ns;
        end
    end

    % ---- console table ----
    fprintf('\n%s\n', repmat('=',1,84));
    fprintf('PER-ANIMAL decoding index | sessions averaged within animal | %s\n', spec.name);
    fprintf('%s\n', repmat('-',1,84));
    for g = 1:nG
        for ts = 1:nTS
            M = out.anmCurves{g,ts};
            if isempty(M), continue; end
            fprintf('  %-7s %-6s : %d animals, %d sessions\n', ...
                spec.groupLabels{g}, spec.testSets(ts).name, size(M,2), sum(out.anmNSess{g,ts}));
            for i = 1:size(M,2)
                fprintf('    %-8s (%d sess) :', out.anmNames{g,ts}{i}, out.anmNSess{g,ts}(i));
                fprintf(' %6.3f', M(:,i));  fprintf('\n');
            end
            fprintf('    %-8s          :', 'MEAN');
            fprintf(' %6.3f', mean(M, 2, 'omitnan'));  fprintf('\n');

        end
    end
    fprintf('%s\n', repmat('=',1,84));

    % ---- the figures: one per group, same rule as everywhere else ----
    if spec.overlayGroups
        aFigs = {1:nG};  aNames = {'all groups'};
    else
        aFigs = num2cell(1:nG);  aNames = spec.groupLabels;
    end

    for fi = 1:numel(aFigs)
        gl = aFigs{fi};
        if all(cellfun(@isempty, out.anmCurves(gl,1))), continue; end

        fh = figure('Name', sprintf('%s%s - %s - PER ANIMAL', passTag, spec.name, aNames{fi}));
        tg = uitabgroup(fh);

        % ---- TAB: decoding index, mean +/- 95% CI across animals ----
        ax = axes('Parent', uitab(tg, 'Title', 'per animal'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = gl
            for ts = 1:nTS
                M = out.anmCurves{g,ts};
                if isempty(M), continue; end
                col = pickColor(spec, g, ts, grpCols, tsCols);

                % every animal as a thin line, named so a datatip identifies it
                for i = 1:size(M,2)
                    hs = plot(ax, lickVec, M(:,i), '-', 'Color', [col fig4.lineAlpha], ...
                        'LineWidth', fig4.lineWidth, 'HandleVisibility','off');
                    hs.DisplayName = sprintf('%s %s (%d sess)', out.anmNames{g,ts}{i}, ...
                        spec.testSets(ts).name, out.anmNSess{g,ts}(i));
                end

                [m, ci] = meanCI(M);
                if fig4.showCI
                    gd = isfinite(m) & isfinite(ci);
                    if any(gd)
                        fill(ax, [lickVec(gd); flipud(lickVec(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                            col, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off');
                    end
                end
                h(end+1) = plot(ax, lickVec, m, '-o', 'Color', col, 'MarkerFaceColor', col, ...
                    'LineWidth', 3);   %#ok<SAGROW>
                lbl{end+1} = sprintf('%s (%d animals, %d sess)', spec.testSets(ts).name, ...
                    size(M,2), sum(out.anmNSess{g,ts}));   %#ok<SAGROW>
            end
        end
        finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{zero}', h, lbl, fig4.yLim);
        title(ax, {sprintf('%s | %s | PER ANIMAL', spec.name, aNames{fi}), ...
            'thin line = one animal per test set (its sessions averaged); bold = mean across animals'}, ...
            'FontSize',9,'FontWeight','normal');

        % ---- TAB: delta from contact 1, per animal ----
        ax = axes('Parent', uitab(tg, 'Title', 'delta from contact 1'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = gl
            for ts = 1:nTS
                M = out.anmCurves{g,ts};
                if isempty(M), continue; end
                D = M - M(1,:);   % every animal starts at 0
                col = pickColor(spec, g, ts, grpCols, tsCols);

                for i = 1:size(D,2)
                    hs = plot(ax, lickVec, D(:,i), '-', 'Color', [col fig4.lineAlpha], ...
                        'LineWidth', fig4.lineWidth, 'HandleVisibility','off');
                    hs.DisplayName = sprintf('%s %s', out.anmNames{g,ts}{i}, spec.testSets(ts).name);
                end

                [m, ci] = meanCI(D);
                if fig4.showCI
                    gd = isfinite(m) & isfinite(ci);
                    if any(gd)
                        fill(ax, [lickVec(gd); flipud(lickVec(gd))], [m(gd)+ci(gd); flipud(m(gd)-ci(gd))], ...
                            col, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off');
                    end
                end
                h(end+1) = plot(ax, lickVec, m, '-o', 'Color', col, 'MarkerFaceColor', col, ...
                    'LineWidth', 3);   %#ok<SAGROW>
                lbl{end+1} = sprintf('%s (%d animals)', spec.testSets(ts).name, size(D,2));   %#ok<SAGROW>
            end
        end
        finishAxes(ax, nLA, '\Delta decoding index (from contact 1)', h, lbl, []);
        title(ax, {sprintf('%s | %s | PER ANIMAL, delta from contact 1', spec.name, aNames{fi}), ...
            'contact 1 is 0 for every animal by construction'}, ...
            'FontSize',9,'FontWeight','normal');
    end
end

%% FIT QUALITY TABLE
logf('\n%s\n', repmat('=',1,101));
logf('%-5s %-8s %-12s %5s %-8s %6s %6s %6s %10s %7s %9s %9s %11s\n', ...
    'sess','animal','date','probe','region','units','nTrn','nTst','lambda','edf', ...
    'R2 train','R2 test','win samp');
logf('%s\n', repmat('-',1,101));
for r = res
    wl = r.winLenSamples;
    logf('%-5d %-8s %-12s %5d %-8s %6d %6d %6d %10.3g %7.1f %9.3f %9.3f %5d-%-5d\n', ...
        r.session, r.anm, r.date, r.probe, r.region, r.nUnits, r.nTrain, r.nTest, ...
        r.lambda, r.edf, r.r2TrSel, r.r2TeSel, min(wl), max(wl));
end
logf('%s\n', repmat('=',1,101));

end   % repIdx -- REPEAT PASSES

%% LOCAL FUNCTIONS
function r = tongueRuns(v, cfg)
    % Protrusions in one trial's tongue length: [onset offset] sample pairs, one row
    % per protrusion (a contiguous run of length > 0). Runs separated by
    % <= cfg.mergeGapSamples zeros are merged first (tracking dropouts); only then
    % are runs shorter than cfg.minContactSamples removed, so the two halves of a
    % dropout-split protrusion are rejoined before the duration rule applies.
    v = v(:);
    above = v > 0;
    e = diff([0; above; 0]);
    r = [find(e == 1), find(e == -1) - 1];
    if isempty(r), return; end

    if cfg.mergeGapSamples > 0 && size(r,1) > 1
        keepRow = true(size(r,1),1);
        for ii = 2:size(r,1)
            prev = find(keepRow(1:ii-1), 1, 'last');
            if r(ii,1) - r(prev,2) - 1 <= cfg.mergeGapSamples
                r(prev,2) = r(ii,2);   % absorb into the previous run
                keepRow(ii) = false;
            end
        end
        r = r(keepRow, :);
    end

    r = r(r(:,2) - r(:,1) + 1 >= cfg.minContactSamples, :);
end

function first = pickFirstRun(r, alignSample)
    % Which protrusion is contact 1: the one containing alignSample (t = 0 is the
    % first post-go-cue port contact, so the tongue is out at that instant). If
    % tracking dropped out across t = 0, fall back to a run that ended up to 3
    % samples earlier, then to the run with the nearest onset.
    first = find(r(:,1) <= alignSample & r(:,2) >= alignSample, 1);
    if ~isempty(first), return; end
    before = find(r(:,2) < alignSample, 1, 'last');   % last one that ended before
    if ~isempty(before) && alignSample - r(before,2) <= 3
        first = before;  return   % dropout right at contact
    end
    [~, first] = min(abs(r(:,1) - alignSample));
end

function [nRun, nInWin] = longestLickRun(lickTimes, goCue, cfg)
    % Longest run of CONSECUTIVE port contacts that are all inside cfg.boutWin_s of
    % the go cue and no more than cfg.boutMaxILI_s apart. nInWin is how many contacts
    % fall in the window at all, printed so a failing session can be diagnosed
    % ("no bouts" vs "bouts with one long gap").
    nRun = 0;  nInWin = 0;
    if isempty(lickTimes), return; end
    t = sort(lickTimes(lickTimes > goCue));
    t = t(t - goCue <= cfg.boutWin_s);
    nInWin = numel(t);
    if nInWin == 0, return; end
    nRun = 1;  cur = 1;
    for ii = 2:numel(t)
        if t(ii) - t(ii-1) <= cfg.boutMaxILI_s, cur = cur + 1; else, cur = 1; end
        nRun = max(nRun, cur);
    end
end

function winIdx = fitWindowForTrial(tlTrial, traceTime, goCueRel_s, contactRel, cfg)
    % Sample indices of one trial's fitting window (see cfg.fitWindowMode):
    %   'contact1'        go cue - cfg.preGoCue_s  ->  retraction of the protrusion
    %                     containing contact 1 + cfg.postTongueRetract_s
    %   'bout'            contact 2 + cfg.winBoutStartPad_s  ->  bout-ending contact
    %                     + cfg.winBoutEndPad_s (contact times only)
    %   'lick2ToBoutEnd'  retraction of lick 2 + cfg.winBoutStartPad_s  ->  same end
    % contactRel holds the post-go-cue contact times relative to contact 1. Returns
    % [] when there is no valid window or it is shorter than cfg.minWindowSamples.
    winIdx = [];
    if ~isfinite(goCueRel_s), return; end

    switch lower(cfg.fitWindowMode)
        case 'contact1'
            r = tongueRuns(tlTrial, cfg);
            if isempty(r), return; end
            [~, alignSample] = min(abs(traceTime));
            first  = pickFirstRun(r, alignSample);
            tStart = goCueRel_s - cfg.preGoCue_s;
            tEnd   = traceTime(r(first,2)) + cfg.postTongueRetract_s;

        case 'bout'
            if numel(contactRel) < 2, return; end
            kEnd = boutEndContact(contactRel, cfg.winBoutGap_s);
            if kEnd < 2, return; end   % the bout is over before contact 2
            tStart = contactRel(2)    + cfg.winBoutStartPad_s;
            tEnd   = contactRel(kEnd) + cfg.winBoutEndPad_s;

        case 'lick2toboutend'
            % The start is the retraction of the protrusion carrying contact 2 (found
            % with the same containment rule as contact 1, pickFirstRun), plus the pad.
            if numel(contactRel) < 2, return; end
            kEnd = boutEndContact(contactRel, cfg.winBoutGap_s);
            if kEnd < 2, return; end   % the bout is over before contact 2
            r = tongueRuns(tlTrial, cfg);
            if isempty(r), return; end
            [~, s2] = min(abs(traceTime - contactRel(2)));   % sample at contact 2
            i2 = pickFirstRun(r, s2);   % its protrusion
            tStart = traceTime(r(i2,2)) + cfg.winBoutStartPad_s;   % its retraction
            tEnd   = contactRel(kEnd)   + cfg.winBoutEndPad_s;

        otherwise
            error('cfg.fitWindowMode must be ''contact1'', ''bout'' or ''lick2ToBoutEnd''.');
    end

    if tEnd <= tStart, return; end
    w = find(traceTime >= tStart & traceTime <= tEnd);
    if numel(w) >= cfg.minWindowSamples, winIdx = w; end
end

function k = boutEndContact(contactRel, gap_s)
    % Which contact ENDS the bout: scanning forward from contact 1, the first one
    % whose gap to the NEXT contact exceeds gap_s. If no gap exceeds it, the bout
    % runs to the last contact of the trial.
    d = diff(contactRel(:));
    k = find(d > gap_s, 1);
    if isempty(k), k = numel(contactRel); end
end

function s = groupOfProbe(spec, sessionIdx, probeIdx)
    % Group label (region, or learning day) this probe carries in this session, read
    % off the group maps. Overlapping maps are joined with '+'; 'no group' if none
    % points here.
    hit = cellfun(@(m) m(sessionIdx) == probeIdx, spec.groupMaps);
    if ~any(hit)
        s = 'no group';
    else
        s = strjoin(spec.groupLabels(hit), '+');
    end
end

function [anm, M, nSess] = collapseByAnimal(names, C)
    % Average a group's per-session curves within each animal (used by FIGURE 4).
    %   names : 1 x nSess cell of animal names, one per column of C
    %   C     : nContacts x nSess
    %   M     : nContacts x nAnimals, animals in first-appearance order
    %   nSess : how many sessions went into each animal's column
    anm = {};  M = [];  nSess = [];
    if isempty(names), return; end
    names = names(:);
    [anm, ~, ic] = unique(names, 'stable');
    M = nan(size(C,1), numel(anm));
    nSess = zeros(1, numel(anm));
    for i = 1:numel(anm)
        cols = (ic == i);
        M(:,i)   = mean(C(:, cols), 2, 'omitnan');
        nSess(i) = sum(cols);
    end
    anm = anm(:)';
end

function [trainTrials, testTrials] = splitTypeStratified(poolTrials, trialTypes, trainFrac)
    % Equal numbers per trial type up to the quota, remainder held out.
    ty = trialTypes(poolTrials);
    present = unique(ty(isfinite(ty) & ty > 0));
    lists = arrayfun(@(t) poolTrials(ty == t), present, 'UniformOutput', false);
    nPer  = min([cellfun(@numel, lists(:))', floor(numel(poolTrials)*trainFrac/max(numel(present),1))]);
    if nPer < 1, trainTrials = []; testTrials = []; return; end
    picked = cellfun(@(L) randsample(L, nPer, false), lists, 'UniformOutput', false);
    trainTrials = vertcat(picked{:});
    trainTrials = trainTrials(randperm(numel(trainTrials)));
    testTrials  = poolTrials(~ismember(poolTrials, trainTrials));
end

function v = stackWindows(M, winCell, cols)
    % Stack M(winCell{t}, cols(t)) for every t, in the same row order laggedDesign
    % produces -- so X and y are paired by construction rather than by convention.
    lens = cellfun(@numel, winCell(:));
    v = nan(sum(lens), 1);
    r0 = 0;
    for t = 1:numel(winCell)
        ki = winCell{t}(:);
        if isempty(ki), continue; end
        v(r0 + (1:numel(ki))) = M(ki, cols(t));
        r0 = r0 + numel(ki);
    end
end

function [Xw, rowTrial] = laggedDesign(A, lags, winCell)
    % Lagged design from A (nT x nUnits x nTrials) over per-trial sample sets.
    % ROW ORDER: trials in order, samples within a trial in order (as stackWindows).
    % COLUMN ORDER: lag-major, unit-minor; column (li-1)*nUnits + u is unit u at
    % lag lags(li). No lag crosses a trial boundary: a lag that leaves the trial
    % gives NaN and that row is dropped before fitting.
    [nT, nU, nTr] = size(A);
    nL = numel(lags);
    lens = cellfun(@numel, winCell(:));
    Xw = nan(sum(lens), nU*nL);
    rowTrial = zeros(sum(lens), 1);
    r0 = 0;
    for t = 1:nTr
        ki = winCell{t}(:);
        n = numel(ki);
        if n == 0, continue; end
        rows = r0 + (1:n);
        rowTrial(rows) = t;
        for li = 1:nL
            src = ki + lags(li);
            good = src >= 1 & src <= nT;
            blk = nan(n, nU);
            if any(good), blk(good,:) = A(src(good), :, t); end
            Xw(rows, (li-1)*nU + (1:nU)) = blk;
        end
        r0 = r0 + n;
    end
end

function P = predictFullTrials(A, lags, b0, beta, nU)
    % Full-trial prediction for every trial in A (nT x nUnits x nTrials), computed
    % as a filter rather than by materialising a lagged design. Numerically
    % identical to  b0 + laggedDesign(...)*beta, just without the memory traffic.
    [nT, ~, nTr] = size(A);
    P = zeros(nT, nTr);
    valid = true(nT, 1);
    for li = 1:numel(lags)
        bl  = double(beta((li-1)*nU + (1:nU)));
        src = (1:nT)' + lags(li);
        good = src >= 1 & src <= nT;
        valid = valid & good;
        for t = 1:nTr
            P(good,t) = P(good,t) + double(A(src(good), :, t)) * bl;
        end
    end
    P = b0 + P;
    P(~valid, :) = NaN;
end

function R = ridgePath(X, y, lambdas, cfg, wantEdf)
    % Ridge coefficients at every lambda on the grid.
    %   cfg.ridgeImpl = 'matlab'  MATLAB's ridge(y, X, lambdas, 0). The trailing 0
    %                             returns coefficients on the original predictor
    %                             scale, (p+1) x nL with the intercept in row 1.
    %   cfg.ridgeImpl = 'eig'     the algebraically identical closed form.
    % Both penalize the standardized coefficients and leave the intercept
    % unpenalized, so 'matlab' requires cfg.standardizeCols = true.
    % R.edf(k) = 1 + sum_i d_i/(d_i + lambda_k) over the eigenvalues d of the
    % standardized Gram matrix (trace of the hat matrix, +1 for the intercept). Only
    % the fit printout uses it; CV folds pass wantEdf = false.
    if nargin < 5 || isempty(wantEdf), wantEdf = true; end
    standardizeCols = cfg.standardizeCols;
    impl = lower(cfg.ridgeImpl);

    X = double(X);  y = double(y(:));
    lambdas = lambdas(:)';
    nL = numel(lambdas);
    dGram = [];   % set by the 'eig' path; makes R.edf free there

    switch impl
        case 'matlab'
            assert(standardizeCols, ...
                ['cfg.ridgeImpl = ''matlab'' calls the MATLAB function ridge, which ' ...
                 'always penalizes standardized coefficients. standardizeCols = false ' ...
                 'has no equivalent there -- set cfg.standardizeCols = true, or use ' ...
                 'cfg.ridgeImpl = ''eig''.']);
            B      = ridge(y, X, lambdas, 0);   % (p+1) x nL, ORIGINAL scale
            R.b0   = B(1,:);
            R.beta = B(2:end,:);

        case 'eig'
            % One eigendecomposition of the Gram matrix gives the whole path, instead of one
            % augmented least-squares solve per lambda as in ridge().
            mx = mean(X,1);  my = mean(y);
            Xc = X - mx;     yc = y - my;
            if standardizeCols
                % SDs below sqrt(eps) are replaced by 1, the same threshold ridge() uses.
                sx = std(X,0,1);
                sx(abs(sx) < sqrt(eps(class(sx))) | ~isfinite(sx)) = 1;
                Xc = Xc ./ sx;
            else
                sx = ones(1, size(X,2));
            end
            G = Xc'*Xc;  G = (G + G')/2;   % symmetry against round-off
            [V, d] = eig(G, 'vector');
            d  = max(real(d), 0);   % PSD by construction
            Vc = V' * (Xc'*yc);
            clear Xc   % free the p-column copy now
            % THE WHOLE PATH AS ONE MATRIX PRODUCT. Vc ./ (d + lambdas) is
            % p x nL by implicit expansion, so this replaces nL separate
            % p x p matrix-vector products with a single BLAS-3 call.
            R.beta = (V * (Vc ./ (d + lambdas))) ./ sx(:);
            R.b0   = my - mx*R.beta;
            % Keep the eigenvalues for R.edf (avoids a separate svd).
            dGram  = d;

        otherwise
            error('cfg.ridgeImpl must be ''matlab'' or ''eig''; got ''%s''.', impl);
    end

    if wantEdf
        if ~isempty(dGram)
            % The Gram eigenvalues are the squared singular values of the centered, scaled
            % design.
            d = dGram;
        else
            mx = mean(X,1);
            sx = std(X,0,1);
            sx(abs(sx) < sqrt(eps(class(sx))) | ~isfinite(sx)) = 1;
            d  = svd((X - mx) ./ sx, 'econ').^2;
        end
        R.edf = 1 + sum(d(:) ./ (d(:) + lambdas), 1);
    else
        R.edf = nan(1, nL);
    end
end

function checkRidgeEquivalence(X, y, lambdas, cfg)
    % Fit the same data with ridge() and the closed form and report the largest
    % relative coefficient difference; warns if it exceeds 1e-6.
    p = size(X,2);
    if p > 2000
        logf(['  [ridge check] skipped: %d predictors would mean %d MATLAB ridge\n' ...
                 '                solves. Check on a smaller session.\n'], p, numel(lambdas));
        return
    end
    Bm = ridge(y, X, lambdas, 0);
    mx = mean(X,1);  my = mean(y);
    sx = std(X,0,1); sx(sx == 0 | ~isfinite(sx)) = 1;
    Xc = (X - mx) ./ sx;  yc = y - my;
    G  = Xc'*Xc;  G = (G+G')/2;
    [V,d] = eig(G,'vector');  d = max(real(d),0);  Vc = V'*(Xc'*yc);
    Be = zeros(p+1, numel(lambdas));
    for k = 1:numel(lambdas)
        bRaw = (V*(Vc./(d+lambdas(k))))./sx(:);
        Be(:,k) = [my - mx*bRaw; bRaw];
    end
    rel = max(abs(Bm(:) - Be(:))) / max(abs(Bm(:)));
    logf('  [ridge check] MATLAB ridge vs closed form: max relative difference %.3g\n', rel);
    if rel > 1e-6
        warning(['MATLAB ridge and the closed form differ by %.3g. Use ' ...
                 'cfg.ridgeImpl = ''matlab''.'], rel);
    end
end

function [Wa, Wp, slotStart, traceLen] = warpedLickMatrix(Yact, Ypred, alignSample, cfg)
    % Every protrusion is resampled onto cfg.warpedSegmentLength points and placed in
    % its own slot on a synthetic axis, slotStart = warpedSlotFirst :
    % warpedSlotSpacing : ... . The spacing exceeds the segment length, so the gaps
    % between slots stay NaN. Resampling uses the finite samples of each segment
    % (cfg.warpedInterp); a segment with one finite sample is placed at the slot
    % center. Predictions are warped through the same sample indices as the actual
    % trace. Lick 1 is found with tongueRuns / pickFirstRun.
    nEx  = cfg.maxLicksToExtract;
    segL = cfg.warpedSegmentLength;
    slotStart = cfg.warpedSlotFirst : cfg.warpedSlotSpacing : ...
                (cfg.warpedSlotFirst + cfg.warpedSlotSpacing*(nEx-1));
    traceLen  = max(slotStart) + segL - 1;

    nTr = size(Yact, 2);
    Wa  = nan(traceLen, nTr);
    Wp  = nan(traceLen, nTr);

    for tt = 1:nTr
        r = tongueRuns(Yact(:,tt), cfg);
        if isempty(r), continue; end
        first = pickFirstRun(r, alignSample);
        r = r(first : min(size(r,1), first + nEx - 1), :);
        for L = 1:size(r,1)
            seg  = r(L,1):r(L,2);
            dest = slotStart(L) : slotStart(L) + segL - 1;
            Wa(dest,tt) = warpSegment(Yact(seg,tt),  segL, cfg.warpedInterp);
            Wp(dest,tt) = warpSegment(Ypred(seg,tt), segL, cfg.warpedInterp);
        end
    end
end

function wt = warpSegment(v, segL, method)
    % One excursion resampled onto segL points; a single finite sample goes to the
    % center.
    wt = nan(segL,1);
    xv = v(isfinite(v));
    if numel(xv) == 1
        wt(ceil(segL/2)) = xv;
    elseif numel(xv) > 1
        wt = interp1(linspace(0,1,numel(xv))', xv(:), linspace(0,1,segL)', method, NaN);
    end
end

function drawWarpedSlotPanel(ax, Wa, Wp, slotStart, nShow, colP, nameStr, showX, cfg)
    % One condition on one axes: actual in black, predicted in colP, labeled with
    % text on the axes rather than a legend.
    hold(ax, 'on')
    drawWarpedSet(ax, Wa, [0 0 0], slotStart, cfg, nShow);
    drawWarpedSet(ax, Wp, colP,    slotStart, cfg, nShow);

    if cfg.warpedTickCenter, tickAt = slotStart + cfg.warpedSegmentLength/2; else, tickAt = slotStart; end
    xr = [slotStart(1)-5, slotStart(nShow)+cfg.warpedSegmentLength+5];
    xlim(ax, xr);
    xticks(ax, tickAt(1:nShow));
    if showX
        xticklabels(ax, arrayfun(@(L) sprintf('%d', L), 1:nShow, 'UniformOutput', false));
        xlabel(ax, 'lick contact (warped, equal width)', 'FontSize', 10);
    else
        xticklabels(ax, []);   % only the bottom panel carries the axis
    end
    ylabel(ax, 'Tongue length (norm)', 'FontSize', 10);

    yl = ylim(ax);
    tx = xr(1) + 0.60*diff(xr);
    text(ax, tx, yl(1) + 0.94*diff(yl), 'Actual', 'Color', [0 0 0], ...
        'FontWeight','bold', 'FontSize',10, 'Clipping','off');
    text(ax, tx, yl(1) + 0.80*diff(yl), sprintf('%s predicted', nameStr), 'Color', colP, ...
        'FontWeight','bold', 'FontSize',10, 'Clipping','off');
    box(ax, 'off')
end

function h = drawWarpedSet(ax, W, col, slotStart, cfg, nShow)
    % Draw one warped set (actual, or predicted) on ax: mean +/- 1.96 SEM, drawn
    % segment by segment so no line crosses the gap between two licks.
    segL = cfg.warpedSegmentLength;
    x  = (1:size(W,1))';
    nF = sum(isfinite(W), 2);
    mW = mean(W, 2, 'omitnan');
    sW = std(W, 0, 2, 'omitnan') ./ sqrt(max(nF,1));
    sW(nF < 2) = NaN;
    loW = mW - 1.96*sW;   hiW = mW + 1.96*sW;

    h = gobjects(1);
    if nargin < 6 || isempty(nShow), nShow = cfg.maxLickShow; end
    for L = 1:min(nShow, numel(slotStart))
        idxw = slotStart(L) : slotStart(L) + segL - 1;
        vf = isfinite(loW(idxw)) & isfinite(hiW(idxw));
        if any(vf)
            patch(ax, [x(idxw(vf)); flipud(x(idxw(vf)))], [loW(idxw(vf)); flipud(hiW(idxw(vf)))], ...
                col, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off');
        end
        hh = plot(ax, x(idxw), mW(idxw), '-', 'Color', col, 'LineWidth', 2.5);
        if L == 1, h = hh; else, set(hh, 'HandleVisibility', 'off'); end
    end
end

function ax = gridAxes(parent, nRow, nCol, k)
    % Axes at the position subplot(nRow, nCol, k) would use, inside any parent
    % (subplot only targets figures).
    r = ceil(k / nCol);  c = k - (r-1)*nCol;
    w = 0.84/nCol;  h = 0.80/nRow;
    x = 0.075 + (c-1)*(0.90/nCol);
    y = 0.90 - r*(0.88/nRow);
    ax = axes('Parent', parent, 'Position', [x y w h]);
end

function vline0(ax, col)
    % A vertical marker at t = 0 that does NOT rescale the axes: the limits are read
    % first and restored after, so the line is clipped to whatever the data set.
    yl = ylim(ax);
    line(ax, [0 0], yl, 'Color', col, 'LineStyle','--', 'HandleVisibility','off');
    ylim(ax, yl);
end

function [decIdx, peakL, nTrialL] = contactDecodingIndex(Yact, Ypred, alignSample, cfg)
    % Per-contact decoding index 1 - RMSE_dec / RMSE_zero, squared errors pooled
    % across trials. Contacts come from tongueRuns, and contact 1 is the protrusion
    % containing alignSample, as in the fitting window.
    nLA = cfg.nLicksAnalyze;
    nTr = size(Yact, 2);
    sseDec = zeros(nLA,1);  sseZero = zeros(nLA,1);
    nTrialL = zeros(nLA,1); peaks = nan(nLA, nTr);

    for tt = 1:nTr
        r = tongueRuns(Yact(:,tt), cfg);
        if isempty(r), continue; end
        first = pickFirstRun(r, alignSample);
        r = r(first : min(size(r,1), first + nLA - 1), :);
        for L = 1:size(r,1)
            seg = r(L,1):r(L,2);
            a = Yact(seg,tt);  p = Ypred(seg,tt);
            g = isfinite(a) & isfinite(p);
            if ~any(g), continue; end
            sseDec(L)  = sseDec(L)  + sum((a(g) - p(g)).^2);
            sseZero(L) = sseZero(L) + sum(a(g).^2);
            nTrialL(L) = nTrialL(L) + 1;
            peaks(L,tt) = max(a(g));
        end
    end

    decIdx = 1 - sqrt(sseDec ./ sseZero);   % n cancels inside the ratio
    decIdx(nTrialL == 0) = NaN;
    peakL = mean(peaks, 2, 'omitnan');
end

function h = ciBand(ax, x, M, col)
    % Mean +/- 95% CI across columns of M, drawn on ax. Returns the mean line handle.
    n  = sum(isfinite(M), 2);
    m  = mean(M, 2, 'omitnan');
    se = std(M, 0, 2, 'omitnan') ./ sqrt(max(n,1));
    hw = se .* tinv(0.975, max(n-1,1));
    hw(n < 2) = NaN;
    g = isfinite(m) & isfinite(hw);
    if any(g)
        fill(ax, [x(g); flipud(x(g))], [m(g)+hw(g); flipud(m(g)-hw(g))], col, ...
            'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility','off');
    end
    h = plot(ax, x, m, '-', 'Color', col, 'LineWidth', 2);
end

function [m, ci] = meanCI(M)
    % Mean and half-width of the 95% CI across the columns of M.
    n  = sum(isfinite(M), 2);
    m  = mean(M, 2, 'omitnan');
    ci = std(M, 0, 2, 'omitnan') ./ sqrt(max(n,1)) .* tinv(0.975, max(n-1,1));
    ci(n < 2) = NaN;
end

function col = pickColor(spec, g, ts, grpCols, tsCols)
    % Color by group when groups share an axes, otherwise by test set.
    if spec.overlayGroups
        col = grpCols(min(g, size(grpCols,1)), :);
    else
        col = tsCols(min(ts, size(tsCols,1)), :);
    end
end

function s = seriesLabel(spec, g, ts, n)
    if spec.overlayGroups && numel(spec.testSets) == 1
        s = sprintf('%s (n = %d)', spec.groupLabels{g}, n);
    elseif spec.overlayGroups
        s = sprintf('%s %s (n = %d)', spec.groupLabels{g}, spec.testSets(ts).name, n);
    else
        s = sprintf('%s (n = %d)', spec.testSets(ts).name, n);
    end
end

function finishAxes(ax, nLA, ylab, h, lbl, ylimVal)
    % Shared formatting for the per-contact axes: zero line, contact ticks, labels,
    % legend. ylimVal = [] leaves the y limits unchanged.
    if nargin < 6, ylimVal = []; end
    line(ax, [0.5 nLA+0.5], [0 0], 'Color','k', 'LineStyle','--', 'HandleVisibility','off');
    xlim(ax, [0.5 nLA+0.5]);  xticks(ax, 1:nLA);
    if ~isempty(ylimVal), ylim(ax, ylimVal); end
    xlabel(ax, 'contact number', 'FontSize', 11);
    ylabel(ax, ylab, 'FontSize', 11);
    if ~isempty(h)
        legend(ax, h, lbl, 'Location','best','FontSize',8);  legend(ax,'boxoff')
    end
    box(ax, 'off');
end

function drawStars(ax, contacts, ~, q, alpha)
    % Asterisks per contact that survives Benjamini-Hochberg at q < alpha,
    % graded by pStars: * q < alpha, ** q < 0.01, *** q < 0.001.
    % The third argument is unused.
    yl  = ax.YLim;
    pad = 0.04 * range(yl);
    for L = contacts(:)'
        if L > numel(q), continue; end
        st = pStars(q(L), [alpha 0.01 0.001]);
        if ~isempty(st)
            text(ax, L, yl(2) - 1.0*pad, st, 'Color','r', 'FontSize',16, ...
                'FontWeight','bold', 'HorizontalAlignment','center', 'Clipping','off');
        end
    end
    ylim(ax, yl);
end

function q = bhFDR(p)
    % Benjamini-Hochberg FDR over the finite entries of p; NaNs preserved in place.
    q = nan(size(p));
    fi = find(isfinite(p));
    if isempty(fi), return; end
    [ps, ord] = sort(p(fi));
    m  = numel(ps);
    qs = ps(:) .* (m ./ (1:m)');
    for ii = m-1:-1:1, qs(ii) = min(qs(ii), qs(ii+1)); end
    qs = min(qs, 1);
    qb = nan(m,1);  qb(ord) = qs;
    q(fi) = qb;
end

function S = tlLoadSession(loader, dateStr, spec, params, cfg)
    % TLLOADSESSION  Load one session, with an on-disk cache.
    % Returns a slim struct holding only the fields the decoder reads:
    %   S.time, S.trialdat, S.nUnitsTotal, S.tongueRaw, S.cluid, S.trialid,
    %   S.anm, S.date, S.Ntrials, S.lickL, S.goCue, S.trialTypes, S.hit
    % The cache key (cacheKey) is stored in the .mat and compared on read, so a
    % change to any loading parameter forces a reload. trialdat is stored as single;
    % ridgePath casts back to double.
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
    k.source  = 'slimExport';   % entries built from the slim export
    key = jsonencode(k);
end

function [kSel, sse] = selectLambdaCV5(X, y, rowTrial, cfg)
    % Lambda by k-fold cross-validation within the training set (k = cfg.cvFolds):
    % fit on k-1 folds, predict the held-out fold, and sum squared errors across
    % folds (longer trials carry proportionally more weight).
    % Folds are over whole trials, not rows: consecutive samples are strongly
    % autocorrelated, so row-wise folds would measure interpolation and pick too
    % small a lambda. rowTrial (from laggedDesign) gives each row's trial.
    % The caller takes the coefficients from the full training-set path at kSel.
    % cfg.lambdaFixed pins lambda to the nearest grid point instead.
    nK = numel(cfg.ridgeGrid);
    if ~isempty(cfg.lambdaFixed)
        sse = nan(nK,1);
        [~, kSel] = min(abs(log10(cfg.ridgeGrid(:)) - log10(cfg.lambdaFixed)));
        return
    end

    trials = unique(rowTrial(:));
    nT     = numel(trials);
    kF     = min(cfg.cvFolds, nT);
    assert(kF >= 2, ['Cross-validation needs at least 2 training trials; this ' ...
                     'session has %d. Lower cfg.cvFolds or check the split.'], nT);

    if isempty(cfg.cvSeed)
        ord = (1:nT)';
    else
        ord = randperm(RandStream('mt19937ar','Seed',cfg.cvSeed), nT)';
    end
    fold      = zeros(nT,1);
    fold(ord) = mod(0:nT-1, kF) + 1;

    sse   = zeros(nK,1);
    nUsed = 0;
    for f = 1:kF
        isTe = ismember(rowTrial(:), trials(fold == f));
        if ~any(isTe) || all(isTe), continue; end
        Rf  = ridgePath(X(~isTe,:), y(~isTe), cfg.ridgeGrid, cfg, false);
        Pf  = Rf.b0 + X(isTe,:)*Rf.beta;   % nTe x nK
        sse = sse + sum((y(isTe) - Pf).^2, 1)';
        nUsed = nUsed + 1;
    end
    assert(nUsed >= 2, 'Only %d usable folds; lambda cannot be cross-validated.', nUsed);
    [~, kSel] = min(sse);
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

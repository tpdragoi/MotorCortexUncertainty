%% FS4D_decodeJaw.m
%  Ridge decoder of jaw position from population spike rates, Delayed Reward Task.
%  Spike rates are z-scored and lagged, then fit per session with a ridge
%  penalty chosen by cross-validation over held-out training trials. The
%  decoding index is 1 - RMSE(decoded) / RMSE(zero baseline), computed for
%  each lick contact and averaged across sessions.
%  Produces Fig. S4D.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; Data = dataRoot in setPaths.m
%  ANALYSIS SETTINGS
%    params.alignEvent  'firstLick'
%    params.dt          1/300
%    params.smooth      10
%    params.quality     {'good','excellent'}
%    params.lowFR       0.01
%    params.window      -2.5 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


%% RUN SETTINGS (edit these two)
% RUN.nRepeats: number of full passes of the analysis (1 = a single run).
% Each pass re-draws the train/test split and the CV fold assignment; passes are
% not averaged or pooled. Each pass gets its own figure windows, titled with the pass number.
RUN.nRepeats = 1;

% CROSS-VALIDATION FOLDS for choosing lambda, inside the training set.
RUN.cvFolds  = 5;   % five-fold cross-validation
%% PATHS

% Sessions are read from the exported obj/kin files in the data folder set in
% setPaths.m, through slimMeta / slimToLegacy / loadSlimSession (in shared/).
% slimToLegacy rebuilds obj/params/kin for this script's params from the exported
% spike times; the pipeline functions it needs are in shared\pipelineCopies.
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
assert(exist(fullfile(repoRoot, 'setPaths.m'), 'file') == 2, ...
    'Cannot find setPaths.m. Run this script from its file, or cd to the repository root first.');
addpath(repoRoot);
cfgPaths = setPaths();   % data locations are set once, in setPaths.m
spec.dataDir = '';   % raw data folder not used
%% SETTINGS
rewardLickB = 4;   % rewardedLick value for conditions 5-7 and 9

%% TIMING: the three windows, all in seconds
% Everything is relative to t = 0, which IS the first lickport contact after the
% go cue (params.alignEvent = 'firstLick').

% FITTING WINDOW -- per trial, variable length.
%   go cue - preGoCue_s  ->  closing trough of contact 1's jaw cycle + postJawClose_s.
% Jaw displacement has no zero, so the window runs to the end of the first jaw
% cycle (its closing trough) rather than to the contact, which occurs mid-cycle.
% The trough comes from localJawSegments (the persistence-based cycle finder), so
% the fitting window, the decoding index and the warped figure share one definition of cycle 1.
cfg.preGoCue_s     = 0.20;   % window opens this long BEFORE the go cue
cfg.postJawClose_s = 0.005;   % window closes this long AFTER contact 1's cycle does
cfg.minWindowSamples = 10;   % trials with a shorter window are dropped

%% WHICH FITTING WINDOW (three options)
% 'contact1'        default: go cue - preGoCue_s  ->  contact 1's jaw cycle closes
% 'bout'            From winBoutStartPad_s AFTER CONTACT 2, to winBoutEndPad_s
%                   after the contact that ENDS the bout -- the first one not
%                   followed by another within winBoutGap_s. Defined ENTIRELY
%                   from contact times, so it matches the tongue script's 'bout'
%                   window sample for sample.
% 'lick2ToBoutEnd'  END OF LICK 2 -> END OF BOUT. Same end as 'bout', but the
%                   start is the CLOSING TROUGH of cycle 2 + winBoutStartPad_s,
%                   not contact 2's time. The jaw analogue of "after lick 2 has
%                   finished": contact 2 happens part way up the opening, so
%                   anchoring on it would put the start inside lick 2.
cfg.fitWindowMode     = 'contact1';   % 'contact1' | 'bout' | 'lick2ToBoutEnd'
cfg.winBoutStartPad_s = 0.02;   % opens this long after contact 2 / cycle 2 closing
cfg.winBoutEndPad_s   = 0.20;   % window closes this long after the bout's last contact
cfg.winBoutGap_s      = 0.75;   % a gap longer than this is what ENDS the bout

% NEURAL Z-SCORING WINDOW
cfg.zscoreWin_s = [-2.0 2.0];   % s, around the first port contact
cfg.minBaselineStd = 1e-6;   % units below this are DROPPED, not divided by

% NEURAL LAG CONTEXT (lag set is -preBins : +postBins, inclusive)
cfg.lagPre_s  = 0.06;
cfg.lagPost_s = 0.04;

% ANALYSIS WINDOW for the per-contact decoding index
cfg.analysisWin_s = [-0.09 3.0];

%% TARGET: jaw displacement
% A NaN in jaw_ydisp_view1 means the tracker lost the jaw marker, not a zero
% displacement, so NaNs are never zero-filled:
%   - training rows with a NaN target are dropped (same okTr mask as NaN predictors);
%   - every score (ridge R^2, decoded RMSE, baseline RMSE) uses one common finite
%     mask over actual and predicted, so numerator and denominator cover the same samples;
%   - predictions are not masked, so the full-trial predicted trace is continuous;
%   - segmentation runs on a gap-filled copy (findpeaks cannot take NaN) and its
%     indices are applied to the NaN-bearing arrays, so interpolated samples are never scored.
cfg.jawFeature  = 'jaw_ydisp_view1';
cfg.normPctLow  = 2;   % lower scaling percentile
cfg.normPctHigh = 99;   % upper scaling percentile

% false: percentile-scaled values may run slightly outside [0,1] (peaks are not clipped).
% true: clip the scaled target to [0,1].
cfg.clipNormalizedRange = false;
cfg.reportDropout       = true;   % print the % of untracked jaw samples per session

%% RIDGE
cfg.ridgeGrid       = logspace(-2, 8, 41);
cfg.standardizeCols = true;   % penalty lambda*sum(beta_j^2 * sigma_j^2);
% required by cfg.ridgeImpl = 'matlab'

cfg.ridgeImpl = 'eig';   % closed-form ridge via eigendecomposition: 'eig' | 'matlab'
% 'eig' is algebraically identical to MATLAB's ridge and much faster when p is
% large ('matlab' does one least-squares solve per lambda).
cfg.ridgeCheckEquivalence = false;   % fit the first session both ways and print
% the max relative difference (expect ~1e-11)

% LAMBDA -- k-fold cross-validation inside the training set, k = cfg.cvFolds.
% Folds are over whole trials, never rows; the test trials are never used in fitting.
cfg.lambdaRule  = 'cv5';   % 'cv5' = k-fold CV (k from cfg.cvFolds, not from this label)
cfg.cvFolds     = RUN.cvFolds;   % set in the RUN block at the top of this file
cfg.cvSeed      = [];   % [] = cross-validation folds in trial order, no shuffle.
%                        A number seeds the shuffle and each pass offsets it by repIdx.
RUN.cvSeed0     = cfg.cvSeed;   % base seed; each pass offsets it by repIdx
cfg.lambdaFixed = [];   % [] = choose lambda by cfg.lambdaRule; a number pins lambda instead

%% DECODING INDEX: what it is measured against
% The baseline is zero: the null model predicts a scaled jaw of 0 at every sample,
% the same null as the tongue decoders.
% The target is scaled as (jaw - pLow)/(pHigh - pLow) over all tracked samples in
% the session, so 0 sits slightly below the resting jaw; absolute index values are
% best compared within a session, across contacts (contact L vs contact 1).
% A constant baseline shrinks when jaw excursion shrinks, so peakL (per-contact
% excursion) is carried alongside the index.

cfg.nLicksAnalyze = 8;

%% JAW CYCLE SEGMENTATION
% Jaw displacement has no resting floor to threshold, so cycles are found by
% persistence: build the alternating max/min sequence of the smoothed trace, then
% repeatedly delete the least prominent adjacent extremum pair until every remaining
% pair clears cyclePersistFrac x the trial's robust (p5-p95) jaw range. Surviving
% minima are the cycle boundaries, with exactly one peak between consecutive minima.
% Cycles are matched to lick contacts containment first, then by nearest peak.
cfg.segAnchor          = 'persistCycle';   % 'persistCycle' | 'contactCycle' | 'peakDetect'
cfg.segMinTime_s       = -0.1;   % ignore peaks before this time
cfg.smoothSamples      = 3;   % movmean applied before detection
cfg.cyclePersistFrac   = 0.3;   % an excursion smaller than this fraction of the
% trial's p5-p95 jaw range is not a lick cycle
cfg.cyclePersistFloor  = 0.03;   % absolute floor, for trials with a tiny range
cfg.cycleMaxHalfILIFrac = 0.75;   % a cycle half may not exceed this x the trial's
% median inter-lick interval (stops the last
% cycle growing a multi-second closed tail)
cfg.rescueMissedCycles   = true;
cfg.rescueMaxSpanILIFrac = 2.2;
cfg.rescueMinAmpFrac     = 0.12;
cfg.rescueEdgeContacts   = true;
cfg.contactMatchTol_s    = 0.10;
cfg.contactMatchILIFrac  = 0.7;   % effective tolerance = max(tol_s, this x median ILI)
cfg.peakSearchFrac       = 0.5;
cfg.troughRiseFrac       = 0.25;
cfg.troughRiseFloor      = 0.01;
cfg.troughMinDepthFrac   = 0.6;
cfg.troughFlatTol        = 0;
cfg.allowEdgeTruncated   = true;
cfg.segMinLenSamples     = 5;
cfg.segMaxLenSamples     = Inf;
cfg.minCycleAmp          = 0;
cfg.segMode              = 'troughToTrough';   % 'troughToTrough' | 'fixedHalfWin'
cfg.segMaxHalfCycleFrac  = 1.5;
% Used only by the 'peakDetect' anchor:
cfg.jawLickIndexRule   = 'contactMatched';   % 'contactMatched' | 'peakOrder'
cfg.peakMinHeight      = 0.08;
cfg.peakProminence     = 0.05;
cfg.peakMinDist        = 5;
cfg.troughProminence   = 0.0001;
cfg.troughMinDist      = 7;
cfg.troughMaxValue     = 0.25;
cfg.segHalfWin         = 5;
cfg.requireTroughBelowMax = false;
% localJawSegments reads this name; it is set from cfg.nLicksAnalyze.
cfg.numLicksToAnalyze = cfg.nLicksAnalyze;
cfg.reportSegmentation = true;   % print how many contacts got a cycle, per session

%% BOUT REQUIREMENT (one toggle, on or off)
% Keep a trial only if the animal produced a real BOUT on it: at least
% boutMinLicks port contacts IN A ROW, consecutive ones no more than
% boutMaxILI_s apart, all of them inside boutWin_s of the go cue.
% cfg.requireBout applies to the whole trial pool (training and test trials).
cfg.requireBout   = true;   % true | false
cfg.boutMinLicks  = 6;
cfg.boutWin_s     = 1.50;
cfg.boutMaxILI_s  = 0.25;

%% STATISTICS
cfg.alpha      = 0.05;

%% SESSION CACHE
% Loading dominates run time and does not depend on cfg, so loaded sessions are
% cached and later runs skip loading.
% The cache key is the set of params fields that affect the loaded data; changing
% any of them recomputes the entry.
% Budget roughly (nTime x nUnits x nTrials x 4 bytes) per session: trialdat is
% stored as single, so a 1500 x 300 x 400 session is about 700 MB. Point
% cacheDir at a scratch disk, and delete the folder any time to force a reload.
cfg.useCache = true;
cfg.cacheDir = fullfile(cfgPaths.cacheRoot, 'tlDecodeCache');

%% HOUSEKEEPING
cfg.rngSeed          = [];
cfg.sessionsToRun    = [];   % [] = all

% DIAGNOSTIC FIGURES. [] means "choose automatically", by spec.diagRule:
%   'random'  -- cfg.nDiagSessions sessions drawn at random from the ones that
%                fit, using cfg.rngSeed (reproducible).
%   'groups'  -- every session in spec.diagGroups.
% Put explicit session numbers in cfg.diagSessions to override either.
% Each diagnostic session is a tab: one figure with the six-panel fit diagnostics
% and one with the warped lick train.
cfg.diagSessions     = [];
cfg.nDiagSessions    = 0;

% cfg.warpAllSessions (the warped figure shows whether the cycle segmentation found the licks):
%   true   every fitted session gets a warped tab, whatever the diagnostic rule picked
%   false  the warped figure follows the same session list as the six-panel one
cfg.warpAllSessions  = false;
cfg.warnPredictorCount = 8000;

% 'slots' warped layout: first slot at 10, slots every 50 (a 20-sample NaN gap
% between licks); up to 10 licks extracted and the first 8 drawn.
cfg.warpedSegmentLength = 30;   % samples each excursion is resampled onto
cfg.warpedSlotFirst     = 10;   % where lick 1's slot starts
cfg.warpedSlotSpacing   = 50;   % slot pitch; spacing - segment = the NaN gap
cfg.maxLicksToExtract   = 10;
cfg.maxLickShow         = 8;

%   warpedInterp     'pchip' | 'linear' ('linear' cannot overshoot at segment ends).
%   warpedTickCenter false puts the "contact k" label at the left edge of the
%                    slot; true centers it under the excursion.
cfg.warpedInterp     = 'pchip';   % 'pchip' | 'linear'

% HOW THE WARPED FIGURE IS BUILT.
%   'continuous'  a piecewise-affine map per inter-lick interval takes this trial's
%                 contact times onto the across-trial template (spacing warpDelta_s),
%                 and the trace is resampled through the inverse map: one continuous
%                 trace on a warped time axis, with a dotted line at each template contact.
%   'slots'       the tongue scripts' layout: each cycle resampled into its own
%                 fixed-width slot with NaN between them.
% 'continuous' keeps the jaw trace between licks, which 'slots' discards.
cfg.warpStyle    = 'continuous';   % 'continuous' | 'slots'
cfg.warpDelta_s  = 0.15;   % fixed inter-lick spacing of the template

% HOW FAR THE WARP REACHES vs HOW MUCH IS SHOWN -- two different numbers.
%   cfg.warpNumLicks   how many contacts the TEMPLATE spans, i.e. how much of the
%                      bout actually gets warped. Set it generously: any part of
%                      the trace past the last mapped interval is handled by
%                      cfg.warpEdgeMode rather than by a real landmark, so a
%                      short template leaves a long extrapolated tail. 12 covers
%                      essentially every bout in this data.
%   spec.warpShowLicks how many contacts are LABELED and kept inside the x
%                      limits (set per study).
cfg.warpNumLicks = 12;

% What happens OUTSIDE the mapped intervals -- before the first and after the
% last landmark: 'extend' carries the nearest affine map outwards, 'nan' leaves
% them empty, 'identity' samples them on the raw clock.
cfg.warpEdgeMode = 'extend';   % 'extend' | 'nan' | 'identity'

% With two test sets, draw both on the same warped axes -- solid actual, dashed
% predicted, each in that test set's color. Ignored with one test set.
cfg.warpSplitTestSets = true;
cfg.warpedTickCenter = false;


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

% Index 8 is the R1 pool and index 9 the second reward-lick pool (rewardLickB).
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
%% STUDY: R1 vs R4
spec.name = 'R1R4';
spec.sessionDates = { ...
    '2023-02-21','2023-02-22','2023-02-23','2023-02-24', ...   % TD1d
    '2023-02-21','2023-02-24','2023-02-25','2023-03-19', ...   % TD4d
    '2024-11-11','2024-11-12','2024-11-13','2024-11-21','2024-11-22','2024-11-24','2024-11-25', ...   % TD13d
    '2024-11-24','2024-11-25','2024-11-26','2024-11-27', ...   % TD15d
    '2024-09-06','2024-09-07','2024-09-08','2024-09-09','2024-09-10','2024-09-22', ...   % TD8d
    '2025-06-17','2025-06-18','2025-06-19','2025-06-20','2025-06-21', ...   % TD22d
    '2025-06-17','2025-06-18','2025-06-19','2025-06-20','2025-06-21' };   % TD23d
spec.sessionLoaders = [ ...
    repmat({@loadTD1_neural},  1,4), repmat({@loadTD4_neural}, 1,4), ...
    repmat({@loadTD13_neural}, 1,7), repmat({@loadTD15_neural},1,4), ...
    repmat({@loadTD8_neural},  1,6), repmat({@loadTD22_neural},1,5), ...
    repmat({@loadTD23_neural}, 1,5) ];
spec.groupMaps = { ...
    [1 1 0 1  2 0 2 1  0 1 2 2 0 0 0  2 1 1 2  0 0 0 0 0 0  0 1 2 2 1  2 2 2 1 0] , ...   % M1
    [2 2 2 2  0 2 0 2  0 0 0 0 0 0 0  0 0 0 0  0 0 0 0 0 0  2 2 1 1 0  1 1 1 2 2] , ...   % ALM
    [0 0 0 0  0 0 0 0  2 0 1 0 1 1 1  0 2 2 1  1 2 2 2 1 1  0 0 0 0 0  0 0 0 0 0] };   % S1
spec.groupLabels = {'M1','ALM','S1'};
% TRAIN pool = conditions 8 and 9 together (R1 + R4), type-stratified quota.
spec.poolCondIdx = [8 9];
spec.splitRule   = 'typeStratified';
spec.trainFrac   = 0.60;   % fraction of the pool used for training
spec.testSets    = struct('name', {'R1','R4'}, 'condIdx', {8, 9});
spec.trialCaps = { 'TD13d','2024-11-11', 278; ...
                   'TD8d', '2024-09-07', 313 ; ...
                   'TD8d', '2024-09-09', 298 ; ...
                   'TD22d','2025-06-20', 170 };
spec.singleProbeSessions = cell(0,2);
spec.lickCountWin_s      = 1.25;

spec.minTestTrials    = 10;   % 'lickCountFrac': never leave fewer long-bout
% trials than this to score on
spec.testFromManyOnly = true;   % 'lickCountFrac': test on long-bout trials only,
% as the plain 'lickCount' rule does

spec.warpShowLicks    = 8;   % contacts LABELED on the warped figure

%% bout filter + diagnostics
spec.applyBoutFilter = true;   % cfg.requireBout is the on/off switch
spec.diagRule        = 'random';   % 'random' = cfg.nDiagSessions sessions at random
spec.diagGroups      = [];   % used only by diagRule 'groups'

%% colors
% testSets are {'R1','R4'} in that order: R1 DARK BLUE, R4 LIGHT BLUE.
spec.tsColors  = [0.00 0.00 0.55;   % R1  dark blue
                   0.35 0.70 0.95];   % R4  light blue
spec.grpColors = [];
spec.predColor = [0.00 0.00 0.55];   % prediction on the diagnostics

%% figure behavior
spec.overlayGroups = false;   % false = never draw two groups on one axes
spec.showFigCI     = true ;   % Figure 1: decoding index with error bars
spec.showFigDelta  = false;   % Figure 3: delta from contact 1
spec.deltaGroups   = [];   % groups on the delta figures ([] = all)
spec.figRef       = 'Fig. S4D';   % the manuscript panel this file produces
spec.vsC1Tail     = 'both';   % 'both' | 'left' | 'right'
spec.vsC1Test     = false;   % per-contact test vs contact 1
spec.betweenTest  = true ;   % per-contact test BETWEEN the two trial types
spec.groupTest    = 'none';   % 'none' | 'ranksum' | 'crossTask'
spec.groupSummaryLicks = 3:8;   % contacts averaged into the single summary value
spec.summaryFile  = fullfile(cfgPaths.cacheRoot, 'decodeSummary_jaw_r14.mat');
spec.crossTaskFile = fullfile(cfgPaths.cacheRoot, 'decodeSummary_jaw_r1.mat');
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

% mask maps to the run set, then work out which probe each session needs
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
        logf('  fitting window [contact1]: go cue -%.0f ms  ->  contact-1 jaw cycle closes +%.0f ms (per trial)\n', ...
            1000*cfg.preGoCue_s, 1000*cfg.postJawClose_s);
    case 'bout'
        logf('  fitting window [bout]: contact 2 +%.0f ms  ->  bout-end contact +%.0f ms (bout ends on a gap > %.2f s)\n', ...
            1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s, cfg.winBoutGap_s);
        fprintf('    CONTACT 1 IS OUTSIDE THE TRAINING WINDOW -- it is the extrapolation, not contacts 2-8.\n');
    case 'lick2toboutend'
        logf('  fitting window [lick2ToBoutEnd]: cycle 2 CLOSES +%.0f ms  ->  bout-end contact +%.0f ms (gap > %.2f s)\n', ...
            1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s, cfg.winBoutGap_s);
        fprintf('    CONTACT 1 AND CYCLE 2 ARE OUTSIDE THE TRAINING WINDOW -- contact 1 is the extrapolation.\n');
    otherwise
        error('cfg.fitWindowMode must be ''contact1'', ''bout'' or ''lick2ToBoutEnd''.');
end
logf('  z-score window: %.1f to %.1f s | lag context %.0f to +%.0f ms\n', ...
    cfg.zscoreWin_s(1), cfg.zscoreWin_s(2), -1000*cfg.lagPre_s, 1000*cfg.lagPost_s);
logf('  decoding index baseline: ZERO (1 - RMSE_dec/RMSE_0, same null as the tongue scripts)\n');
logf('  jaw cycles: %s segmentation\n', cfg.segAnchor);
if applyBout
    logf('  BOUT FILTER ON: >= %d contacts in a row, ILI <= %.0f ms, all within %.2f s of the go cue\n', ...
        cfg.boutMinLicks, 1000*cfg.boutMaxILI_s, cfg.boutWin_s);
elseif cfg.requireBout
    logf('  bout filter: not applied in this script (spec.applyBoutFilter = false)\n');
else
    logf('  bout filter: OFF (cfg.requireBout = false)\n');
end
switch lower(cfg.lambdaRule)
    case 'cv5'
        logf(['  lambda: %d-fold cross-validation over TRAINING TRIALS, ' ...
                 'refit on all training rows at the winning value\n'], cfg.cvFolds);
end
logf('  diagnostics: sessions %s\n', mat2str(diagList));
if cfg.warpAllSessions
    warpList = runList(~cellfun(@isempty, probesOf(runList)));
else
    warpList = diagList;
end
logf('  warped figure: sessions %s\n', mat2str(warpList));

% Colors come from the study block so each task keeps its own scheme; the shared
% code never hard-codes one.
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
assert(size(S.jawRaw,1) == nT, 'sess %d: kinematics %d time samples, time %d.', ...
    sessionIdx, size(S.jawRaw,1), nT);

% ---------------- TARGET: JAW DISPLACEMENT ----------------
% Percentile-scaled, not clipped by default, NaN left as NaN (tracker dropout; see
% TARGET settings). Every score uses a finite mask shared by actual and predicted.
jawRaw  = S.jawRaw;
visMask = ~isnan(jawRaw);
visVals = jawRaw(visMask);
pLo = prctile(visVals, cfg.normPctLow);
pHi = prctile(visVals, cfg.normPctHigh);
assert(pHi > pLo, 'sess %d: jaw percentiles are degenerate (p%g = p%g).', ...
    sessionIdx, cfg.normPctLow, cfg.normPctHigh);
jawNorm = (jawRaw - pLo) ./ (pHi - pLo);
if cfg.clipNormalizedRange
    jawNorm = min(max(jawNorm, 0), 1);   % flattens the largest openings
end
if cfg.reportDropout
    logf('  [jaw]   %.1f%% of samples untracked (NaN) | scaling p%g = %.4g, p%g = %.4g\n', ...
        100*mean(~visMask(:)), cfg.normPctLow, pLo, cfg.normPctHigh, pHi);
end

% ---------------- JAW CYCLES, ONCE PER SESSION ----------------
% localJawSegments is run ONCE over every usable trial, on the analysis window,
% and its output defines the per-trial fitting window, the per-contact decoding
% index and the warped figure, so "contact 1" means one thing throughout.
% It runs on the ANALYSIS WINDOW, not the full trace, because the persistence
% threshold is a percentile of the trace it is handed: over the full -2.5 to 5 s
% trace most samples are quiet, the p5-p95 range collapses, and the threshold
% drops with it.
maxUsable = min([size(S.jawRaw,2), size(S.trialdat,3), S.Ntrials]);
capRow = find(strcmp(anm, spec.trialCaps(:,1)) & strcmp(dte, spec.trialCaps(:,2)), 1);
if ~isempty(capRow), maxUsable = min(maxUsable, spec.trialCaps{capRow,3}); end

aIdx = find(traceTime >= cfg.analysisWin_s(1) & traceTime <= cfg.analysisWin_s(2));
aTime = traceTime(aIdx);
allTrials = (1:maxUsable)';
tSeg = tic;
[segCell, nSegDisagree, segDiag] = localJawSegments(jawNorm(aIdx, allTrials), ...
    allTrials, S.goCue, S.lickL, aTime, cfg);
if cfg.reportSegmentation
    haveC = segDiag(:,1) > 0;
    logf('  [seg]   %d trials in %.1f s | contacts %.1f/trial, cycles found %.1f/trial (%.0f%%) | %d rescued | %d edge-truncated\n', ...
        maxUsable, toc(tSeg), mean(segDiag(haveC,1)), mean(segDiag(haveC,2)), ...
        100*sum(segDiag(:,2))/max(sum(min(segDiag(:,1), cfg.nLicksAnalyze)),1), ...
        sum(segDiag(:,6)), sum(segDiag(:,5)));
    if nSegDisagree > 0
        logf('          %d trial(s) where cycle order and contact order disagree (see localJawSegments)\n', nSegDisagree);
    end
end

% ---------------- PER-TRIAL FITTING WINDOWS ----------------
% goCueRel is negative: the cue precedes contact 1.
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
    winCell{tr}   = jawFitWindow(segCell{tr}, aIdx, traceTime, goCueRel(tr), post - post(1), cfg);
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

% ---------------- MEAN JAW INSIDE THE FITTING WINDOW ----------------
% Every trial's own window, laid back on the real time axis and averaged across
% trials at each sample. Windows are variable length (the go-cue latency and the
% contact-1 duration both vary), so a sample is averaged over however many trials
% actually cover it -- jawWinN records that count.
Mjw = nan(nT, numel(poolTrials));
for tt = 1:numel(poolTrials)
    ki = winCell{poolTrials(tt)};
    Mjw(ki,tt) = jawNorm(ki, poolTrials(tt));
end
jawWinMean = mean(Mjw, 2, 'omitnan');
jawWinN    = sum(isfinite(Mjw), 2);
clear Mjw

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
% Hard partition on contact count: every short-bout trial trains, every
% long-bout trial tests; nothing is sampled, so the split is deterministic.
        testTrials  = poolTrials(lickCount(poolTrials) >= spec.lickCountMin);
        trainTrials = poolTrials(lickCount(poolTrials) <  spec.lickCountMin);

    case 'lickCountFrac'
% Training quota = trainFrac x pool, filled few-lick trials first: all short-bout
% trials go in, then the quota is topped up with random long-bout trials.
% E.g. 100 pool trials, trainFrac 0.7, 40 short-bout trials -> all 40 short + 30
% long train, 30 long test.
% Two guards:
%   spec.minTestTrials    the top-up is capped so at least this many
%                         long-bout trials are LEFT to score. Without it
%                         a high trainFrac on a session with few short
%                         trials eats the test set entirely.
%   spec.testFromManyOnly true keeps the test set to long-bout trials, so the
%                         per-contact n at contacts 4-8 is not diluted by
%                         trials that never had a contact 4; leftover short-bout
%                         trials are unused. false puts every non-training
%                         trial in the test set.
% Unlike 'lickCount' (a pure domain shift), a non-zero top-up makes this partly a
% within-distribution split, which raises the decoding index. The split is random;
% set cfg.rngSeed for reproducible runs.
        isMany   = lickCount(poolTrials) >= spec.lickCountMin;
        fewPool  = poolTrials(~isMany);
        manyPool = poolTrials(isMany);
        nWant    = round(spec.trainFrac * numel(poolTrials));
        nFewTake = min(nWant, numel(fewPool));
        nTopUp   = min(nWant - nFewTake, max(numel(manyPool) - spec.minTestTrials, 0));
        trainFew  = [];  trainMany = [];
        if nFewTake > 0, trainFew  = fewPool(randperm(numel(fewPool),  nFewTake)); end
        if nTopUp   > 0, trainMany = manyPool(randperm(numel(manyPool), nTopUp));  end
        trainTrials = [trainFew(:); trainMany(:)];
        trainTrials = trainTrials(randperm(numel(trainTrials)));
        if spec.testFromManyOnly
            testTrials = manyPool(~ismember(manyPool, trainMany));
        else
            testTrials = poolTrials(~ismember(poolTrials, trainTrials));
        end
        logf(['  [split] quota %d of %d (%.0f%%) | few-lick %d of %d used | ' ...
                 'many-lick top-up %d of %d | %d left to test\n'], ...
            nWant, numel(poolTrials), 100*spec.trainFrac, nFewTake, numel(fewPool), ...
            nTopUp, numel(manyPool), numel(testTrials));
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
yTrain = stackWindows(jawNorm, winTrain, trainTrials);
yTest  = stackWindows(jawNorm, winTest,  testTrials);
okTr = all(isfinite(Xtrain),2) & isfinite(yTrain);
okTe = all(isfinite(Xtest), 2) & isfinite(yTest);

winLen = cellfun(@numel, winTrain);
logf('  [design] %d units x %d taps = %d predictors | window %d-%d samples (%.0f-%.0f ms) | %d train rows\n', ...
    nUnits, numel(lags), nPred, min(winLen), max(winLen), ...
    1000*params.dt*min(winLen), 1000*params.dt*max(winLen), sum(okTr));
logf('  [target] %.1f%% of window samples untracked (dropped from the fit)\n', 100*mean(~isfinite(yTrain)));

% ---------------- FIT: ONE RIDGE PATH ON ALL THE TRAINING ROWS -------------
% One ridge path over all the training rows gives the coefficients at every
% candidate lambda from a single eigendecomposition. Lambda itself is then
% chosen by k-fold cross-validation (k = cfg.cvFolds) over TRAINING TRIALS (cfg.lambdaRule =
% 'cv5'), and the coefficients are read off this full-training-set path at the
% winning lambda -- so "refit on all training trials at the lambda minimizing
% error" costs nothing extra. The test set is not touched by any of this.
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
% across folds, then refit on all training rows at the winning lambda --
% which is what R already holds, so only the index changes.
        [kSel, cvSSE] = selectLambdaCV5(Xtrain(okTr,:), yTrain(okTr), rowTrain(okTr), cfg);
        lamHow = sprintf('%d-fold trial-wise CV', cfg.cvFolds);
% the curve the diagnostic figure draws under panel 1
        lamCurve = cvSSE;
        lamCurveName = sprintf('%d-fold CV SSE', cfg.cvFolds);
        lamShort = sprintf('%d-fold CV', cfg.cvFolds);
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
% A lagged linear model applied over a whole trial is a filter, so the full-trial
% prediction needs one small mat-vec per lag per trial, not a materialized
% nT x (nUnits*nLags) design matrix per trial.
predFull = predictFullTrials(spikes(:,:,testTrials), lags, b0, beta, nUnits);

% ---------------- DECODING INDEX vs the ZERO baseline, per test set --------
% Per contact L:  1 - RMSE(actual, predicted) / RMSE(actual, 0)
% over exactly the samples localJawSegments assigned to cycle L, and over ONE
% finite mask shared by actual and predicted so the two RMSEs always cover the
% same samples. The baseline is the constant 0 for every trial and every
% session -- see the DECODING INDEX block in SETTINGS for what that 0 means.
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
    [decIdx(:,ts), peakL(:,ts), nTrialL(:,ts)] = jawDecodingIndex( ...
        jawNorm(aIdx, testTrials(sel)), predFull(aIdx, sel), ...
        segCell(testTrials(sel)), cfg);
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
res(k).tlTime    = traceTime;   % for the jaw-in-window panel
res(k).tlWinMean = jawWinMean;
res(k).tlWinN    = jawWinN;

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
gn = lamCurve;  gn(~isfinite(gn)) = NaN;
gn = (gn - min(gn)) ./ max(range(gn(isfinite(gn))), eps);   % rescaled for shape only
plot(ax, logK, gn, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
ylim(ax, [-0.2 1]);
line(ax, [logK(kSel) logK(kSel)], [-0.2 1], 'Color','k', 'HandleVisibility','off');
line(ax, [logK(1) logK(end)], [0 0], 'Color','k', 'LineStyle',':', 'HandleVisibility','off');
xlabel(ax, 'log_{10} \lambda'); ylabel(ax, sprintf('R^2   (%s rescaled)', lamCurveName));
legend(ax, {'train','test',[lamCurveName ' (scaled)']}, 'Location','southwest','FontSize',7);
if isempty(cfg.lambdaFixed), lamSrc = lamShort; else, lamSrc = 'FIXED'; end
title(ax, sprintf('\\lambda by %s = %.3g | edf %.0f of %d', lamSrc, cfg.ridgeGrid(kSel), ...
    R.edf(kSel), nPred+1), 'FontSize',9,'FontWeight','normal');
box(ax,'off')

% ---- 2: how long the per-trial windows came out ----
ax = gridAxes(tb, 2, 3, 2);
histogram(ax, winLen * params.dt * 1000, 20, 'FaceColor',[0.35 0.35 0.35], 'EdgeColor','none');
xlabel(ax, 'fitting window length (ms)'); ylabel(ax, 'train trials');
switch lower(cfg.fitWindowMode)
    case 'bout',           winTxt = sprintf('contact 2 +%.0f ms -> bout end +%.0f ms', 1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s);
    case 'lick2toboutend', winTxt = sprintf('cycle 2 closes +%.0f ms -> bout end +%.0f ms', 1000*cfg.winBoutStartPad_s, 1000*cfg.winBoutEndPad_s);
    otherwise,             winTxt = sprintf('go cue -%.0f ms -> cycle 1 closes +%.0f ms', 1000*cfg.preGoCue_s, 1000*cfg.postJawClose_s);
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
    Mact(ki,tt)  = jawNorm(ki, trainTrials(tt));
    Mpred(ki,tt) = Ptr(rowTrain == tt, kSel);
end
covFrac = 0.2;
covN    = sum(isfinite(Mact), 2);
hA = ciBand(ax, traceTime, Mact,  [0 0 0]);
hP = ciBand(ax, traceTime, Mpred, predCol);
gd = find(covN >= covFrac*nTrain);
if ~isempty(gd), xlim(ax, [traceTime(gd(1))-0.02, traceTime(gd(end))+0.02]); end
mAct = mean(Mact, 2, 'omitnan');
yl = [min(mAct), max(mAct)];
if all(isfinite(yl)) && diff(yl) > 0
    ylim(ax, yl + [-0.35 0.35]*max(diff(yl), 0.1));   % jaw has no [0 1] bound
end
vline0(ax, 'r');
xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'jaw displacement (norm)');
legend(ax, [hA hP], {'actual','predicted'}, 'Location','northwest','FontSize',7);
title(ax, {'TRAIN, in the fitting window only', ...
    sprintf('>= %.0f%% trial coverage | window mode: %s', 100*covFrac, cfg.fitWindowMode)}, ...
    'FontSize',9,'FontWeight','normal');
box(ax,'off')

% ---- 4: the whole trial, TEST ----
ax = gridAxes(tb, 2, 3, 4);  hold(ax,'on')
h1 = plot(ax, traceTime, mean(jawNorm(:,testTrials), 2, 'omitnan'), 'k', 'LineWidth', 2);
h2 = plot(ax, traceTime, mean(predFull, 2, 'omitnan'), 'Color', predCol, 'LineWidth', 2);
xlim(ax, [-0.5 1.5]);  vline0(ax, 'r');
xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'jaw displacement (norm)');
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
    h1 = plot(ax, traceTime, jawNorm(:, testTrials(tt)), '-', 'Color',[0 0 0], 'LineWidth',1.25);
    h2 = plot(ax, traceTime, predFull(:,tt), '-', 'Color', predCol, 'LineWidth',1.75);
    xlim(ax, [-0.5 1.5]);  vline0(ax, 'r');
    xlabel(ax, 'time from contact 1 (s)'); ylabel(ax, 'jaw displacement (norm)');
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
% 'slots' style: each cycle resampled to cfg.warpedSegmentLength points in its own slot, mean
% +/- 1.96 SEM, drawn segment by segment so nothing is joined across the gaps.
if wantWarp
if isempty(warpFig) || ~ishandle(warpFig)
    warpFig = figure('Name', sprintf('%s%s - warped lick train', passTag, spec.name), ...
        'Position', [110 110 1050 620]);
    warpTG  = uitabgroup(warpFig);
end
wt = uitab(warpTG, 'Title', tabTitle);

if strcmpi(cfg.warpStyle, 'continuous')
% One affine map per inter-lick interval carries this trial's contacts onto
% the across-trial median template, and the trace is resampled through the
% inverse of that map, giving one continuous trace on a warped time axis.
    nHere = numel(testTrials);
    gcRelList   = nan(nHere,1);
    contactCell = cell(nHere,1);
    for tr = 1:nHere
        [gcRelList(tr), contactCell{tr}] = localJawWarpLicks(testTrials(tr), S.goCue, S.lickL, cfg.warpNumLicks);
    end
    medTemplate = localJawWarpMedian(contactCell, cfg.warpNumLicks, cfg.warpDelta_s);
    if any(isnan(medTemplate))
        ax = axes('Parent', wt, 'Position', [0.09 0.15 0.86 0.72]);
        title(ax, {hdr, 'insufficient lick data to build the warp template'}, ...
            'FontSize',10,'FontWeight','normal');
        axis(ax,'off')
    else
% localJawWarpFits takes the median go cue from whatever list it is
% handed, so fitting the two conditions separately would anchor them on
% slightly different pre-lick offsets and the panels would be comparing
% two shifted axes. Built once, split afterwards, same x for both.
        warpFits = localJawWarpFits(contactCell, gcRelList, medTemplate);
        gridT = traceTime;
        Wa = nan(numel(gridT), nHere);  Wp = nan(numel(gridT), nHere);
        for tr = 1:nHere
            Wa(:,tr) = localJawWarpResample(traceTime, jawNorm(:,testTrials(tr)), gcRelList(tr), ...
                contactCell{tr}, medTemplate, warpFits(tr,:), gridT, cfg.warpEdgeMode);
            Wp(:,tr) = localJawWarpResample(traceTime, predFull(:,tr), gcRelList(tr), ...
                contactCell{tr}, medTemplate, warpFits(tr,:), gridT, cfg.warpEdgeMode);
        end
        nShow = min(spec.warpShowLicks, cfg.warpNumLicks);

% ---- ONE PANEL PER CONDITION, STACKED ----
% With two test sets (R1 vs R4, R1 vs R6) the tab holds two axes, R1 on
% top, in the order spec.testSets declares them. Each panel is one
% condition: actual black, predicted in that condition's color, labeled
% in-plot. Only the bottom panel carries the x axis.
        if nTS >= 2 && cfg.warpSplitTestSets, panels = 1:nTS; else, panels = 1; end
        nPan = numel(panels);
        nStr = '';
        for pp = 1:nPan
            ts = panels(pp);
            if nPan == 1 || isempty(spec.testSets(ts).condIdx)
                selTS = true(nHere,1);
            else
                memberTS = unique(cell2mat(S.trialid(spec.testSets(ts).condIdx)'));
                selTS = ismember(testTrials, memberTS);
            end
            axP = gridAxes(wt, nPan, 1, pp);
            if ~any(selTS)
                axis(axP, 'off');  continue
            end
            if nPan == 1
                colP = predCol;  nameStr = 'Model';
            else
                colP = tsCols(min(ts, size(tsCols,1)), :);  nameStr = spec.testSets(ts).name;
            end
            drawWarpedPanel(axP, gridT, Wa(:,selTS), Wp(:,selTS), medTemplate, ...
                nShow, colP, nameStr, pp == nPan, cfg);
            nStr = [nStr sprintf('%s n = %d   ', nameStr, sum(selTS))];   %#ok<AGROW>
        end
        annotation(wt, 'textbox', [0 0.955 1 0.04], 'String', ...
            sprintf('%s  |  jaw, time-warped on real lick contacts  |  mean \\pm 1.96 SEM  |  %s', hdr, strtrim(nStr)), ...
            'EdgeColor','none', 'HorizontalAlignment','center', ...
            'FontWeight','bold', 'FontSize',10, 'Interpreter','tex');
    end
else
    ax = axes('Parent', wt, 'Position', [0.09 0.15 0.86 0.72]);  hold(ax,'on')
    [Wa, Wp, slotStart] = jawWarpedMatrix(jawNorm(aIdx, testTrials), predFull(aIdx,:), segCell(testTrials), cfg);
    hA = drawWarpedSet(ax, Wa, [0 0 0], slotStart, cfg);
    hP = drawWarpedSet(ax, Wp, predCol,  slotStart, cfg);
    xlabel(ax, 'lick contact');  ylabel(ax, 'Jaw y-displacement (norm)');
    if cfg.warpedTickCenter, tickAt = slotStart + cfg.warpedSegmentLength/2; else, tickAt = slotStart; end
    xticks(ax, tickAt(1:cfg.maxLickShow));
    xticklabels(ax, arrayfun(@(L) sprintf('contact %d', L), 1:cfg.maxLickShow, 'UniformOutput', false));
    xtickangle(ax, 45);
    xlim(ax, [slotStart(1)-5, slotStart(cfg.maxLickShow)+cfg.warpedSegmentLength+5]);
    legend(ax, [hA hP], {'actual','predicted'}, 'Location','northeast','FontSize',8);
    title(ax, {hdr, sprintf('warped jaw cycles 1-%d, TEST | n = %d trials', ...
        cfg.maxLickShow, numel(testTrials))}, 'FontSize',10,'FontWeight','normal');
    box(ax,'off')
end

end   % wantWarp

clear spikes Xtrain Xtest Ptr Pte predFull Mact Mpred Wa Wp
end   % probe
end   % session
%% REPORT
nLA  = cfg.nLicksAnalyze;
nG   = numel(spec.groupMaps);
nTS  = numel(spec.testSets);
lickVec = (1:nLA)';

% Colors come from the study block.
%   spec.tsColors   one row per TEST SET   (R1 vs R4, R1 vs R6, ...)
%   spec.grpColors  one row per GROUP, or [] for the default gradient.
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

for g = 1:nG
    rows = find(arrayfun(@(r) spec.groupMaps{g}(r.session) == r.probe, res));
    out.sessIdx{g} = [res(rows).session];
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
% Paired Wilcoxon signed-rank at every contact from 2 up (tail set by spec.vsC1Tail).
% BH-FDR across contacts 2..nLicks is also computed; significance marks use q.
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
logf('  (* = BH-FDR q < %.2f; the figure stars these same q values)\n', ...
    cfg.alpha);
fprintf('%s\n', repmat('=',1,84));

%% BETWEEN TEST SETS, per contact
% For the two-pool tasks (R1 vs R4, R1 vs R6): paired at each contact, TWO-tailed
% because there is no prior direction.
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

%% THE TEST THE PANEL REPORTS, AND THE SUMMARY FILE
% One value per session, relative to the within-session reference (contact 1):
% so  value(session) = mean over spec.groupSummaryLicks of
%                      ( decIdx(L) - decIdx(1) ).
% Contact 1 is zero by construction, so the absolute level does not enter.
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
% The learning panels (Fig. 5G / S4J) use a within-animal Day 1 vs Day 5 test,
% not a session-level test, so nothing is tested here.

case 'crosstask'
% Fig. 3H / S4F: C1 trials in the Simple Reward Task against C1 trials in
% the Double Reward Task. The two sides come from DIFFERENT SCRIPTS, so the
% Simple Reward script has to have been run first -- it writes its summary
% to spec.summaryFile, and this reads it from spec.crossTaskFile.
% One-tailed rank-sum, one test per region, not corrected.
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

% ---- save the per-session summary (read by the cross-task comparison) ----
summaryByGroup = cell(1, nG);
for g = 1:nG, summaryByGroup{g} = out.summary{g,1}; end
figRef   = spec.figRef;   %#ok<NASGU>
specName = spec.name;   %#ok<NASGU>
sessionsByGroup = out.sessIdx;   %#ok<NASGU>
groupLabels  = spec.groupLabels;        %#ok<NASGU>   % saved so a reader can check the regions
summaryLicks = spec.groupSummaryLicks;  %#ok<NASGU>   % and which contacts were averaged
save(spec.summaryFile, 'summaryByGroup', 'sessionsByGroup', 'groupLabels', ...
     'summaryLicks', 'figRef', 'specName');
logf('summary written to %s\n\n', spec.summaryFile);

%% FIGURE 1: DECODING INDEX, TABBED
if spec.showFigCI
% Two tabs: the mean decoding index with its CI band, and the mean jaw in the fitting window.
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
        finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{0}', h, lbl, [0 1]);
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

% ---- TAB: mean jaw displacement inside the fitting window ----
% Each trial contributes only the samples inside its own (variable-length) window,
% so the trace is averaged over however many trials cover each sample; that count
% is shown on the right axis and thins out at the edges.
        ax = axes('Parent', uitab(tg, 'Title', 'jaw in window'));  hold(ax,'on')
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
            ylabel(ax, 'jaw displacement (norm)', 'FontSize', 11);   % not clamped to [0 1]: jaw is percentile-scaled
            yyaxis(ax,'right')
            plot(ax, tt, mean(N, 2, 'omitnan'), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.25);
            ylabel(ax, 'trials covering the sample', 'FontSize', 10);
            yyaxis(ax,'left')
            good = find(mean(N,2,'omitnan') > 0);
            if ~isempty(good)
                xlim(ax, [tt(good(1))-0.05, tt(good(end))+0.05]);
            end
        end
        yl = ylim(ax);  plot(ax, [0 0], yl, 'r--', 'HandleVisibility','off');  ylim(ax, yl);
        xlabel(ax, 'time from contact 1 (s)', 'FontSize', 11);  box(ax,'off')
        title(ax, {sprintf('%s | %s | mean jaw displacement inside the fitting window', spec.name, figNames{fi}), ...
                   'mean \pm 95% CI across sessions; red line = contact 1 (t = 0)'}, ...
            'FontSize',9,'FontWeight','normal');
    end
end

%% FIGURE 2: MEANS ONLY (no error bars)
% The same means with nothing else drawn, so the group ordering is readable.
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
    finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{0}', h, lbl, [0 1]);
    title(ax, sprintf('%s | mean decoding index', spec.name), ...
        'FontSize',10,'FontWeight','normal');
end

%% FIGURE 3: DELTA FROM CONTACT 1, TABBED
if spec.showFigDelta
% Every session starts at exactly 0 by construction, so the band at contact 1
% collapses to a point. This is the quantity the contact-L-vs-contact-1 test is
% run on, so the stars here are the same test as on Figure 1.
% One figure per group unless spec.overlayGroups.
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
end

%% FIT QUALITY TABLE
logf('\n%s\n', repmat('=',1,110));
logf('%-5s %-8s %-12s %5s %-8s %6s %6s %6s %10s %7s %9s %9s %11s\n', ...
    'sess','animal','date','probe','region','units','nTrn','nTst','lambda','edf', ...
    'R2 train','R2 test','win samp');
logf('%s\n', repmat('-',1,110));
for r = res
    wl = r.winLenSamples;
    logf('%-5d %-8s %-12s %5d %-8s %6d %6d %6d %10.3g %7.1f %9.3f %9.3f %5d-%-5d\n', ...
        r.session, r.anm, r.date, r.probe, r.region, r.nUnits, r.nTrain, r.nTest, ...
        r.lambda, r.edf, r.r2TrSel, r.r2TeSel, min(wl), max(wl));
end
logf('%s\n', repmat('=',1,110));

end   % repIdx -- REPEAT PASSES

%% LOCAL FUNCTIONS

function winIdx = jawFitWindow(segs, aIdx, traceTime, goCueRel_s, contactRel, cfg)
% Sample indices of ONE trial's fitting window, in FULL-TRACE coordinates.
% segs is that trial's cell array from localJawSegments, indexed by CONTACT
% NUMBER and holding sample indices into the ANALYSIS WINDOW, so segs{k} is
% contact k's jaw cycle, segs{k}(end) is where that cycle closes, and aIdx maps
% it back onto the full trace. contactRel holds the post-go-cue contact times
% relative to contact 1.
% 'contact1'        go cue - preGoCue_s  ->  cycle 1 closes + postJawClose_s.
%                   The end is the jaw CLOSING, not the contact: the contact
%                   happens part way up the opening while the jaw is still
%                   moving, so the closing trough is what ends the movement.
% 'bout'            contact 2 + winBoutStartPad_s  ->  bout-end contact +
%                   winBoutEndPad_s. Contact times only, so the window matches
%                   the tongue script's 'bout' window.
% 'lick2ToBoutEnd'  cycle 2 CLOSES + winBoutStartPad_s  ->  same end as 'bout'.
    winIdx = [];
    if ~isfinite(goCueRel_s), return; end
    switch lower(cfg.fitWindowMode)
        case 'contact1'
            if isempty(segs) || isempty(segs{1}), return; end
            tStart = goCueRel_s - cfg.preGoCue_s;
            tEnd   = traceTime(aIdx(segs{1}(end))) + cfg.postJawClose_s;
        case 'bout'
            if numel(contactRel) < 2, return; end
            kEnd = boutEndContact(contactRel, cfg.winBoutGap_s);
            if kEnd < 2, return; end   % the bout is over before contact 2
            tStart = contactRel(2)    + cfg.winBoutStartPad_s;
            tEnd   = contactRel(kEnd) + cfg.winBoutEndPad_s;
        case 'lick2toboutend'
            if numel(contactRel) < 2, return; end
            kEnd = boutEndContact(contactRel, cfg.winBoutGap_s);
            if kEnd < 2, return; end
            if numel(segs) < 2 || isempty(segs{2}), return; end   % no cycle 2
            tStart = traceTime(aIdx(segs{2}(end))) + cfg.winBoutStartPad_s;
            tEnd   = contactRel(kEnd)              + cfg.winBoutEndPad_s;
        otherwise
            error('cfg.fitWindowMode must be ''contact1'', ''bout'' or ''lick2ToBoutEnd''.');
    end
    if tEnd <= tStart, return; end
    w = find(traceTime >= tStart & traceTime <= tEnd);
    if numel(w) >= cfg.minWindowSamples, winIdx = w; end
end

function k = boutEndContact(contactRel, gap_s)
% Which contact ENDS the bout: scanning forward from contact 1, the first one
% whose gap to the NEXT contact exceeds gap_s. If no gap ever exceeds it, the
% bout runs to the last contact of the trial.
% This is a stopping rule on the continuous train, so a stray lick seconds after
% the bout ended is not swept in.
    d = diff(contactRel(:));
    k = find(d > gap_s, 1);
    if isempty(k), k = numel(contactRel); end
end

function [decIdx, peakL, nTrialL] = jawDecodingIndex(Yact, Ypred, segs, cfg)
% Per-contact index 1 - RMSE_dec/RMSE_0 over the samples of each contact's cycle, pooled across trials.
% RMSE_0 is the error made by predicting a scaled jaw of 0 and never moving --
% the same null the tongue scripts use. It is a CONSTANT, identical for every
% trial and every session, so nothing about the denominator varies with the
% animal's posture or with how the bout ended.
% A fixed baseline does shrink when the excursion shrinks: a
% purely behavioral decline in jaw opening across the bout lowers the index on
% its own, with no change in how well the neural data predicts. peakL, the
% per-contact peak above 0, is returned alongside so that can be checked.
% One finite mask, shared by actual and predicted, is used for both RMSEs, so the
% index cannot rise with the tracking-failure rate.
% Cycles come in precomputed from localJawSegments -- the same segmentation the
% fitting window and the warped figure use.
    nLA = cfg.nLicksAnalyze;
    nTr = size(Yact, 2);
    sseDec = zeros(nLA,1);  sseBase = zeros(nLA,1);
    nTrialL = zeros(nLA,1); peaks = nan(nLA, nTr);
    for tt = 1:nTr
        sg = segs{tt};
        if isempty(sg), continue; end
        for L = 1:min(numel(sg), nLA)
            seg = sg{L};
            if isempty(seg), continue; end
            a = Yact(seg,tt);  p = Ypred(seg,tt);
            g = isfinite(a) & isfinite(p);
            if ~any(g), continue; end
            sseDec(L)  = sseDec(L)  + sum((a(g) - p(g)).^2);
            sseBase(L) = sseBase(L) + sum(a(g).^2);   % baseline = 0
            nTrialL(L) = nTrialL(L) + 1;
            peaks(L,tt) = max(a(g));   % excursion above 0
        end
    end
    decIdx = 1 - sqrt(sseDec ./ sseBase);   % n cancels inside the ratio
    decIdx(nTrialL == 0 | sseBase == 0) = NaN;
    peakL = mean(peaks, 2, 'omitnan');
end

function [Wa, Wp, slotStart, traceLen] = jawWarpedMatrix(Yact, Ypred, segs, cfg)
% THE WARPED LICK TRAIN, same slot geometry the tongue scripts use, fed by the
% jaw CYCLES instead of tongue protrusions.
% Every cycle is resampled onto cfg.warpedSegmentLength points and dropped into
% its own slot, slotStart = warpedSlotFirst : warpedSlotSpacing : ... . The
% spacing is wider than the segment so the gaps stay NaN and the licks read as
% separate events. Predictions are warped through the SAME sample indices as the
% actual trace, never re-detected, so the two curves align cycle for cycle.
% WHY WARP. Later cycles in a bout are shorter as well as shallower. On a real
% time axis that smears them together and flattens the mean for reasons that have
% nothing to do with decoding; warping takes duration out so what is left to
% compare is amplitude and shape.
    nEx  = cfg.maxLicksToExtract;
    segL = cfg.warpedSegmentLength;
    slotStart = cfg.warpedSlotFirst : cfg.warpedSlotSpacing : ...
                (cfg.warpedSlotFirst + cfg.warpedSlotSpacing*(nEx-1));
    traceLen  = max(slotStart) + segL - 1;
    nTr = size(Yact, 2);
    Wa  = nan(traceLen, nTr);
    Wp  = nan(traceLen, nTr);
    for tt = 1:nTr
        sg = segs{tt};
        if isempty(sg), continue; end
        for L = 1:min(numel(sg), nEx)
            seg = sg{L};
            if numel(seg) < 2, continue; end
            dest = slotStart(L) : slotStart(L) + segL - 1;
            Wa(dest,tt) = warpSegment(Yact(seg,tt),  segL, cfg.warpedInterp);
            Wp(dest,tt) = warpSegment(Ypred(seg,tt), segL, cfg.warpedInterp);
        end
    end
end

function [nRun, nInWin] = longestLickRun(lickTimes, goCue, cfg)
% Longest run of CONSECUTIVE port contacts that are all inside cfg.boutWin_s of
% the go cue and no more than cfg.boutMaxILI_s apart. Identical to the tongue
% scripts -- it reads lick CONTACTS, not kinematics, so it does not care which
% variable is being decoded.
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

function s = groupOfProbe(spec, sessionIdx, probeIdx)
% Which GROUP this probe carries in THIS session -- the region, for the region
% scripts, and the learning day for the learning script.
% Read off the group maps rather than assumed, because a probe index means
% different things in different sessions: probe 2 is ALM in one session and M1 in
% the next. Normally exactly one group points here; if two maps overlap they are
% joined with '+', and if none points here it says so rather than guessing.
    hit = cellfun(@(m) m(sessionIdx) == probeIdx, spec.groupMaps);
    if ~any(hit)
        s = 'no group';
    else
        s = strjoin(spec.groupLabels(hit), '+');
    end
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
% Lagged design from A (nT x nUnits x nTrials) over PER-TRIAL sample sets.
% ROW ORDER: trials in order, samples within a trial in order -- matching
% stackWindows above.
% COLUMN ORDER: lag-major, unit-minor; column (li-1)*nUnits + u is unit u at lag lags(li).
% Each trial is embedded against ONLY its own samples, so no lag ever crosses a
% trial boundary. A lag that would leave the trial gives NaN and that row is
% dropped before fitting.
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
% cfg.ridgeImpl = 'matlab' calls MATLAB's ridge; 'eig' uses the algebraically
% identical closed form below, which is far faster with many predictors. Their
% agreement can be checked at runtime (cfg.ridgeCheckEquivalence).
% ---- ridge(y, X, k, 0) ----------------------------------------------------
% The trailing 0 asks for coefficients ON THE SCALE OF THE ORIGINAL PREDICTORS,
% and the returned matrix is (p+1) x numel(lambdas) with the INTERCEPT IN ROW 1.
% Without the trailing 0, ridge returns standardized-scale coefficients and no intercept.
% ---- WHAT IS PENALIZED --------------------------------------------------
% ridge always penalizes the STANDARDIZED coefficients and always leaves the
% intercept unpenalized. That is the same estimator as the closed form below
% with standardizeCols = true, so standardizeCols = false has no counterpart in
% ridge and is rejected rather than quietly ignored.
% ---- EFFECTIVE DEGREES OF FREEDOM ---------------------------------------
% R.edf(k) = 1 + sum_i d_i/(d_i + lambda_k), over the eigenvalues d of the
% standardized Gram matrix; the 1 is the free intercept. It is the trace of the
% hat matrix. ridge does not return it, so it is computed from the singular
% values of the centered, scaled design -- ONE svd, whatever the number of
% lambdas. It is needed only by the fit printout,
% so the cross-validation folds pass wantEdf = false and skip it.
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
% Algebraically identical to MATLAB's ridge, from ONE
% eigendecomposition of the Gram matrix rather than one augmented
% least-squares solve per lambda. ridge() forms [Z; sqrt(k)*I] and
% calls mldivide separately for EVERY k, so a 50-point grid with p
% in the thousands is 50 QR factorisations of an (n+p) x p matrix
% per path -- and there is one path per fit plus one per CV fold.
% This is one eig and then the whole grid as a single product.
% Checked against ridge(): the same lambda is chosen and predictions agree to ~1e-14.
            mx = mean(X,1);  my = mean(y);
            Xc = X - mx;     yc = y - my;
            if standardizeCols
% Standard deviations below sqrt(eps) are set to 1, the same threshold ridge()
% uses, so near-constant columns are scaled identically by both routes.
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
% The edf below needs exactly these eigenvalues. Hand them over
% rather than recomputing them with an svd of the full n x p
% design, which costs about as much as the entire fit.
            dGram  = d;

        otherwise
            error('cfg.ridgeImpl must be ''matlab'' or ''eig''; got ''%s''.', impl);
    end

    if wantEdf
        if ~isempty(dGram)
% 'eig' already has them: the squared singular values of the
% centered, scaled design ARE the Gram eigenvalues. The svd below
% recomputes the same p numbers from the full n x p matrix at
% roughly the cost of the whole fit, and allocates another copy
% of X to do it.
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
% Fit the same data both ways once and report the largest relative difference.
% This exists so that switching cfg.ridgeImpl to 'eig' for speed is a checkable
% claim rather than a promise. The two routes are algebraically the same
% estimator; in practice they agree to around 1e-11 relative, which is round-off
% at double precision. If this prints something larger, use cfg.ridgeImpl =
% 'matlab' for that data.
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

function wt = warpSegment(v, segL, method)
% One excursion resampled onto segL points (a single finite sample goes to the middle).
    wt = nan(segL,1);
    xv = v(isfinite(v));
    if numel(xv) == 1
        wt(ceil(segL/2)) = xv;
    elseif numel(xv) > 1
        wt = interp1(linspace(0,1,numel(xv))', xv(:), linspace(0,1,segL)', method, NaN);
    end
end

function h = drawWarpedSet(ax, W, col, slotStart, cfg)
% Draw one warped set (actual, or predicted) on ax: mean +/- 1.96 SEM, drawn
% SEGMENT BY SEGMENT so no line is drawn across the gap between two licks.
% Uses the normal 1.96 rather than a t quantile.
    segL = cfg.warpedSegmentLength;
    x  = (1:size(W,1))';
    nF = sum(isfinite(W), 2);
    mW = mean(W, 2, 'omitnan');
    sW = std(W, 0, 2, 'omitnan') ./ sqrt(max(nF,1));
    sW(nF < 2) = NaN;
    loW = mW - 1.96*sW;   hiW = mW + 1.96*sW;
    h = gobjects(1);
    for L = 1:min(cfg.maxLickShow, numel(slotStart))
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

function drawWarpedPanel(ax, gridT, Wa, Wp, medTemplate, nShow, colP, nameStr, showX, cfg)
% ONE CONDITION ON ONE AXES: actual in black, predicted in that condition's
% color, labeled with text on the axes rather than a legend box.
% The labels sit upper-right because that is where the trace has room once the
% bout decays; if a session's licks stay large to the end they may land on data,
% and the two text() lines below are the place to move them.
    hold(ax, 'on')
    localBandPlot(ax, gridT, Wa, [0 0 0], 2.5);
    localBandPlot(ax, gridT, Wp, colP,   2.0);
    axis(ax, 'tight');
    yl = ylim(ax);
    for lk = 1:nShow   % dotted marker at each template contact
        line(ax, [medTemplate(lk) medTemplate(lk)], yl, 'Color',[0.5 0.5 0.5], ...
            'LineStyle',':', 'HandleVisibility','off');
    end
    ylim(ax, yl);
    xr = [medTemplate(1)-2*cfg.warpDelta_s, medTemplate(nShow)+1.5*cfg.warpDelta_s];
    xlim(ax, xr);
    xticks(ax, medTemplate(1:nShow));
    if showX
        xticklabels(ax, arrayfun(@(L) sprintf('%d', L), 1:nShow, 'UniformOutput', false));
        xlabel(ax, 'lick contact (warped, uniform spacing)', 'FontSize', 10);
    else
        xticklabels(ax, []);   % only the bottom panel carries the axis
    end
    ylabel(ax, 'Jaw y-disp. (norm)', 'FontSize', 10);
    tx = xr(1) + 0.60*diff(xr);
    text(ax, tx, yl(1) + 0.94*diff(yl), 'Actual', 'Color', [0 0 0], ...
        'FontWeight','bold', 'FontSize',10, 'Clipping','off');
    text(ax, tx, yl(1) + 0.80*diff(yl), sprintf('%s predicted', nameStr), 'Color', colP, ...
        'FontWeight','bold', 'FontSize',10, 'Clipping','off');
    box(ax, 'off')
end

function ax = gridAxes(parent, nRow, nCol, k)
% Axes at the position subplot(nRow,nCol,k) would use, parented to any container
% (subplot only targets a figure, so it cannot place axes inside a uitab).
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
% Color by GROUP when groups share an axes (the learning script), otherwise by
% test set (R1 vs R4 and friends, one figure per region).
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
% Shared formatting for per-contact axes: zero line, contact ticks, labels, legend.
% ylimVal = [] leaves the y limits automatic; [0 1] is used for the decoding index.
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
% Returns a SLIM struct holding only the fields the decoder actually reads:
%   S.time, S.trialdat, S.nUnitsTotal, S.jawRaw, S.cluid, S.trialid,
%   S.anm, S.date, S.Ntrials, S.lickL, S.goCue, S.trialTypes, S.hit
% Loading dominates run time and is deterministic given (loader, date, params),
% so the result is cached; with the cache warm, a re-run skips straight to the fit.
% CACHE KEY. The full set of params fields that can change what gets loaded is
% written into the .mat and compared on read, so a stale cache can never be used
% silently: changing params.smooth, params.quality or the condition list reloads.
% trialdat is stored as SINGLE. It is the bulk of the file (nTime x nUnits x
% nTrials) and single halves it; ridgePath casts back to double before forming
% the Gram matrix, so nothing downstream loses precision.

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

featCol = find(strcmp(kin.featLeg, cfg.jawFeature), 1);
assert(~isempty(featCol), '%s not found in kin.featLeg for %s', cfg.jawFeature, dateStr);

S.time        = obj.time(:);
S.trialdat    = single(obj.trialdat);
S.nUnitsTotal = size(obj.psth, 2);
S.jawRaw      = kin.dat(:,:,featCol);
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
k.feature = 'jaw';   % keeps jaw and tongue caches from colliding
k.source  = 'slimExport';   % never reuse a cache entry built by the pipeline
key = jsonencode(k);
end

function [lickRuns, nDisagree, segDiag] = localJawSegments(obsMat, testTrialNums, goCueAll, lickLAll, timeSub, cfg)
% Jaw lick segmentation: one trough-to-trough cycle per lick contact, returned as
% sample indices into the analysis window. cfg.segAnchor selects the method:
%   'persistCycle'  persistence simplification of the extremum sequence (default;
%                   see the branch below).
%   'contactCycle'  contact-anchored:
%   1. ONE PEAK PER CONTACT. For contact k, take the MAXIMUM of the smoothed
%      trace over the window bounded by the midpoints to contacts k-1 and k+1
%      (scaled by cfg.peakSearchFrac). A maximum over a window is not a local
%      feature -- there is exactly one and it is the highest point -- so a
%      plateau wiggle can never be mistaken for the peak. Midpoint bounds also
%      make the windows disjoint, so peaks come out strictly ordered and two
%      contacts can never claim the same excursion.
%   2. TROUGHS BY PROMINENCE WALK, not by argmin of a window. Step outward from
%      the peak tracking the running minimum; stop when the trace has climbed
%      back above that minimum by more than max(cfg.troughRiseFrac x depth,
%      cfg.troughRiseFloor). The running minimum at that moment is the trough.
%      A plateau wiggle cannot stop the walk (near the peak the depth is small,
%      so the absolute floor governs and the wiggle is under it); a genuine
%      climb into the next cycle always does.
%   3. BROAD BASINS. If the bottom is a flat stretch rather than a sharp V, the
%      boundary moves to the end of the CONTIGUOUS run of samples within
%      cfg.troughFlatTol of the trough -- the moment opening begins on the way
%      in, bottom reached on the way out. The default tolerance of 0 takes the
%      actual trough; wiggle protection comes from the walk.
%   4. NUMBERING IS EXACT BY CONSTRUCTION. Segment k IS contact k; there is no
%      peak-to-contact matching step, so no lick can be skipped or mislabeled.
%   5. AMPLITUDE GATE. peak minus the deeper of its two troughs must reach
%      cfg.minCycleAmp, so a contact with no real jaw excursion (dropout, or a
%      lick made with almost no opening) does not get a window drawn around
%      noise. Set cfg.minCycleAmp = 0 to segment every contact regardless.
% findpeaks still runs, but under 'contactCycle' it only supplies a count for
% the audit and decides no boundary.
%   'peakDetect'    findpeaks-based (uses cfg.jawLickIndexRule, cfg.requireTroughBelowMax and cfg.segMode).
% Segments are one full cycle (trough -> peak -> trough), so their length varies
% between licks and trials.
% NaN handling: findpeaks and min/max cannot be trusted across
% NaN, so detection runs on a gap-filled COPY and the returned indices are
% applied to the real NaN-bearing trace by the caller. No interpolated sample is
% ever scored.
% OUTPUTS
%   nDisagree  trials where contact-anchored windows differ from what plain
%              peak-order would have produced (a scale for how much the anchor
%              matters on this data)
%   segDiag    nTrials x 6: [nContacts nSegmented nPeaksFindpeaks nFailedAmp nEdgeTruncated nRescued]

% ---- defaults for any cfg field not set ----
    if ~isfield(cfg,'segAnchor'),              cfg.segAnchor = 'persistCycle'; end
    if ~isfield(cfg,'cyclePersistFrac'),       cfg.cyclePersistFrac = 0.3; end
    if ~isfield(cfg,'cyclePersistFloor'),      cfg.cyclePersistFloor = 0.03; end
    if ~isfield(cfg,'rescueMissedCycles'),     cfg.rescueMissedCycles = true; end
    if ~isfield(cfg,'rescueMaxSpanILIFrac'),   cfg.rescueMaxSpanILIFrac = 2.2; end
    if ~isfield(cfg,'rescueMinAmpFrac'),       cfg.rescueMinAmpFrac = 0.12; end
    if ~isfield(cfg,'rescueEdgeContacts'),     cfg.rescueEdgeContacts = true; end
    if ~isfield(cfg,'contactMatchILIFrac'),    cfg.contactMatchILIFrac = 0.7; end
    if ~isfield(cfg,'cycleMaxHalfILIFrac'),    cfg.cycleMaxHalfILIFrac = 0.75; end
    if ~isfield(cfg,'segMode'),                cfg.segMode = 'troughToTrough'; end
    if ~isfield(cfg,'requireTroughBelowMax'),  cfg.requireTroughBelowMax = false; end
    if ~isfield(cfg,'segMinLenSamples'),       cfg.segMinLenSamples = 5; end
    if ~isfield(cfg,'segMaxLenSamples'),       cfg.segMaxLenSamples = Inf; end
    if ~isfield(cfg,'peakSearchFrac'),         cfg.peakSearchFrac = 0.5; end
    if ~isfield(cfg,'segMaxHalfCycleFrac'),    cfg.segMaxHalfCycleFrac = 1.5; end
    if ~isfield(cfg,'troughRiseFrac'),         cfg.troughRiseFrac = 0.25; end
    if ~isfield(cfg,'troughRiseFloor'),        cfg.troughRiseFloor = 0.01; end
    if ~isfield(cfg,'troughMinDepthFrac'),     cfg.troughMinDepthFrac = 0.6; end
    if ~isfield(cfg,'allowEdgeTruncated'),     cfg.allowEdgeTruncated = true; end
    if ~isfield(cfg,'troughFlatTol'),          cfg.troughFlatTol = 0; end
    if ~isfield(cfg,'minCycleAmp'),            cfg.minCycleAmp = 0; end

    [nSamp, nTrials] = size(obsMat);
    lickRuns  = cell(nTrials,1);
    nDisagree = 0;
    segDiag   = zeros(nTrials, 6);   % [nContacts nSegmented nPeaksFindpeaks nFailedAmp nEdgeTruncated nRescued]

    timeSub = timeSub(:);
    dtSub   = median(diff(timeSub));
    [~, minLimIdx] = min(abs(timeSub - cfg.segMinTime_s));

    for tr = 1:nTrials
        trialIdx = testTrialNums(tr);
        gc = goCueAll(trialIdx);
        licks = lickLAll{trialIdx, 1};
        if isempty(licks) || all(isnan(licks)), continue; end
        licks = sort(licks(licks > gc));
        if isempty(licks), continue; end
        contactRel = licks(:) - licks(1);   % traceTime == 0 is contact 1

        raw = obsMat(:,tr);
        if all(isnan(raw)), continue; end
        filled = fillmissing(raw, 'linear', 'EndValues','nearest');   % detection copy only
        sig    = movmean(filled, cfg.smoothSamples);

% findpeaks: boundary-setting under 'peakDetect', reporting under 'contactCycle'
        [~, fpLocs] = findpeaks(sig, 'MinPeakProminence',cfg.peakProminence, ...
            'MinPeakDistance',cfg.peakMinDist, 'MinPeakHeight',cfg.peakMinHeight);
        fpLocs = fpLocs(:);
        nFP    = numel(fpLocs);

        runs       = cell(cfg.numLicksToAnalyze,1);
        nAmpFail   = 0;
        nEdgeTrunc = 0;
        nRescued   = 0;

        switch lower(cfg.segAnchor)

        case 'persistcycle'
% PERSISTENCE-BASED CYCLE SEGMENTATION.
% A plateau wiggle and a real trough have the same local shape and differ only
% in scale, which is a global property of the waveform. So the full alternating
% max/min sequence of the trace is built, and the least prominent adjacent
% extremum PAIR is repeatedly deleted until every remaining extremum has
% prominence at least persistThresh (topological persistence simplification):
%   - a plateau wiggle is a max/min pair separated by a tiny amplitude.
%     It has near-zero persistence and is merged away early, so it can
%     never become a boundary.
%   - a real lick cycle has persistence on the order of the jaw's
%     working range and always survives.
%   - the threshold is one interpretable number: a fraction of the
%     trial's own robust jaw range. It says "an excursion smaller than
%     this fraction of the range is not a lick cycle", which is a
%     statement about jaw physiology, not about detector tuning.
% The surviving minima are the cycle boundaries, the surviving maxima
% are the peaks, and between consecutive minima there is exactly one
% peak by construction. Cycles are then assigned to lick contacts
% CONTAINMENT FIRST (a contact falling inside a cycle gets that cycle),
% then nearest-peak for whatever is left. Containment first is what
% stops adjacent contacts from being skipped when the licking is fast
% and cfg.contactMatchTol_s would otherwise be the binding constraint.
% Trace endpoints count as boundaries, so lick 1 keeps its cycle when
% its opening trough is off the left edge of the analysis window.

% ---- STEP 1: alternating extremum sequence over the whole trace ----
            dsig = diff(sig);
            sTrend = sign(dsig);
% carry the last nonzero slope across flat stretches, so a plateau
% does not register as a pile of spurious extrema
            lastS = 0;
            for i = 1:numel(sTrend)
                if sTrend(i) == 0
                    sTrend(i) = lastS;
                else
                    lastS = sTrend(i);
                end
            end
            firstNZ = find(sTrend ~= 0, 1);
            if isempty(firstNZ)
                segDiag(tr,:) = [numel(contactRel), 0, nFP, 0, 0, 0];   % flat trace
                continue
            end
            sTrend(1:firstNZ-1) = sTrend(firstNZ);

            E = unique([1; find(diff(sTrend) ~= 0) + 1; nSamp]);
            if numel(E) < 3
                segDiag(tr,:) = [numel(contactRel), 0, nFP, 0, 0, 0];
                continue
            end

% classify each as max or min
            isMax = false(numel(E),1);
            for i = 1:numel(E)
                if i == 1
                    isMax(i) = sig(E(1)) > sig(E(2));
                elseif i == numel(E)
                    isMax(i) = sig(E(end)) > sig(E(end-1));
                else
                    isMax(i) = sig(E(i)) >= sig(E(i-1)) && sig(E(i)) >= sig(E(i+1));
                end
            end

% enforce strict alternation: within any run of same-type extrema
% keep only the most extreme one
            keepE = true(numel(E),1);
            i = 1;
            while i <= numel(E)
                j = i;
                while j+1 <= numel(E) && isMax(j+1) == isMax(i), j = j + 1; end
                if j > i
                    blockIdx = E(i:j);
                    if isMax(i)
                        [~, b] = max(sig(blockIdx));
                    else
                        [~, b] = min(sig(blockIdx));
                    end
                    keepE(i:j) = false;
                    keepE(i + b - 1) = true;
                end
                i = j + 1;
            end
            E = E(keepE); isMax = isMax(keepE);

% ---- STEP 2: persistence simplification ----
% Robust range of THIS trial, so the threshold scales with how far
% the animal actually moves its jaw rather than being a fixed number
% that is too strict on one session and too loose on another.
            rngLo = prctile(sig, 5); rngHi = prctile(sig, 95);
            jawRange = rngHi - rngLo;
            if ~isfinite(jawRange) || jawRange <= 0
                jawRange = max(sig) - min(sig);
            end
            persistThresh = max(cfg.cyclePersistFrac * jawRange, cfg.cyclePersistFloor);

            while numel(E) >= 3
                v  = sig(E);
                d1 = abs(v(2:end-1) - v(1:end-2));
                d2 = abs(v(2:end-1) - v(3:end));
                pr = min(d1, d2);   % persistence of each interior extremum

% ---- ONE-SIDED PERSISTENCE AT THE WINDOW EDGES ----
% Persistence normally takes the SMALLER
% of the two drops on either side of an extremum, which is right
% when both neighbors are real troughs. The first element of E
% is not a real trough -- it is wherever the analysis window
% happened to be cut (t = -0.09). When lick 1 has no trough
% before it inside the window, the jaw is already part-way open
% at that cut, so the "drop" on the left is small and would
% wrongly give lick 1 low persistence.
% So at a window edge, persistence is measured on the side that
% has a real trough, and only that side.
                if E(1) == 1 && ~isempty(pr)
                    pr(1) = d2(1);   % left neighbor is the window start
                end
                if E(end) == nSamp && ~isempty(pr)
                    pr(end) = d1(end);   % right neighbor is the window end
                end

                [mn, jj] = min(pr);
                if mn >= persistThresh, break; end
                j = jj + 1;   % index into E

% Remove the low-persistence extremum together with the
% neighbor it merges into. Removing a PAIR is what keeps the
% sequence alternating; removing one alone would leave two
% maxima adjacent and corrupt every later step.
                takeLeft = (d1(jj) <= d2(jj));
% ...but never delete a window edge. It is the boundary that
% lets an incomplete first or last cycle still have a window,
% and merging it away would drop lick 1.
                if takeLeft && (j-1) == 1 && E(1) == 1
                    if (j+1) > numel(E), break; end
                    takeLeft = false;
                end
                if ~takeLeft && (j+1) == numel(E) && E(end) == nSamp
                    if (j-1) < 1, break; end
                    takeLeft = true;
                end
                if takeLeft
                    rm = [j-1, j];
                else
                    rm = [j, j+1];
                end
                if rm(1) < 1 || rm(2) > numel(E), break; end
                E(rm) = []; isMax(rm) = [];
            end

            minIdx = E(~isMax);
            maxIdx = E(isMax);
            if isempty(maxIdx)
                segDiag(tr,:) = [numel(contactRel), 0, nFP, 0, 0, 0];
                continue
            end

% trace endpoints act as boundaries, so a first or last cycle whose
% trough is off-screen still gets a window
            bnd = minIdx(:);
            if isempty(bnd) || maxIdx(1)   < bnd(1),   bnd = [1; bnd];      end
            if isempty(bnd) || maxIdx(end) > bnd(end), bnd = [bnd; nSamp];  end
            bnd = unique(bnd);

% this trial's own lick rate: sets how long a cycle is allowed to be
            if numel(contactRel) >= 2
                iliHere = median(diff(contactRel));
            else
                iliHere = cfg.contactMatchTol_s;
            end
            if ~isfinite(iliHere) || iliHere <= 0, iliHere = cfg.contactMatchTol_s; end
            maxHalfSamp = max(2, round(cfg.cycleMaxHalfILIFrac * iliHere / dtSub));

% ---- STEP 3: one cycle per surviving peak ----
            nM = numel(maxIdx);
            cycS = nan(nM,1); cycE = nan(nM,1);
            for i = 1:nM
% STRICTLY before / strictly after, with the window edge as the
% fallback when there is no boundary on that side.
% A peak sitting exactly on the left window edge must not take
% bnd = 1 on both sides, or its cycle would span the whole window.
                b0 = bnd(find(bnd < maxIdx(i), 1, 'last'));
                if isempty(b0), b0 = 1;     end   % no trough before -> window start
                b1 = bnd(find(bnd > maxIdx(i), 1, 'first'));
                if isempty(b1), b1 = nSamp; end   % no trough after  -> window end
                if b1 <= b0, continue; end

% ---- TRIM LONG TAILS ----
% A cycle runs from one trough to the next, which is right while
% the animal is licking. After the LAST lick of a bout the jaw
% closes and simply stays closed, so the next trough can be
% seconds away and the final cycle grows a long flat tail that
% would dominate the cycle.
% A lick cannot last much longer than the animal's own lick
% period, so that is the bound: at most
% cfg.cycleMaxHalfILIFrac x the trial's median inter-lick
% interval on each side of the peak. When the real trough is
% further out than that, the boundary becomes the LOWEST point
% within the allowance rather than a flat cutoff.
                if maxIdx(i) - b0 > maxHalfSamp
                    segw = (maxIdx(i) - maxHalfSamp):maxIdx(i);
                    [~, rTrim] = min(sig(segw));
                    b0 = segw(rTrim);
                end
                if b1 - maxIdx(i) > maxHalfSamp
                    segw = maxIdx(i):(maxIdx(i) + maxHalfSamp);
                    [~, rTrim] = min(sig(segw));
                    b1 = segw(rTrim);
                end
                if b1 <= b0, continue; end

                cycS(i) = b0; cycE(i) = b1;
            end

            okCyc = isfinite(cycS) & isfinite(cycE) & ...
                    maxIdx(:) >= minLimIdx & ...
                    (cycE - cycS + 1) >= cfg.segMinLenSamples & ...
                    (cycE - cycS + 1) <= cfg.segMaxLenSamples;

% amplitude check, ignoring any side that is a window edge rather than a trough
            for i = find(okCyc(:).')
                aB = sig(maxIdx(i)) - sig(cycS(i));
                aF = sig(maxIdx(i)) - sig(cycE(i));
                edgeS = (cycS(i) == 1); edgeF = (cycE(i) == nSamp);
                if edgeS && ~edgeF
                    ampI = aF;
                elseif edgeF && ~edgeS
                    ampI = aB;
                else
                    ampI = min(aB, aF);
                end
                if ampI < cfg.minCycleAmp
                    okCyc(i) = false; nAmpFail = nAmpFail + 1;
                end
                if edgeS || edgeF, nEdgeTrunc = nEdgeTrunc + 1; end
            end

            kMax = maxIdx(okCyc); kS = cycS(okCyc); kE = cycE(okCyc);
            nK = numel(kMax);
            nCc = min(numel(contactRel), cfg.numLicksToAnalyze);
            if nK == 0 || nCc == 0
                segDiag(tr,:) = [numel(contactRel), 0, nFP, nAmpFail, nEdgeTrunc, nRescued];
                continue
            end

% ---- STEP 4: assign cycles to contacts, CONTAINMENT FIRST ----
            cT = contactRel(1:nCc);
            pkT = timeSub(kMax);  pkT = pkT(:).';
            s0T = timeSub(kS);    s0T = s0T(:).';
            s1T = timeSub(kE);    s1T = s1T(:).';

            D  = abs(repmat(cT, 1, nK) - repmat(pkT, nCc, 1));
            IN = repmat(cT, 1, nK) >= repmat(s0T, nCc, 1) & ...
                 repmat(cT, 1, nK) <= repmat(s1T, nCc, 1);

            usedC = false(nCc,1); usedK = false(nK,1);

% pass 1: contacts that fall INSIDE a cycle take that cycle. Closest
% to the peak commits first, so if two contacts sit in one cycle the
% better-matched one wins and the other goes to pass 2.
            Dp = D; Dp(~IN) = Inf;
            while true
                [mv, li] = min(Dp(:));
                if ~isfinite(mv), break; end
                [ci, pj] = ind2sub(size(Dp), li);
                runs{ci} = (kS(pj):kE(pj)).';
                usedC(ci) = true; usedK(pj) = true;
                Dp(ci,:) = Inf; Dp(:,pj) = Inf;
            end

% pass 2: whatever is left, nearest peak within tolerance. The
% tolerance also scales with the trial's own lick rate (iliHere,
% computed above), so a fast train is not held to a tolerance
% calibrated on a slow one.
            matchTol = max(cfg.contactMatchTol_s, cfg.contactMatchILIFrac * iliHere);

            Dq = D; Dq(usedC, :) = Inf; Dq(:, usedK) = Inf; Dq(Dq > matchTol) = Inf;
            while true
                [mv, li] = min(Dq(:));
                if ~isfinite(mv), break; end
                [ci, pj] = ind2sub(size(Dq), li);
                runs{ci} = (kS(pj):kE(pj)).';
                Dq(ci,:) = Inf; Dq(:,pj) = Inf;
            end

% ---- PASS 3: RESCUE A CONTACT THAT SITS BETWEEN TWO GOOD CYCLES ----
% cfg.cyclePersistFrac is a per-trial threshold on the WHOLE trial's jaw range,
% so a lick whose opening is genuinely smaller than its neighbors can fall
% under it and be merged away. If contact k has no cycle but contacts k-1 and
% k+1 both do, the trace between k-1's closing trough and k+1's opening trough
% is exactly one inter-lick interval with contact k inside it, so the excursion
% in there is lick k's. The rescue only has to find the peak inside that span;
% no global threshold is loosened for the rest of the trial.
% GUARDS, so this cannot invent cycles out of flat trace:
%   - the span must be no wider than cfg.rescueMaxSpanILIFrac x the
%     trial's own inter-lick interval. A wide span means a contact is
%     missing on one side too and the bracket is not one cycle.
%   - the contact itself must fall inside the span.
%   - the peak must clear BOTH ends by cfg.rescueMinAmpFrac of the
%     trial's jaw range. A flat stretch has no peak and is refused.
%   - the peak must be strictly interior; a maximum sitting on a span
%     edge means the excursion belongs to the neighbor, not here.
% EDGE CONTACTS (cfg.rescueEdgeContacts). Lick 1 has no left
% neighbor and the last lick has no right one, so the window edge
% stands in for the missing boundary (e.g. lick 1's opening is off-screen or
% shallow, and the span from the window start to lick 2's opening trough
% contains exactly one excursion).
% TWO PASSES, because a rescued cycle becomes a valid anchor for the
% next gap: with 3 and 4 both missing, rescuing 3 against 2 and 5
% gives pass 2 the boundary it needs for 4.
            if cfg.rescueMissedCycles
              maxSpanSamp = round(cfg.rescueMaxSpanILIFrac * iliHere / dtSub);
              for rescuePass = 1:2
                for k = 1:nCc
                  if ~isempty(runs{k}), continue; end

% nearest assigned contact on each side
                  kl = find(~cellfun(@isempty, runs(1:k-1)), 1, 'last');
                  krRel = find(~cellfun(@isempty, runs(k+1:nCc)), 1, 'first');

                  if ~isempty(kl)
                      spanLo = runs{kl}(end);
                  elseif cfg.rescueEdgeContacts
                      spanLo = 1;   % window start stands in
                  else
                      continue
                  end
                  if ~isempty(krRel)
                      spanHi = runs{k + krRel}(1);
                  elseif cfg.rescueEdgeContacts
                      spanHi = nSamp;   % window end stands in
                  else
                      continue
                  end
                  if spanHi <= spanLo, continue; end
                  if (spanHi - spanLo + 1) < cfg.segMinLenSamples, continue; end
                  if (spanHi - spanLo) > maxSpanSamp, continue; end

% the contact has to be inside the bracket for this to be its cycle
                  ciSamp = round(interp1(timeSub, (1:nSamp).', contactRel(k), 'linear','extrap'));
                  ciSamp = max(1, min(nSamp, ciSamp));
                  if ciSamp < spanLo || ciSamp > spanHi, continue; end

% peak inside the span, but not further from the contact than
% a half interval -- keeps a shoulder at the far end of a wide
% bracket from being taken as this lick's opening
                  nearLo = max(spanLo, ciSamp - round(0.6*iliHere/dtSub));
                  nearHi = min(spanHi, ciSamp + round(0.6*iliHere/dtSub));
                  if nearHi <= nearLo, nearLo = spanLo; nearHi = spanHi; end
                  segw = nearLo:nearHi;
                  [~, relPk] = max(sig(segw));
                  pkHere = segw(relPk);
                  if pkHere <= spanLo || pkHere >= spanHi, continue; end

% Troughs: the lowest point on each side of the peak within the span. Usually
% that is the span edge itself: the previous cycle's close is this cycle's
% open.
                  [~, rl] = min(sig(spanLo:pkHere)); loHere = spanLo + rl - 1;
                  [~, rr] = min(sig(pkHere:spanHi)); hiHere = pkHere + rr - 1;
                  if hiHere <= loHere, continue; end

                  ampHere = sig(pkHere) - max(sig(loHere), sig(hiHere));
                  if ampHere < cfg.rescueMinAmpFrac * jawRange, continue; end

% same tail trim the normal path uses, so a rescued cycle
% cannot be longer than a real one
                  if pkHere - loHere > maxHalfSamp
                      sw = (pkHere - maxHalfSamp):pkHere;
                      [~, rT] = min(sig(sw)); loHere = sw(rT);
                  end
                  if hiHere - pkHere > maxHalfSamp
                      sw = pkHere:(pkHere + maxHalfSamp);
                      [~, rT] = min(sig(sw)); hiHere = sw(rT);
                  end
                  lenHere = hiHere - loHere + 1;
                  if lenHere < cfg.segMinLenSamples || lenHere > cfg.segMaxLenSamples, continue; end
                  if pkHere <= minLimIdx, continue; end

                  runs{k} = (loHere:hiHere).';
                  nRescued = nRescued + 1;
                end
              end
            end
% how different is this from plain peak order?
            nTake = min([nK, nCc]);
            for kk2 = 1:nTake
                if isempty(runs{kk2}) || runs{kk2}(1) ~= kS(kk2)
                    nDisagree = nDisagree + 1; break
                end
            end
        case 'contactcycle'
            nC = min(numel(contactRel), cfg.numLicksToAnalyze);

% Typical inter-lick interval for THIS trial. Used only as the
% fallback gap at the two ends of the train (lick 1 has no previous
% contact, the last lick has no next one); the trough-search backstop
% uses the local gap.
            if numel(contactRel) >= 2
                ILI = median(diff(contactRel));
            else
                ILI = cfg.contactMatchTol_s;
            end
            if ~isfinite(ILI) || ILI <= 0, ILI = cfg.contactMatchTol_s; end

% ---- STEP 1: one peak per contact, as a windowed MAXIMUM ----
            pkIdx = nan(nC,1);
            for k = 1:nC
                tk = contactRel(k);
                if tk < timeSub(1) || tk > timeSub(end), continue; end

% half-distance to each neighboring contact bounds the window;
% at the ends fall back to half the median interval
                if k > 1
                    backHalf = (tk - contactRel(k-1)) / 2;
                else
                    backHalf = ILI / 2;
                end
                if k < numel(contactRel)
                    fwdHalf = (contactRel(k+1) - tk) / 2;
                else
                    fwdHalf = ILI / 2;
                end
                w0t = tk - 2*cfg.peakSearchFrac*backHalf;
                w1t = tk + 2*cfg.peakSearchFrac*fwdHalf;

                w0 = round(interp1(timeSub, (1:nSamp).', w0t, 'linear','extrap'));
                w1 = round(interp1(timeSub, (1:nSamp).', w1t, 'linear','extrap'));
                w0 = max(1, min(nSamp, w0));
                w1 = max(1, min(nSamp, w1));
                if w1 <= w0, continue; end

                [~, relMax] = max(sig(w0:w1));
                pkIdx(k) = w0 + relMax - 1;
            end

% ---- STEP 2 + 3: troughs by PROMINENCE WALK ----
% Step outward from the peak tracking the running minimum, and stop as soon
% as the trace has climbed back above that minimum by more than a threshold,
% i.e. it is out of the trough and into the next rise. The running minimum at
% that moment is the trough.
% The threshold is cfg.troughRiseFrac of the excursion depth so far,
% with cfg.troughRiseFloor as an absolute noise floor, so:
%   - a plateau wiggle near the peak cannot stop the walk (depth is
%     small there, so the floor governs and the wiggle is under it),
%   - a real climb out of the trough always does (depth is large by
%     then, so a genuine rise clears the fraction),
%   - and nothing depends on a window being the right width (edge licks and
%     long inter-lick gaps are handled the same way).
            segStart = nan(nC,1); segEnd = nan(nC,1);
            for k = 1:nC
                pk = pkIdx(k);
                if ~isfinite(pk), continue; end

% hard bounds: the neighboring CONTACT peaks, so a cycle can
% never cross into an adjacent lick's excursion
                prevPk = 1;
                kk = find(isfinite(pkIdx(1:k-1)), 1, 'last');
                if ~isempty(kk), prevPk = pkIdx(kk); end
                nextPk = nSamp;
                kk = find(isfinite(pkIdx(k+1:end)), 1, 'first');
                if ~isempty(kk), nextPk = pkIdx(k+kk); end

% Backstop from the LOCAL gap: a lick with a wide gap on one side gets a
% correspondingly wide allowance on that side. It is only a backstop --
% the walk normally stops well inside it.
                if k > 1
                    gapBack = contactRel(k) - contactRel(k-1);
                else
                    gapBack = ILI;
                end
                if k < numel(contactRel)
                    gapFwd = contactRel(k+1) - contactRel(k);
                else
                    gapFwd = ILI;
                end
                capBack = max(2, round(cfg.segMaxHalfCycleFrac * gapBack / dtSub));
                capFwd  = max(2, round(cfg.segMaxHalfCycleFrac * gapFwd  / dtSub));

                lo0 = max([prevPk, pk - capBack, 1]);
                hi0 = min([nextPk, pk + capFwd,  nSamp]);
                if lo0 >= pk || hi0 <= pk, continue; end

% MINIMUM DEPTH BEFORE THE WALK IS ALLOWED TO STOP.
% When the jaw reaches maximum opening it can hang high and jitter for tens of ms
% before closing. A jitter is a rise above the running minimum
% and could stop the walk up on the plateau, calling a high point the trough.
% Rises alone cannot distinguish the two, because the shapes are
% the same; what separates them is HOW FAR DOWN the trace has
% come. So: look ahead over the whole permitted interval, see how
% deep this side can possibly go, and refuse to let the walk stop
% until it has covered cfg.troughMinDepthFrac of that. A plateau
% wiggle sits at a small fraction of the full depth and is
% ignored; a rise near the bottom still stops the walk
% before it can run into the next cycle.
                needBack = cfg.troughMinDepthFrac * (sig(pk) - min(sig(lo0:pk)));
                needFwd  = cfg.troughMinDepthFrac * (sig(pk) - min(sig(pk:hi0)));

% --- trough BEFORE: walk backward from the peak ---
                runMin = sig(pk); runIdx = pk;
                for i = (pk-1):-1:lo0
                    if sig(i) < runMin
                        runMin = sig(i); runIdx = i;
                    else
                        depth = sig(pk) - runMin;
                        if depth >= needBack && ...
                           (sig(i) - runMin) >= max(cfg.troughRiseFrac*depth, cfg.troughRiseFloor)
                            break   % deep enough AND climbing out; runIdx is it
                        end
                    end
                end
                lo = runIdx;
% did the walk run all the way to its bound without turning?
% That means the trough is off-screen or past the neighboring
% peak, not that we found one.
                truncBack = (lo <= lo0);
% contiguous flat basin only: if the bottom is a flat stretch
% rather than a sharp V, start at the end of THAT stretch, the
% moment opening begins. Only the contiguous stretch counts, so a
% shoulder further up the flank cannot win.
                while lo < pk && (sig(lo+1) - sig(runIdx)) <= cfg.troughFlatTol
                    lo = lo + 1;
                end

% --- trough AFTER: walk forward from the peak ---
                runMin = sig(pk); runIdx = pk;
                for i = (pk+1):hi0
                    if sig(i) < runMin
                        runMin = sig(i); runIdx = i;
                    else
                        depth = sig(pk) - runMin;
                        if depth >= needFwd && ...
                           (sig(i) - runMin) >= max(cfg.troughRiseFrac*depth, cfg.troughRiseFloor)
                            break
                        end
                    end
                end
                hi = runIdx;
                truncFwd = (hi >= hi0);
                while hi > pk && (sig(hi-1) - sig(runIdx)) <= cfg.troughFlatTol
                    hi = hi - 1;
                end

                if hi <= lo, continue; end

% is this boundary sitting on the edge of the analysis window
% rather than on a real trough?
                edgeTrunc = (truncBack && lo0 == 1) || (truncFwd && hi0 == nSamp);
                if edgeTrunc && ~cfg.allowEdgeTruncated
                    continue
                end

% --- amplitude gate, side-aware ---
% For lick 1 the opening trough is often off the left edge of the analysis
% window (the jaw is already open at the cut), so a boundary that hit the
% window edge without turning is not evidence about amplitude and is excluded
% from the test. If both sides are truncated, the shallower side stands.
                ampBack = sig(pk) - sig(lo);
                ampFwd  = sig(pk) - sig(hi);
                if truncBack && ~truncFwd
                    amp = ampFwd;
                elseif truncFwd && ~truncBack
                    amp = ampBack;
                else
                    amp = min(ampBack, ampFwd);
                end
                if amp < cfg.minCycleAmp
                    nAmpFail = nAmpFail + 1;
                    continue
                end

                segStart(k) = lo;
                segEnd(k)   = hi;
                nEdgeTrunc  = nEdgeTrunc + double(edgeTrunc);
            end

% ---- optional fixed half-window shape, same anchors ----
            if strcmpi(cfg.segMode, 'fixedhalfwin')
                keep = isfinite(pkIdx) & isfinite(segStart);
                segStart(keep) = pkIdx(keep) - cfg.segHalfWin;
                segEnd(keep)   = pkIdx(keep) + cfg.segHalfWin;
            end

% ---- validity, then straight assignment: segment k IS contact k ----
            for k = 1:nC
                lo = segStart(k); hi = segEnd(k);
                if ~isfinite(lo) || ~isfinite(hi), continue; end
                if lo < 1 || hi > nSamp, continue; end
                if pkIdx(k) <= minLimIdx, continue; end
                len = hi - lo + 1;
                if len < cfg.segMinLenSamples || len > cfg.segMaxLenSamples, continue; end
                runs{k} = (lo:hi).';
            end

% how different is this from plain peak-order? a scale for how much
% the anchor is doing on this dataset
            nTake = min([nFP, nC, cfg.numLicksToAnalyze]);
            for kk2 = 1:nTake
                if ~isfinite(pkIdx(kk2)) || abs(pkIdx(kk2) - fpLocs(kk2)) > cfg.peakMinDist
                    nDisagree = nDisagree + 1; break
                end
            end

        case 'peakdetect'
% findpeaks peaks, each bounded by the NEAREST trough on either side, then
% matched to contacts. Nearest-trough boundaries can land on plateau
% wiggles; this path is kept for comparison.
            peakLocs = fpLocs;
            if isempty(peakLocs)
                segDiag(tr,:) = [numel(contactRel), 0, 0, 0, 0, 0];
                continue
            end

            [troughVals, troughLocs] = findpeaks(-sig, 'MinPeakProminence',cfg.troughProminence, ...
                'MinPeakDistance',cfg.troughMinDist);
            troughVals = -troughVals;
            troughLocs = troughLocs(:);
            if cfg.requireTroughBelowMax
                troughLocs = troughLocs(troughVals <= cfg.troughMaxValue);
            end

            nP = numel(peakLocs);
            segStart = nan(nP,1); segEnd = nan(nP,1);
            for i = 1:nP
                pk = peakLocs(i);
                prevT = troughLocs(troughLocs < pk);
                if ~isempty(prevT)
                    lo = prevT(end);
                else
                    loBound = 1; if i > 1, loBound = peakLocs(i-1); end
                    segw = loBound:pk; [~, rel] = min(sig(segw)); lo = segw(rel);
                end
                nextT = troughLocs(troughLocs > pk);
                if ~isempty(nextT)
                    hi = nextT(1);
                else
                    hiBound = nSamp; if i < nP, hiBound = peakLocs(i+1); end
                    segw = pk:hiBound; [~, rel] = min(sig(segw)); hi = segw(rel);
                end
                if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, continue; end
                segStart(i) = lo; segEnd(i) = hi;
            end

            if strcmpi(cfg.segMode, 'fixedhalfwin')
                segStart = peakLocs - cfg.segHalfWin;
                segEnd   = peakLocs + cfg.segHalfWin;
            end

            okPeak = isfinite(segStart) & isfinite(segEnd) & ...
                     segStart >= 1 & segEnd <= nSamp & ...
                     peakLocs > minLimIdx & ...
                     (segEnd - segStart + 1) >= cfg.segMinLenSamples & ...
                     (segEnd - segStart + 1) <= cfg.segMaxLenSamples;
            keepPeaks = peakLocs(okPeak);
            keepStart = segStart(okPeak);
            keepEnd   = segEnd(okPeak);
            if isempty(keepPeaks)
                segDiag(tr,:) = [numel(contactRel), 0, nFP, 0, 0, 0];
                continue
            end

            switch cfg.jawLickIndexRule
                case 'peakOrder'
                    nTake = min(numel(keepPeaks), cfg.numLicksToAnalyze);
                    for kk = 1:nTake
                        runs{kk} = (keepStart(kk):keepEnd(kk)).';
                    end
                case 'contactMatched'
% global, distance-ordered assignment
                    peakTimes = timeSub(keepPeaks);
                    nC = min(numel(contactRel), cfg.numLicksToAnalyze);
                    nK = numel(keepPeaks);
                    D = abs(repmat(contactRel(1:nC), 1, nK) - repmat(peakTimes(:).', nC, 1));
                    D(D > cfg.contactMatchTol_s) = Inf;
                    while true
                        [minVal, linIdx] = min(D(:));
                        if ~isfinite(minVal), break; end
                        [ci, pj] = ind2sub(size(D), linIdx);
                        runs{ci} = (keepStart(pj):keepEnd(pj)).';
                        D(ci, :) = Inf;
                        D(:, pj) = Inf;
                    end
                    nTake = min(nK, cfg.numLicksToAnalyze);
                    for kk = 1:nTake
                        orderSeg = (keepStart(kk):keepEnd(kk)).';
                        if isempty(runs{kk}) || ~isequal(runs{kk}, orderSeg)
                            nDisagree = nDisagree + 1; break
                        end
                    end
                otherwise
                    error('cfg.jawLickIndexRule must be ''contactMatched'' or ''peakOrder''.');
            end

        otherwise
            error('cfg.segAnchor must be ''persistCycle'', ''contactCycle'' or ''peakDetect''.');
        end

        lickRuns{tr} = runs;
        segDiag(tr,:) = [numel(contactRel), sum(~cellfun(@isempty, runs)), nFP, nAmpFail, nEdgeTrunc, nRescued];
    end
end

function [gcRel, licksRel] = localJawWarpLicks(trialIdx, goCueAll, lickLAll, nLicks)
% Real lick-contact times for this trial (not jaw-oscillation peaks), from the
% goCue and lickL arrays of the slim session struct.
% goCue/lickL are on each trial's own ABSOLUTE clock, while
% traceTime is resampled per trial so THAT trial's first post-go-cue contact
% sits at 0 (params.alignEvent = 'firstLick'). So both are converted to be
% relative to that contact before use. By construction licksRel(1) == 0.
    gcRel = NaN;
    licksRel = [];
    gcAbs = goCueAll(trialIdx);
    licksAbs = lickLAll{trialIdx, 1};
    if isempty(licksAbs) || all(isnan(licksAbs)), return; end
    licksPostGC = sort(licksAbs(licksAbs > gcAbs));
    if isempty(licksPostGC), return; end
    firstLickAbs = licksPostGC(1);
    licksRel = licksPostGC - firstLickAbs;
    gcRel    = gcAbs - firstLickAbs;
    if numel(licksRel) > nLicks, licksRel = licksRel(1:nLicks); end
end

function medStart = localJawWarpMedian(contactCell, nLicks, delta)
% Warp template from real contact times: per lick index, the median across
% trials that have that many licks, then replaced by uniform spacing delta
% from lick 1. licksRel(1) is always exactly 0, so the template
% always starts at 0.
    medStart = nan(nLicks,1);
    for i = 1:nLicks
        hasI = cellfun(@(x) numel(x) >= i, contactCell);
        if any(hasI)
            ith = cellfun(@(x) x(i), contactCell(hasI));
            medStart(i) = median(ith);
        end
    end
    if isfinite(medStart(1))
        t0 = medStart(1);
        medStart = (t0 : delta : t0 + delta*(nLicks-1))';
    end
end

function pfit = localJawWarpFits(contactCell, goCueList, medStart)
% Warp fits on real contact times: one affine map per
% inter-lick interval (goCue->lick1, lick1->lick2, ..., last lick->end), taking
% this trial's boundaries onto the across-trial median template. The go-cue
% anchor is nanmedian(goCue) across trials.
% Interior boundaries use the exact shared anchor: the
% resampler below inverse-maps OUTPUT grid points per interval rather than
% masking input samples, so there is no overlap to guard against, and sharing
% the endpoint keeps the warp continuous instead of leaving a one-sample gap
% (and a visible tick) at every lick.
    nTrials = numel(contactCell);
    nLicks  = numel(medStart);
    pfit = cell(nTrials, nLicks);
    if any(isnan(medStart)), return; end
    mls = medStart;
    gcMed = nanmedian(goCueList);
    for trix = 1:nTrials
        ls = contactCell{trix};
        gc = goCueList(trix);
        for lix = 1:numel(ls)
            if lix == 1 && lix == numel(ls)
                x = [gc, ls(1)+median(diff(mls))];
                y = [gcMed, mls(1)+median(diff(mls))];
            elseif lix == 1
                x = [gc, ls(2)];
                y = [gcMed, mls(2)];
            elseif lix == numel(ls)
                x = [ls(lix), ls(lix)+median(diff(mls))];
                y = [mls(lix), mls(lix)+median(diff(mls))];
            else
                x = [ls(lix), ls(lix+1)];
                y = [mls(lix), mls(lix+1)];
            end
            if numel(x)==2 && x(2) > x(1) && all(isfinite(x)) && all(isfinite(y))
                pfit{trix,lix} = polyfit(x,y,1);
            end
        end
    end
end

function yOut = localJawWarpResample(timeVec, yTrial, goCue, contacts, medStart, pfitRow, gridT, edgeMode)
% Resamples yTrial (on the always-monotonic timeVec) onto the template grid via
% the INVERSE of the piecewise-affine warp: per interval, ask "which input
% sample lands at this OUTPUT grid point" rather than "where does this input
% sample land". The forward direction would need one globally monotonic combined
% time vector, which real lick timing against a fixed template spacing can
% violate. Grid points outside every mapped interval are handled by edgeMode:
    ls = contacts;
    mls = medStart;
    nL = numel(ls);
%   'identity' keep the UNWARPED value. Samples before the first interval are then
%              on the raw clock and samples after it on the warped clock, which
%              leaves a step discontinuity at the median go cue in both the
%              actual and the predicted trace.
%   'extend'   the FIRST interval's affine map is extrapolated
%              backwards to cover everything before it, and the LAST interval's
%              forwards. Same map on both sides of the seam, so the trace is
%              continuous by construction and the pre-cue baseline is still real
%              data, just warped by the nearest available map.
%   'nan'      leave the unmapped edges empty.
    if nargin < 8 || isempty(edgeMode), edgeMode = 'identity'; end
    if strcmpi(edgeMode, 'identity')
        yOut = interp1(timeVec, yTrial, gridT, 'linear', NaN);
    else
        yOut = nan(size(gridT));
    end
    if nL == 0, return; end
    for lix = 1:nL
        pCur = pfitRow{lix};
        if isempty(pCur), continue; end
        if lix == 1 && lix == nL
            targetLo = polyval(pCur, goCue);
            targetHi = polyval(pCur, ls(1)+median(diff(mls)));
        elseif lix == 1
            targetLo = polyval(pCur, goCue);
            targetHi = polyval(pCur, ls(2));
        elseif lix == nL
            targetLo = polyval(pCur, ls(lix));
            targetHi = polyval(pCur, ls(lix)+median(diff(mls)));
        else
            targetLo = polyval(pCur, ls(lix));
            targetHi = polyval(pCur, ls(lix+1));
        end
        idx = find(gridT >= targetLo & gridT <= targetHi);
        if isempty(idx), continue; end
        tOrig = (gridT(idx) - pCur(2)) / pCur(1);   % inverse of the affine map
        yOut(idx) = interp1(timeVec, yTrial, tOrig, 'linear', NaN);
    end

% ---- the edges, once the interior is done ----
    if strcmpi(edgeMode, 'extend')
        fi = find(~cellfun(@isempty, pfitRow(1:nL)), 1, 'first');
        li = find(~cellfun(@isempty, pfitRow(1:nL)), 1, 'last');
        if ~isempty(fi)
            pCur = pfitRow{fi};
            if fi == 1, tLo = polyval(pCur, goCue); else, tLo = polyval(pCur, ls(fi)); end
            idx = find(gridT < tLo);
            if ~isempty(idx)
                yOut(idx) = interp1(timeVec, yTrial, (gridT(idx) - pCur(2))/pCur(1), 'linear', NaN);
            end
        end
        if ~isempty(li)
            pCur = pfitRow{li};
            if li == nL, tHi = polyval(pCur, ls(li)+median(diff(mls)));
            else,        tHi = polyval(pCur, ls(li+1)); end
            idx = find(gridT > tHi);
            if ~isempty(idx)
                yOut(idx) = interp1(timeVec, yTrial, (gridT(idx) - pCur(2))/pCur(1), 'linear', NaN);
            end
        end
    end
end

function h = localBandPlot(axHandle, xVec, Ymat, col, lw)
% Mean +/- 1.96 SEM band with the mean line.
    nF = sum(isfinite(Ymat), 2);
    mu = nanmean(Ymat, 2);
    se = nanstd(Ymat, 0, 2) ./ sqrt(max(nF,1));
    se(nF < 2) = NaN;
    lo = mu - 1.96*se; hi = mu + 1.96*se;
    good = isfinite(lo) & isfinite(hi);
    if any(good)
        xv = xVec(good);
        patch(axHandle, [xv(:); flipud(xv(:))], [lo(good); flipud(hi(good))], ...
            col, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off');
    end
    h = plot(axHandle, xVec, mu, '-', 'Color', col, 'LineWidth', lw);
end

function [kSel, sse] = selectLambdaCV5(X, y, rowTrial, cfg)
% Lambda by k-fold cross-validation WITHIN THE TRAINING SET: fit on k-1 folds,
% predict the held-out fold, and pool squared errors across folds.
% Squared errors are POOLED across folds -- summed, not averaged per fold --
% so folds that happen to contain longer trials carry proportionally more
% weight, which is what "pooling squared errors" means and what makes the
% criterion a plain held-out SSE.
% FOLDS ARE OVER WHOLE TRIALS, NOT ROWS. Each row of X is one time sample, and
% consecutive samples from the same trial are strongly autocorrelated. Splitting
% rows at random would leave a held-out sample's immediate neighbors sitting in
% the training set, so the held-out error would measure interpolation rather
% than generalisation and would favour a lambda that is too small. rowTrial comes
% from laggedDesign and says which trial each row belongs to.
% The caller takes the coefficients from the FULL training-set path at kSel, so
% the fold fits are used only to score lambda, never to produce a weight.
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

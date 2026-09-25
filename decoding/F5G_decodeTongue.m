%% F5G_decodeTongue.m
%  Ridge decoder of tongue length from population spike rates, Simple Reward Task, first five days of training.
%  Spike rates are z-scored and lagged, then fit per session with a ridge
%  penalty chosen by cross-validation over held-out training trials. The
%  decoding index is 1 - RMSE(decoded) / RMSE(zero baseline), computed for
%  each lick contact and averaged across sessions.
%  Produces Fig. 5G.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'firstLick'
%    params.dt          1/300
%    params.smooth      10
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2.5 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


%% RUN SETTINGS (edit these two)
% HOW MANY TIMES TO RUN THE WHOLE ANALYSIS.
%   1 = one pass, the normal case.
%       sets -- nothing is averaged or pooled, so the statistics in each pass
%       are exactly the statistics of a single run. Compare the passes by eye.
% WHAT CHANGES BETWEEN PASSES: the train/test trial split AND the
% cross-validation fold assignment are both re-drawn, so a pass answers "how
% much of this result depends on which trials happened to land where". Nothing
% Each pass gets its own set of figure windows, titled with the pass number.
% A pass is a full run, so N passes take N times as long.
% various points in these files -- the count and the wrap point differed from
% script to script, so the scripts were not doing the same number of passes as
% each other.)
RUN.nRepeats = 1;

% CROSS-VALIDATION FOLDS for choosing lambda, inside the training set.
% Methods currently say five; this is the single place to change it. Fewer
% folds train on less data per fold, which biases held-out error up and tends
% to pick a slightly larger lambda.
RUN.cvFolds  = 5;
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
%% SETTINGS (shared, identical in all five)
rewardLickB = 4;   % second reward-lick condition string only; touches no timing

%% TIMING: the three windows, all in seconds
% Everything is relative to t = 0, which IS the first lickport contact after the
% go cue (params.alignEvent = 'firstLick').

% FITTING WINDOW -- genuinely per trial, and variable length.
%            produced contact 1) + postTongueRetract_s
% Contact 1 is the moment the tongue touches the port, and that happens PART WAY
% THROUGH the protrusion -- the tongue stays out, often for tens of milliseconds,
% and only then retracts. The window is closed by the RETRACTION, not by the
% contact: take the run of visible (non-zero) tongue length that CONTAINS t = 0,
% back in the mouth -- and pad from there.
% the parameter is renamed so the code says what it does.
cfg.preGoCue_s          = 0.20;   % window opens this long BEFORE the go cue
cfg.postTongueRetract_s = 0.005;   % window closes this long AFTER the tongue is
% last visible in that first protrusion
cfg.minWindowSamples    = 10;   % trials with a shorter window are dropped

%% WHICH FITTING WINDOW (two options)
% 'contact1'  THE DEFAULT, and exactly what the block above describes:
%             go cue - preGoCue_s  ->  the tongue retracts after contact 1,
%             + postTongueRetract_s. Nothing about it has changed.
% 'bout'      THE SECOND OPTION. From winBoutStartPad_s AFTER CONTACT 2, to
%             winBoutEndPad_s AFTER THE CONTACT THAT ENDS THE BOUT -- the first
%             contact that is not followed by another one within winBoutGap_s.
%             Contact 1 and its protrusion are OUTSIDE this window entirely.
% 'lick2ToBoutEnd'  END OF LICK 2 -> END OF BOUT. Same end as 'bout', but the
%             start is the RETRACTION of lick 2's protrusion (its last visible
%             frame) + winBoutStartPad_s, not contact 2's time. This is the one
%             that reads "from the end of lick 2 to the end of the bout".
% 'contact1' the model only ever sees the FIRST protrusion, so contacts 2-8 are
% extrapolation: the decline across the bout could be the model leaving the
% explanations make opposite predictions here. Under 'bout' the training window
% covers contacts 2 onwards and contact 1 is the extrapolation, so:
%   - if the decline across contacts SURVIVES with the window moved to the late
%     contacts, it is not a domain-shift artefact;
% postTongueRetract_s; 'bout' reads the three winBout* values and ignores the
% other two. Nothing else in the script changes between modes -- the same
% z-scoring, the same lags, the same lambda rule, the same decoding index over
% contacts 1..nLicksAnalyze.
cfg.fitWindowMode     = 'contact1';   % EARLY model = Fig. 5G (was 'contact1'). Options: 'contact1' | 'bout' | 'lick2ToBoutEnd'
cfg.winBoutStartPad_s = 0.02;   % window opens this long AFTER contact 2
cfg.winBoutEndPad_s   = 0.20;   % window closes this long AFTER the bout's last contact
cfg.winBoutGap_s      = 0.75;   % a gap longer than this is what ENDS the bout

% A one- or two-frame DLC dropout in the middle of a protrusion would otherwise
% split it into two runs, and the window would end at the DROPOUT instead of at
% the retraction -- ending the fitting window early and, on the diagnostic panel,
% making one lick look like two. Runs separated by this many zero samples or
% fewer are merged into a single protrusion before anything else happens.
% 2 samples = 10 ms at dt = 5 ms. Set to 0 to disable the merge.
cfg.mergeGapSamples = 2;

% NEURAL Z-SCORING WINDOW
cfg.zscoreWin_s = [-2.0 2.0];   % Methods: -2 s to +2 s around the first port contact
cfg.minBaselineStd = 1e-6;   % units below this are DROPPED, not divided by

% NEURAL LAG CONTEXT (lag set is -preBins : +postBins, inclusive)
cfg.lagPre_s  = 0.06;
cfg.lagPost_s = 0.04;

% ANALYSIS WINDOW for the per-contact decoding index
cfg.analysisWin_s = [-0.09 3.0];

%% TARGET: tongue length
% means length zero. A real observation, not missing data, so it is filled with
% zero and never masked or excluded anywhere.
cfg.tongueZeroMode = 'p_low';   % 'p_low' | 'anatomical'
cfg.normPctLow  = 2;   % Methods: 1st percentile. Unified across all
cfg.normPctHigh = 99;   % Methods: 99th percentile

%% RIDGE
cfg.ridgeGrid       = logspace(-2, 5, 50);
cfg.standardizeCols = true;   % penalty lambda*sum(beta_j^2 * sigma_j^2);
% required by cfg.ridgeImpl = 'matlab'

cfg.ridgeImpl = 'eig';   % the inline closed form -- NO call to MATLAB's
% verified head to head (see ridgePath); ~30x faster.
% 'eig' is the algebraically identical closed
% form and is much faster when p is large --
% MATLAB's ridge does one augmented
% least-squares solve PER lambda, so 41 grid
% points with thousands of predictors is 41
% QR factorisations per fit, per CV fold.
cfg.ridgeCheckEquivalence = false;   % fit the first session both ways and print
% the max relative difference (expect ~1e-11)

% LAMBDA -- k-FOLD CROSS-VALIDATION INSIDE THE TRAINING SET, k = cfg.cvFolds,
% which is currently 3. NOTE: the Methods say FIVE-fold; see the divergence
% note in the header block at the top of this file. Folds are over whole TRIALS, never rows (see
% What follows describes that older rule and is kept for reference. There are two
% sets of trials in these scripts and only two: TRAIN and TEST. Nothing is held
% out of the training trials.
% Lambda is chosen by generalized cross-validation on the TRAINING FIT ITSELF:
% edf is the effective number of parameters the penalty leaves standing (+1 for
% the unpenalized intercept); it falls from p towards 0 as lambda grows, so the
% denominator penalizes a fit that is only good because it is flexible. This is
% the rotation-invariant limit of leave-one-out and it costs nothing extra -- the
% eigenvalues are already computed for the lambda path. No row of data is ever
% set aside, and the TEST trials are never touched during fitting.
cfg.lambdaRule  = 'cv5';   % 'cv5' is a SWITCH LABEL ONLY -- the fold count
% comes from cfg.cvFolds below (now 3), not from
% this string. (k-fold CV in the training
% set); 'gcv' = the previous behaviour
cfg.cvFolds     = RUN.cvFolds;   % set at the TOP of this file (RUN block).
% this, so every printout and axis label
% follows automatically.
cfg.cvSeed      = [];   % base seed; each pass offsets it by repIdx
RUN.cvSeed0     = cfg.cvSeed;   % remembered so the offset is from a fixed base            % [] = folds in trial order, no shuffle
cfg.lambdaFixed = [];   % [] = GCV. Set a number to pin lambda instead.

%% CONTACT DETECTION + DECODING INDEX
% ONLY the zero baseline:  decodingIndex(L) = 1 - RMSE_decoded(L)/RMSE_zero(L)
% RMSE_zero is the error a constant-zero predictor makes over contact L's
% samples. Being a CONSTANT it shrinks whenever the excursion shrinks, so a
% purely behavioural amplitude decline lowers the index on its own -- which is
% why the per-contact actual amplitude is plotted beside it everywhere.
cfg.nLicksAnalyze     = 8;
cfg.minContactSamples = 6;   % 6 samples at 1/300 = 20 ms, as at 1/200 with 4 (and as tongue_r1)

%% BOUT REQUIREMENT (one toggle, on or off)
% Keep a trial only if the animal produced a real BOUT on it: at least
% boutMinLicks port contacts IN A ROW, consecutive ones no more than
% boutMaxILI_s apart, all of them inside boutWin_s of the go cue.
% cfg.requireBout is the switch. It is applied to the WHOLE trial pool, so it
% removes training trials as well as test trials -- that is what "exclude trials"
% means, and it also keeps the model trained on the same kind of trial it is
% scored on. Flip it off and re-run to see what it is doing to any number.
% down from 8 to 6. Six contacts is five intervals, and five intervals at the
% 250 ms ceiling span 1.25 s -- inside the 1.5 s window. So it is now
% boutMaxILI_s that does the work, and boutWin_s only rejects a bout that starts
% many trials are dropped, loosen boutMaxILI_s first -- the fit loop prints the
% count for every session.
% NOT APPLIED IN THE VTA SCRIPT (spec.applyBoutFilter = false there). Its split
% rule already divides trials by lick count -- fewer than 4 contacts trains,
% 4 or more tests -- so a 6-lick floor would delete the training set outright.
cfg.requireBout   = true;   % <-- the toggle
cfg.boutMinLicks  = 4;
cfg.boutWin_s     = 1.50;
cfg.boutMaxILI_s  = 0.25;

%% STATISTICS
cfg.earlyLicks = 1:1;
cfg.lateLicks  = 2:8;
cfg.alpha      = 0.05;
% per-animal, trial-level Day 1 vs Day 5 test
cfg.perAnimalTest     = true;   % false = skip it
cfg.perAnimalLicks    = 2:8;   % drop = mean(index over these licks) - index at lick 1
cfg.perAnimalMinLicks = 3;   % a trial needs lick 1 and >= this many of perAnimalLicks

%% SESSION CACHE
% loadSessionData + getKinematics + loadMotionEnergy is where essentially all
% the wall-clock time goes; the ridge fit itself is under a second. None of them
% depend on anything in cfg, so with the cache on, every re-run after a window /
% z-score / lag / penalty change skips loading entirely and goes straight to the
% fit. The first run fills the cache and is no slower than before.
% The cache key is the full set of params fields that can change what gets
% ...). Change any of those and the entry is recomputed automatically -- a stale
% Budget roughly (nTime x nUnits x nTrials x 4 bytes) per session: trialdat is
% stored as single, so a 1500 x 300 x 400 session is about 700 MB. Point
% cacheDir at a scratch disk, and delete the folder any time to force a reload.
cfg.useCache = true;
cfg.cacheDir = fullfile(tempdir, 'tlDecodeCache_clean');

%% HOUSEKEEPING
cfg.rngSeed          = [];
cfg.sessionsToRun    = [];   % [] = all

% DIAGNOSTIC FIGURES. [] means "choose automatically", by spec.diagRule:
%   'random' (R1, R14, R16, VTA)  -- cfg.nDiagSessions sessions drawn at random
%                                    from the ones that actually fit. The draw
%                                    uses cfg.rngSeed, so it is reproducible;
%                                    change the seed to see a different five.
%   'groups' (Learning)           -- every session in spec.diagGroups, which is
%                                    day 1 and day 5, so 8 figures.
% Put explicit session numbers in cfg.diagSessions to override either.
% Every diagnostic session becomes a TAB, not a separate window: one figure with
% the six-panel fit diagnostics (tab per session) and a second figure with the
% warped lick train (tab per session). Click along the tab strip to go through
% them.
cfg.diagSessions     = [];
cfg.nDiagSessions    = 5;

% The WARPED figure can be decoupled from the six-panel one. It is the panel that
%   true   every fitted session gets a warped tab, whatever the diagnostic rule picked
%   false  the warped figure follows the same session list as the six-panel one
cfg.warpAllSessions  = false;
cfg.warnPredictorCount = 8000;

% first slot at 10, slots every 50 (so a 20-sample NaN gap between licks), up to
% 10 licks extracted and the first 8 drawn.
cfg.warpedSegmentLength = 30;   % samples each excursion is resampled onto
cfg.warpedSlotFirst     = 10;   % where lick 1's slot starts
cfg.warpedSlotSpacing   = 50;   % slot pitch; spacing - segment = the NaN gap
cfg.maxLicksToExtract   = 10;
cfg.maxLickShow         = 8;

% Two cosmetic choices that were NOT in the original and default to its behaviour.
%   warpedInterp     'pchip' is what the original used. 'linear' is flatter and
%                    cannot overshoot at a segment's ends.
%   warpedTickCentre false puts the "contact k" label at the LEFT EDGE of the
%                    slot, as the original did. true centres it under the
%                    excursion, which reads better when the gaps are narrow.
cfg.warpedInterp     = 'pchip';   % 'pchip' | 'linear'
cfg.warpedTickCentre = true;

% With TWO test sets (R1 vs R4, R1 vs R6), draw each on its OWN stacked axes --
% actual black, predicted in that test set's colour, labelled in-plot -- instead
% of pooling them into one pair of traces. Ignored where there is one test set.
cfg.warpSplitTestSets = true;

% THE DENSER LAYOUT (what the first tabbed version drew, before the port): a 25-
% sample segment with a 1-sample gap, linear interpolation, centred labels. Paste
% these five lines over the ones above to get it back.
% old repeat loop: for kkkkk = 1:5  -> now RUN.nRepeats at the top
%% PARAMS (field names read by the loaders)
params.alignEvent = 'firstLick';
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 8;
params.lowFR  = 0.01;   % minimum mean firing rate, Hz
params.quality    = {'good'};   % good units only (findClusters trims blanks and ignores case)
params.tmin   = -2.5;
params.tmax   = 5;
params.dt     = 1/300;   % 300 Hz as every other decoder and the Methods
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

% Condition list is identical in every script; index 8 is the R1 pool and index
% 9 the second reward-lick pool (4 or 6 depending on the task).
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
%% STUDY: learning days
spec.name = 'Learning';
spec.sessionDates = { ...
    '2025-02-03','2025-02-04','2025-02-05','2025-02-06','2025-02-07', ...   % TD3l
    '2025-02-03','2025-02-04','2025-02-05','2025-02-06','2025-02-07', ...   % TD2l
    '2025-06-03','2025-06-04','2025-06-05','2025-06-06','2025-06-07', ...   % TD4l
    '2025-08-21','2025-08-22','2025-08-23','2025-08-24','2025-08-25' };   % TD5l
spec.sessionLoaders = { ...
    @loadTD3l_neur, @loadTD3l_neur,    @loadTD3l_neur, @loadTD3l_neur, @loadTD3l_neur, ...
    @loadTD2l_neur, @loadTD2l_neur219, @loadTD2l_neur, @loadTD2l_neur, @loadTD2l_neur, ...
    @loadTD4l_neur, @loadTD4l_neur,    @loadTD4l_neur, @loadTD4l_neur, @loadTD4l_neur, ...
    @loadTD5l_neur, @loadTD5l_neur,    @loadTD5l_neur, @loadTD5l_neur, @loadTD5l_neur };
% Groups are LEARNING DAYS, not regions: day d picks session d of each animal.
spec.groupMaps = { ...
    [1 0 0 0 0   1 0 0 0 0   2 0 0 0 0   2 0 0 0 0] , ...   % Day 1
    [0 2 0 0 0   0 2 0 0 0   0 1 0 0 0   0 1 0 0 0] , ...   % Day 2
    [0 0 1 0 0   0 0 1 0 0   0 0 1 0 0   0 0 2 0 0] , ...   % Day 3
    [0 0 0 1 0   0 0 0 2 0   0 0 0 1 0   0 0 0 1 0] , ...   % Day 4
    [0 0 0 0 1   0 0 0 0 1   0 0 0 0 1   0 0 0 0 1] };   % Day 5
spec.groupLabels = {'Day 1','Day 2','Day 3','Day 4','Day 5'};
spec.poolCondIdx    = 1;
% spec.trainTrialType is unused by this rule and kept only so the 'trialType'
% THE DECODING INDEX IS BACK ON THE OLD SPLIT. 'trialType' trains on a quota of
% the trialTypes == 3 trials and tests on every other condition-1 trial, exactly
% as it did originally -- every number in the decoding-index figures, the deltas
% and the statistics comes from that split and nothing below changes it.
% 'lateTrainEarlyTest' (train on the last trainFrac, test on the first, with
% spec.testNumTrials for a fixed count) is still implemented and one word away if
spec.splitRule      = 'trialType';   % 'trialType' | 'lateTrainEarlyTest' | 'typeStratified'
spec.trainTrialType = 3;
spec.trainFrac      = 0.70;   % Methods: ~70% for the Learning task
spec.testNumTrials  = [];   % only read by 'lateTrainEarlyTest'

% ---- THE WARPED FIGURE ONLY ----
% [] = the warped figure shows the TEST trials, like every other figure.
% A NUMBER (25, say) = it instead shows the FIRST N trials of the pool, in
% acquisition order. This is a DISPLAY OVERRIDE and nothing else: the model is
% already fitted, the decoding index, the deltas and every statistic are
% !! AND IT IS NOT A HELD-OUT NUMBER !! The first N trials of the pool will
% mostly be TRAINING trials under the 'trialType' split, so the predicted trace
% there is an IN-SAMPLE fit and will track the actual better than the honest
% held-out prediction does. The figure says so in its header, and the console
% prints how many of the N were trained on. Do not quote amplitudes off it.
spec.warpTestFirstN = [];   % <-- set 25 (or whatever) to override the warped figure
spec.testSets       = struct('name','all','condIdx',{[]});
spec.trialCaps           = { 'TD4l','2025-06-05', 123 };   % TDl4 -> TD4l typo fixed
spec.singleProbeSessions = cell(0,2);
spec.lickCountWin_s      = 1.25;

spec.warpShowLicks    = 8;   % contacts LABELLED on the warped figure

%% bout filter + diagnostics
spec.applyBoutFilter = true;   % cfg.requireBout is the on/off switch
% DIAGNOSTICS: every day-1 and day-5 session, so 8 figures -- the first and the
% last day of learning, which is where the comparison is made.
spec.diagRule        = 'groups';
spec.diagGroups      = [];

%% colours
% Groups are learning days and they DO share an axes here, so the group colours
% are the ones that matter: DAY 1 BLACK through to DAY 5 RED, with days 2-4
% interpolated so the ordering reads off the figure without the legend.
spec.tsColours  = [0.85 0.10 0.10];   % one test set
spec.grpColours = [0.00 0.00 0.00 ;   % Day 1  black
                   0.30 0.08 0.08 ;   % Day 2
                   0.50 0.09 0.09 ;   % Day 3
                   0.68 0.10 0.10 ;   % Day 4
                   0.85 0.10 0.10];   % Day 5  red
spec.predColour = [0.85 0.10 0.10];   % prediction on the diagnostics

%% figure behaviour
%% Day 1 vs Day 5: first and last day of learning. (The earlier 'day 6' was a
%% typo -- these maps define five days and day 5 is the last one.)
spec.overlayGroups = true ;   % false = never draw two groups on one axes
spec.showFigCI     = false;   % Figure 1: decoding index with error bars
spec.showFigDelta  = false;   % Figure 3: delta from contact 1
spec.deltaGroups   = [1 5];   % groups on the delta figures ([] = all)
spec.figRef       = 'Fig. 5G';   % the manuscript panel this file produces
spec.vsC1Tail     = 'both';   % 'both' | 'left' | 'right', from the legend
spec.groupTest    = 'ranksum';   % 'none' | 'ranksum' | 'crossTask'
spec.groupSummaryLicks = 3:8;   % Day 1 vs Day 5: mean change from contact 1 over contacts 3-8
% the animal-level signed-rank now test the same quantity.
spec.summaryFile  = fullfile(tempdir, 'decodeSummary_tongue_learning_clean.mat');
spec.crossTaskFile = fullfile(tempdir, 'decodeSummary_tongue_r1_clean.mat');
spec.compareGroups = [1 5];   % [a b] = per-contact paired test between two groups
%% MAIN LOOP (shared)
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
switch lower(cfg.lambdaRule)
    case 'cv5'
        logf(['  lambda: %d-fold cross-validation over TRAINING TRIALS (Methods), ' ...
                 'refit on all training rows at the winning value\n'], cfg.cvFolds);
    case 'gcv'
        logf('  lambda: GCV on the training fit -- no folds (pre-Methods behaviour)\n');
end
logf('  diagnostics: sessions %s\n', mat2str(diagList));
% Which sessions get which figure. The warped train can cover every session
if cfg.warpAllSessions
    warpList = runList(~cellfun(@isempty, probesOf(runList)));
else
    warpList = diagList;
end
logf('  warped figure: sessions %s\n', mat2str(warpList));

% Colours come from the study block so each task keeps its own scheme; the shared
% code never hard-codes one.
tsCols  = spec.tsColours;
predCol = spec.predColour;

%% REPEAT PASSES
for repIdx = 1:RUN.nRepeats
fprintf('\n\n%s\n', repmat('#',1,78));
logf('#  PASS %d of %d   (%s)\n', repIdx, RUN.nRepeats, spec.name);
fprintf('%s\n\n', repmat('#',1,78));

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
logf('[pass %d] split seed %d | cv fold seed %d | %d-fold\n\n', ...
    repIdx, 90210 + 1000*repIdx, cfg.cvSeed, cfg.cvFolds);
if RUN.nRepeats > 1
    passTag = sprintf('pass %d/%d - ', repIdx, RUN.nRepeats);
else
    passTag = '';
end

% The two diagnostic figures are created ONCE and get a tab per session.
diagFig = [];  diagTG = [];
warpFig = [];  warpTG = [];

res = struct([]);  nFitted = 0;

for sessionIdx = runList
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
tongueLen(~visMask) = 0;   % tongue not visible -> length 0. Data, not imputation.

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
    case 'lateTrainEarlyTest'
% TRAIN ON THE LAST spec.trainFrac OF THE SESSION, TEST ON THE FIRST.
% poolTrials are trial NUMBERS, so sorting them is acquisition order.
% THIS IS A CHRONOLOGICAL SPLIT, NOT A RANDOM ONE, and that changes what
% the held-out number means. Train and test now differ systematically in
% time-on-task -- satiety, electrode drift, arousal and within-session
% learning all separate the two sets -- so a low test score can be the
% session moving under the model rather than the model being poor. That
% is the point here: the learning script asks whether a decoder fitted on
% practised. Read it as that comparison, not as a generic held-out score,
% and do not put it beside the randomly-split numbers from the other
% tasks without saying so.
% No sampling anywhere, so this split is identical on every re-run
% whatever cfg.rngSeed is.
% HOW BIG THE TEST SET IS -- two ways to say it, and spec.testNumTrials
% wins whenever it is set:
%                             trains. A FIXED COUNT, identical on every
%                             per-contact n has to be comparable across
%                             sessions -- a fraction gives a 20-trial test
%                             set on a short session and a 60-trial one on
%                             a long session, and the CI widths then differ
%                             decoding.
%                             (1 - spec.trainFrac) of the pool.
% A session that cannot supply the requested count, or that would be left
% with too few training trials, is SKIPPED rather than quietly given a
% smaller test set -- silently shrinking it would put sessions with
% different n on the same curve, which is the thing a fixed count exists
% to prevent. The skip line below names the session and the shortfall.
        ordTrials = sort(poolTrials(:));
        nAll      = numel(ordTrials);
        if isempty(spec.testNumTrials)
            nTestQ = nAll - round(spec.trainFrac * nAll);
            howStr = sprintf('%.0f%% of pool', 100*(1-spec.trainFrac));
        else
            nTestQ = spec.testNumTrials;
            howStr = sprintf('fixed %d', spec.testNumTrials);
        end
        nTestQ = max(nTestQ, 0);
        if nTestQ >= nAll
            logf('  [split] chronological | %d pool | asked for %d test trials -- not enough trials, session skipped\n', ...
                nAll, nTestQ);
            trainTrials = [];  testTrials = [];
        else
            testTrials  = ordTrials(1 : nTestQ);   % the FIRST  n, in order
            trainTrials = ordTrials(nTestQ+1 : nAll);   % everything AFTER them
            logf('  [split] chronological (%s) | %d pool | test = first %d (trials %d-%d) | train = last %d (trials %d-%d)\n', ...
                howStr, nAll, numel(testTrials), min(testTrials), max(testTrials), ...
                numel(trainTrials), min(trainTrials), max(trainTrials));
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
% k-fold CV inside the training set (k = cfg.cvFolds; the Methods say
% five, this runs three), squared errors pooled
% across folds, then refit on ALL training rows at the winning lambda --
% which is what R already holds, so only the index changes.
        [kSel, cvSSE] = selectLambdaCV5(Xtrain(okTr,:), yTrain(okTr), rowTrain(okTr), cfg);
        lamHow = sprintf('%d-fold trial-wise CV', cfg.cvFolds);
% the curve the diagnostic figure draws under panel 1
        lamCurve = cvSSE;
        lamCurveName = sprintf('%d-fold CV SSE', cfg.cvFolds);
        lamShort = sprintf('%d-fold CV', cfg.cvFolds);
    case 'gcv'
        [kSel, gcv] = selectLambdaGCV(yTrain(okTr), Ptr(okTr,:), R.edf, cfg);
        lamHow = 'GCV on the training fit (pre-Methods behaviour)';
        lamCurve = gcv;
        lamCurveName = 'GCV';
        lamShort = 'GCV on TRAIN';
    otherwise
        error('cfg.lambdaRule must be ''cv5'' or ''gcv''; got ''%s''.', cfg.lambdaRule);
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
% A lagged linear model applied over a whole trial IS a filter:
% so it needs one small mat-vec per lag per trial, NOT a materialised
% nT x (nUnits*nLags) design matrix per trial. For 100 test trials and 21 taps
predFull = predictFullTrials(spikes(:,:,testTrials), lags, b0, beta, nUnits);

% ---------------- DECODING INDEX vs the ZERO baseline, per test set --------
aIdx = find(traceTime >= cfg.analysisWin_s(1) & traceTime <= cfg.analysisWin_s(2));
[~, aZero] = min(abs(traceTime(aIdx)));
decIdx  = nan(nLA, nTS);
peakL   = nan(nLA, nTS);
nTrialL = zeros(nLA, nTS);
decTrial = cell(1, nTS);   % per-trial index, nLA x nTestTrials(ts)
for ts = 1:nTS
    if isempty(spec.testSets(ts).condIdx)
        sel = true(nTest,1);
    else
        member = unique(cell2mat(S.trialid(spec.testSets(ts).condIdx)'));
        sel = ismember(testTrials, member);
    end
    if ~any(sel), continue; end
    [decIdx(:,ts), peakL(:,ts), nTrialL(:,ts), decTrial{ts}] = ...   % 4th output
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
res(k).decTrial = decTrial;
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
% resampled to cfg.warpedSegmentLength points and dropped into its own slot, mean
% +/- 1.96 SEM, drawn segment by segment so nothing is joined across the gaps.
if wantWarp
if isempty(warpFig) || ~ishandle(warpFig)
    warpFig = figure('Name', sprintf('%s%s - warped lick train', passTag, spec.name), ...
        'Position', [110 110 1050 620]);
    warpTG  = uitabgroup(warpFig);
end
% ---- WHICH TRIALS THIS FIGURE SHOWS (display override, see spec.warpTestFirstN) ----
% Everything else on this tab and in the whole report uses testTrials. Only the
% warped figure honours the override, and only to make a clean example: the
% ridge fit above is already done and is not touched here.
if isempty(spec.warpTestFirstN)
    warpTrials = testTrials(:);
    warpPred   = predFull;
    warpSrc    = 'TEST trials';
else
    warpTrials = sort(poolTrials(:));
    warpTrials = warpTrials(1 : min(spec.warpTestFirstN, numel(warpTrials)));
    warpPred   = predictFullTrials(spikes(:,:,warpTrials), lags, b0, beta, nUnits);
    nInTrain   = sum(ismember(warpTrials, trainTrials));
    warpSrc    = sprintf('FIRST %d pool trials -- %d of them TRAINED ON (in-sample)', ...
        numel(warpTrials), nInTrain);
    logf('  [warp]  figure override: first %d pool trials (%d train, %d test) -- display only\n', ...
        numel(warpTrials), nInTrain, numel(warpTrials)-nInTrain);
end

wt = uitab(warpTG, 'Title', tabTitle);
[Wa, Wp, slotStart] = warpedLickMatrix(tongueLen(aIdx, warpTrials), warpPred(aIdx,:), aZero, cfg);
nShow = min(spec.warpShowLicks, cfg.maxLicksToExtract);

% ---- ONE PANEL PER CONDITION, STACKED ----
% With two test sets (R1 vs R4, R1 vs R6) the tab holds two axes, R1 on top, in
% the order spec.testSets declares them. Each panel is one condition: actual
% black, predicted in that condition's colour, labelled in-plot. Only the bottom
% panel carries the x axis. Both panels are cut from the SAME warped matrix, so
% the slot geometry is identical and the two are directly comparable.
if nTS >= 2 && cfg.warpSplitTestSets, panels = 1:nTS; else, panels = 1; end
nPan = numel(panels);
nStr = '';
for pp = 1:nPan
    ts = panels(pp);
    if nPan == 1 || isempty(spec.testSets(ts).condIdx)
        selTS = true(numel(warpTrials),1);
    else
        memberTS = unique(cell2mat(S.trialid(spec.testSets(ts).condIdx)'));
        selTS = ismember(warpTrials, memberTS);
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
    sprintf('%s  |  warped licks 1-%d  |  %s  |  mean \\pm 1.96 SEM  |  %s', hdr, nShow, warpSrc, strtrim(nStr)), ...
    'EdgeColor','none', 'HorizontalAlignment','center', ...
    'FontWeight','bold', 'FontSize',10, 'Interpreter','tex');

end   % wantWarp

clear spikes Xtrain Xtest Ptr Pte predFull Mact Mpred Wa Wp warpPred
end   % probe
end   % session
%% REPORT (shared)
nLA  = cfg.nLicksAnalyze;
nG   = numel(spec.groupMaps);
nTS  = numel(spec.testSets);
lickVec = (1:nLA)';

% COLOURS come from the study block, never from here -- each task has its own
% scheme and the shared code must not overrule it.
%   spec.tsColours   one row per TEST SET   (R1 vs R4, R1 vs R6, ...)
%   spec.grpColours  one row per GROUP, or [] for the default gradient. Only the
%                    learning script draws groups on one axes, and it sets this.
tsCols = spec.tsColours;
if isempty(spec.grpColours)
    grpCols = [linspace(0.10,0.85,max(nG,2))', zeros(max(nG,2),1), linspace(0.85,0.10,max(nG,2))'];
else
    grpCols = spec.grpColours;
end
assert(size(tsCols,1) >= nTS, 'spec.tsColours has %d rows for %d test sets.', size(tsCols,1), nTS);

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
% Every session with a finite contact-1 value is kept. There is no screen.
    for ts = 1:nTS, out.keep{g,ts} = isfinite(out.curves{g,ts}(1,:)); end
end

%% PER-CONTACT TEST: contact L vs contact 1
% at every contact from 2 up. Both a t-test and a Wilcoxon signed rank, exactly
% BH-FDR across contacts 2..nLicks is ALSO computed and printed. With 7 tests on
% the same sessions the uncorrected stars are optimistic; the q column is what to
out.pVsC1_t  = nan(nLA, nG, nTS);
out.pVsC1_sr = nan(nLA, nG, nTS);
out.qVsC1_sr = nan(nLA, nG, nTS);
fprintf('\n%s\n', repmat('=',1,84));
fprintf('CONTACT L vs CONTACT 1 | %s paired signed-rank | %s | %s\n', spec.vsC1Tail, spec.figRef, spec.name);
fprintf('%s\n', repmat('-',1,84));
for g = 1:nG
    for ts = 1:nTS
        C = out.curves{g,ts}(:, out.keep{g,ts});
        if size(C,2) < 2, continue; end
        for L = 2:nLA
            a = C(L,:);  b = C(1,:);
            ok = isfinite(a) & isfinite(b);
            if sum(ok) < 2, continue; end
            [~, out.pVsC1_t(L,g,ts)] = ttest(a(ok), b(ok), 'Tail','left');
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
            fprintf('C%d p=%.3f(q=%.3f)%s  ', L, out.pVsC1_sr(L,g,ts), out.qVsC1_sr(L,g,ts), mk);
        end
        fprintf('| n = %d\n', size(C,2));
    end
end
fprintf('  (* = BH-FDR q < %.2f; the FIGURE now stars these same q values. t-test column is a reference only.)\n', ...
    cfg.alpha);
fprintf('%s\n', repmat('=',1,84));

%% BETWEEN TEST SETS, per contact
% For the two-pool tasks (R1 vs R4, R1 vs R6): paired at each contact, TWO-tailed
% because there is no prior direction.
out.pBetween = nan(nLA, nG);
out.qBetween = nan(nLA, nG);
if nTS >= 2
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
            fprintf('C%d p=%.3f(q=%.3f)%s  ', L, out.pBetween(L,g), out.qBetween(L,g), mk);
        end
        fprintf('| n = %d, m = %d\n', sum(k), nLA);
    end
    fprintf('%s\n', repmat('-',1,84));
end

%% EARLY vs LATE
fprintf('\n%s\n', repmat('=',1,84));
fprintf('EARLY (contacts %s) vs LATE (contacts %s), one-tailed paired: LATE < EARLY\n', ...
    mat2str(cfg.earlyLicks), mat2str(cfg.lateLicks));
fprintf('%s\n', repmat('-',1,84));
out.pEarlyLate = nan(nG, nTS);
for g = 1:nG
    for ts = 1:nTS
        C = out.curves{g,ts}(:, out.keep{g,ts});
        if isempty(C), continue; end
        e = mean(C(cfg.earlyLicks,:), 1, 'omitnan')';
        l = mean(C(cfg.lateLicks, :), 1, 'omitnan')';
        ok = isfinite(e) & isfinite(l);
        if sum(ok) < 2, continue; end
        out.pEarlyLate(g,ts) = signrank(l(ok), e(ok), 'Tail','left');
        fprintf('  %-7s %-9s : early %+.3f  late %+.3f  diff %+.3f | signrank p = %.4f | n = %d\n', ...
            spec.groupLabels{g}, spec.testSets(ts).name, mean(e(ok)), mean(l(ok)), ...
            mean(l(ok)-e(ok)), out.pEarlyLate(g,ts), sum(ok));
    end
end
fprintf('%s\n', repmat('=',1,84));

%% GROUP A vs GROUP B, per contact
% Only when spec.compareGroups is set (the learning script sets it to the first
% and last day). Paired across sessions, one-tailed: group B lower than group A.
out.pGroupAB = nan(nLA, nTS);
out.pSummary = nan(1, max(nG, nTS));
out.qGroupAB = nan(nLA, nTS);
if ~isempty(spec.compareGroups)
    gA = spec.compareGroups(1);  gB = spec.compareGroups(2);
    fprintf('\n%s\n', repmat('=',1,84));
    fprintf('%s vs %s at each contact | one-tailed paired, H1: %s LOWER\n', ...
        spec.groupLabels{gA}, spec.groupLabels{gB}, spec.groupLabels{gB});
    for ts = 1:nTS
        A = out.curves{gA,ts};  B = out.curves{gB,ts};
        if size(A,2) ~= size(B,2)
            fprintf('  %s and %s have different session counts -- not paired.\n', ...
                spec.groupLabels{gA}, spec.groupLabels{gB});
            break
        end
        for L = 1:nLA
            ok = isfinite(A(L,:)) & isfinite(B(L,:));
            if sum(ok) < 2, continue; end
            out.pGroupAB(L,ts) = signrank(B(L,ok), A(L,ok), 'Tail','left');   % secondary readout only
        end
        out.qGroupAB(:,ts) = bhFDR(out.pGroupAB(:,ts));
        fprintf('  %-9s : ', spec.testSets(ts).name);
        for L = 1:nLA
            mk = ' ';
            if isfinite(out.qGroupAB(L,ts)) && out.qGroupAB(L,ts) < cfg.alpha, mk = '*'; end
            fprintf('C%d p=%.3f%s ', L, out.pGroupAB(L,ts), mk);
        end
        fprintf('| n = %d\n', size(A,2));
    end
    fprintf('%s\n', repmat('=',1,84));
end

%% THE TEST THE PANEL REPORTS, AND THE SUMMARY FILE
% One value per session, exactly as the Methods define it:
%   "we tested the change in the measured quantity relative to a within-session
%    reference rather than its absolute value, as baseline levels varied
%    somewhat across sessions and individuals. Each session contributed a single
%    from the first lick, averaged over licks 3 through 8 for the comparison
%    between the Simple and Double Reward Tasks and between the first and fifth
%    day of training."
% so  value(session) = mean over spec.groupSummaryLicks of
%                      ( decIdx(L) - decIdx(1) ).
% Contact 1 is the within-session reference, which is why it is zero for every
% session by construction and why the absolute level does not enter.
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
% Fig. 5G / S4J: day 1 against day 5. INDEPENDENT groups, one-tailed for a
% reduction on day 5, ONE test per region, NOT corrected -- because each
% session contributes one value, so there is one comparison to make.
% Rank-sum rather than signed-rank, and the reason is the sample size: an
% exact signed-rank on 4 pairs has 2^4 = 16 sign assignments and so cannot
% return a two-sided p below 0.125. The rank-sum pools the two groups,
% giving C(n1+n2, n1) arrangements, and the Methods already state that
% day-1/day-5 sessions are treated as independent groups for this reason.
    gA = spec.compareGroups(1);  gB = spec.compareGroups(2);
    for ts = 1:nTS
        a = out.summary{gA,ts};  a = a(isfinite(a));
        b = out.summary{gB,ts};  b = b(isfinite(b));
        if numel(a) < 2 || numel(b) < 2
            fprintf('  test not run: n = %d vs %d.\n', numel(a), numel(b));
            continue
        end
        p1 = ranksum(a, b, 'tail', 'right');   % day 1 > day 5
        p2 = ranksum(a, b);
        out.pSummary(ts) = p1;
        fprintf('\n  Wilcoxon rank-sum, %s vs %s, ONE test, no correction\n', ...
            spec.groupLabels{gA}, spec.groupLabels{gB});
        fprintf('    one-tailed (%s > %s) p = %.4f%s\n', spec.groupLabels{gA}, ...
            spec.groupLabels{gB}, p1, repmat('  *', 1, p1 < cfg.alpha));
        fprintf('    two-tailed reference       p = %.4f\n', p2);
        fprintf('    exact floor for %d vs %d = %.4f\n', numel(a), numel(b), ...
            1 / nchoosek(numel(a)+numel(b), numel(a)));
    end

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
            p2 = ranksum(a, b);
            out.pSummary(g) = p1;
            fprintf('  %-8s: Simple n = %2d (%+0.4f) vs Double n = %2d (%+0.4f)\n', ...
                spec.groupLabels{g}, numel(a), mean(a), numel(b), mean(b));
            fprintf('            one-tailed rank-sum p = %.4f%s | two-tailed %.4f | floor %.4f\n', ...
                p1, repmat('  *', 1, p1 < cfg.alpha), p2, ...
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

% ---- write this script's summary for whoever needs it ----
summaryByGroup = cell(1, nG);
for g = 1:nG, summaryByGroup{g} = out.summary{g,1}; end
figRef   = spec.figRef;   %#ok<NASGU>
specName = spec.name;   %#ok<NASGU>
sessionsByGroup = out.sessIdx;   %#ok<NASGU>
save(spec.summaryFile, 'summaryByGroup', 'sessionsByGroup', 'figRef', 'specName');
logf('summary written to %s\n\n', spec.summaryFile);

%% DAY A vs DAY B INSIDE EACH ANIMAL, trials as units
% One test per animal: that animal's Day-A test trials against its Day-B test
% trials, Wilcoxon rank-sum on the per-trial drop from lick 1 (cfg.perAnimalLicks),
% one-tailed, H1: Day B drops MORE (Day A > Day B). Test set 1 ('all'), the same
% set the figures and the session-level test use. Not paired: the trials of the
% two days are different trials. The session-level rank-sum above is unchanged.
out.PA = struct('anm',{{}}, 'nA',[], 'nB',[], 'medA',[], 'medB',[], 'p',[], 'k',0, 'n',0);
if cfg.perAnimalTest && ~isempty(spec.compareGroups)
    gA = spec.compareGroups(1);  gB = spec.compareGroups(2);
    rowsA = find(arrayfun(@(r) spec.groupMaps{gA}(r.session) == r.probe, res));
    rowsB = find(arrayfun(@(r) spec.groupMaps{gB}(r.session) == r.probe, res));
    anmA  = arrayfun(@(r) char(string(r.anm)), res(rowsA), 'UniformOutput', false);
    anmB  = arrayfun(@(r) char(string(r.anm)), res(rowsB), 'UniformOutput', false);
    anmAll = unique([anmA(:); anmB(:)], 'stable');
    nAn = numel(anmAll);
    PA  = struct('anm',{anmAll}, 'nA',zeros(nAn,1), 'nB',zeros(nAn,1), ...
                 'medA',nan(nAn,1), 'medB',nan(nAn,1), 'p',nan(nAn,1), 'k',0, 'n',0);
    for ai = 1:nAn
        iA = rowsA(strcmp(anmA, anmAll{ai}));
        iB = rowsB(strcmp(anmB, anmAll{ai}));
        dA = zeros(0,1);  dB = zeros(0,1);
        if ~isempty(iA), dA = trialDrop(res(iA(1)).decTrial{1}, cfg.perAnimalLicks, cfg.perAnimalMinLicks); end
        if ~isempty(iB), dB = trialDrop(res(iB(1)).decTrial{1}, cfg.perAnimalLicks, cfg.perAnimalMinLicks); end
        PA.nA(ai) = numel(dA);   PA.nB(ai) = numel(dB);
        if ~isempty(dA), PA.medA(ai) = median(dA); end
        if ~isempty(dB), PA.medB(ai) = median(dB); end
        if numel(dA) >= 2 && numel(dB) >= 2
            PA.p(ai) = ranksum(dA, dB, 'tail', 'right');   % Day A > Day B
        end
    end
    PA.k = sum(PA.p < cfg.alpha);
    PA.n = sum(isfinite(PA.p));
    out.PA = PA;

    fprintf('\n%s\n', repmat('=',1,84));
    fprintf('%s vs %s INSIDE EACH ANIMAL | Wilcoxon rank-sum on test trials\n', ...
        spec.groupLabels{gA}, spec.groupLabels{gB});
    fprintf('per-trial drop = mean(decoding index, licks %s) - index at lick 1 | >= %d of those licks\n', ...
        mat2str(cfg.perAnimalLicks), cfg.perAnimalMinLicks);
    fprintf('one-tailed, H1: %s drops MORE than %s\n', spec.groupLabels{gB}, spec.groupLabels{gA});
    fprintf('%s\n', repmat('-',1,84));
    fprintf('  %-12s %8s %8s %12s %12s %10s\n', 'animal', ['n ' spec.groupLabels{gA}], ...
        ['n ' spec.groupLabels{gB}], ['med ' spec.groupLabels{gA}], ['med ' spec.groupLabels{gB}], 'p');
    for ai = 1:nAn
        st = '';  if isfinite(PA.p(ai)) && PA.p(ai) < cfg.alpha, st = '  *'; end
        fprintf('  %-12s %8d %8d %+12.4f %+12.4f %10.4g%s\n', PA.anm{ai}, PA.nA(ai), PA.nB(ai), ...
            PA.medA(ai), PA.medB(ai), PA.p(ai), st);
    end
    fprintf('  -> significant in %d of %d animals (p < %.3g)\n', PA.k, PA.n, cfg.alpha);
    if PA.n < nAn
        fprintf('  (%d animal(s) not tested: fewer than 2 usable trials on one of the days)\n', nAn - PA.n);
    end
    fprintf('%s\n', repmat('=',1,84));
end

%% FIGURE 1: DECODING INDEX, TABBED
if spec.showFigCI
% Two tabs on one panel: the mean with its CI band, and the per-session spread.
% The amplitude and trial-count panels are gone; the amplitude caveat still
% trial counts are printed by the fit loop.
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
                col = pickColour(spec, g, ts, grpCols, tsCols);
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
            drawStars(ax, 2:nLA, [], out.qVsC1_sr(:,gl,1), cfg.alpha);
        end
        ttl = sprintf('%s | %s | zero baseline | mean \\pm 95%% CI', spec.name, figNames{fi});
        sub = sprintf('* = BH-FDR q < %.2f, %s signed-rank, contact L vs contact 1', cfg.alpha, spec.vsC1Tail);
        title(ax, {ttl, sub}, 'FontSize',9, 'FontWeight','normal');

% ---- TAB: mean tongue length inside the fitting window ----
% only the samples inside ITS OWN window (variable length), so the trace is
% averaged over however many trials cover each sample; the count is shown on
% the right axis and thins out at the edges.
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
% Requested for the learning script and harmless everywhere else: the same means
% with nothing else drawn, so the group ordering is readable.
if spec.overlayGroups
    fh = figure('Name', sprintf('%s%s - mean traces only', passTag, spec.name));
    ax = axes('Parent', fh); hold(ax,'on')
    [h, lbl] = deal(gobjects(0), {});
    for g = 1:nG
        for ts = 1:nTS
            C = out.curves{g,ts};  k = out.keep{g,ts};
            if isempty(C), continue; end
            col = pickColour(spec, g, ts, grpCols, tsCols);
            h(end+1) = plot(ax, lickVec, mean(C(:,k),2,'omitnan'), '-o', 'Color', col, ...
                'MarkerFaceColor', col, 'LineWidth', 2.5);   %#ok<SAGROW>
            lbl{end+1} = seriesLabel(spec, g, ts, sum(k));   %#ok<SAGROW>
        end
    end
    finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{zero}', h, lbl, [0 1]);
    title(ax, sprintf('%s | mean decoding index, no error bars', spec.name), ...
        'FontSize',10,'FontWeight','normal');
end

%% FIGURE 3: DELTA FROM CONTACT 1, TABBED
if spec.showFigDelta
% Every session starts at exactly 0 by construction, so the band at contact 1
% collapses to a point. This is the quantity the contact-L-vs-contact-1 test is
% run on, so the stars here are the same test as on Figure 1.
% ONE FIGURE PER GROUP unless spec.overlayGroups -- regions are never drawn on
% the same axes, so M1 and ALM never share a plot. The learning script overlays
% because there the groups are learning days and comparing them IS the point.
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
                col = pickColour(spec, g, ts, grpCols, tsCols);
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
            drawStars(ax, 1:nLA, [], out.qGroupAB(:,1), cfg.alpha);   % BH q, per-contact secondary
            if cfg.perAnimalTest && out.PA.n > 0
                text(ax, 0.02, 0.98, sprintf('licks %d-%d, rank-sum within animal: %d/%d animals', ...
                    cfg.perAnimalLicks(1), cfg.perAnimalLicks(end), out.PA.k, out.PA.n), ...
                    'Units','normalized', 'VerticalAlignment','top', 'FontSize',8, 'Color',[0.3 0.3 0.3]);
            end
        end
        title(ax, {sprintf('%s | %s | delta from contact 1 | mean \\pm 95%% CI', spec.name, dNames{di}), ...
                   'contact 1 is 0 for every session by construction'}, ...
            'FontSize',9,'FontWeight','normal');

        ax = axes('Parent', uitab(tg, 'Title', 'per session'));  hold(ax,'on')
        [h, lbl] = deal(gobjects(0), {});
        for g = dg
            C = out.curves{g,1};  k = out.keep{g,1};
            if isempty(C), continue; end
            D = C - C(1,:);
            col = pickColour(spec, g, 1, grpCols, tsCols);
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
fprintf('\n%s\n', repmat('=',1,101));
fprintf('%-5s %-8s %-12s %5s %-8s %6s %6s %6s %10s %7s %9s %9s %11s\n', ...
    'sess','animal','date','probe','region','units','nTrn','nTst','lambda','edf', ...
    'R2 train','R2 test','win samp');
fprintf('%s\n', repmat('-',1,101));
for r = res
    wl = r.winLenSamples;
    fprintf('%-5d %-8s %-12s %5d %-8s %6d %6d %6d %10.3g %7.1f %9.3f %9.3f %5d-%-5d\n', ...
        r.session, r.anm, r.date, r.probe, r.region, r.nUnits, r.nTrain, r.nTest, ...
        r.lambda, r.edf, r.r2TrSel, r.r2TeSel, min(wl), max(wl));
end
fprintf('%s\n', repmat('=',1,101));
%% DAY A vs DAY B, animal level (this script only)
% The per-contact version is computed in the report (spec.compareGroups) and its
% stars are on the delta figure. This adds the ANIMAL-LEVEL summary: one mean
% over cfg.lateLicks per animal, so n = 4 independent pairs.
% At n = 4 a one-tailed signrank floors at 1/2^4 = 0.0625 -- it can corroborate
% but never independently confirm, so the direction count beside it is the honest
% summary at this sample size. Run on the DELTAS, which is what the figure shows.
gA = spec.compareGroups(1);  gB = spec.compareGroups(2);
animalsA = {res(arrayfun(@(r) spec.groupMaps{gA}(r.session) == r.probe, res)).anm};
animalsB = {res(arrayfun(@(r) spec.groupMaps{gB}(r.session) == r.probe, res)).anm};
pairedOK = numel(animalsA) == numel(animalsB) && all(strcmp(animalsA(:), animalsB(:)));
if ~pairedOK
% possibility rather than a bug: day 1 is early learning, and an animal that
% never produced 8 contacts in a row on its first day loses every trial, the
% session is skipped for too few trials, and the two days stop lining up
% animal for animal. Say so and skip the test instead of aborting a run that
% has already produced every figure.
    fprintf('\n%s\n', repmat('=',1,84));
    fprintf('%s and %s do not line up animal-for-animal -- animal-level test SKIPPED.\n', ...
        spec.groupLabels{gA}, spec.groupLabels{gB});
    fprintf('  %s: %s\n  %s: %s\n', spec.groupLabels{gA}, strjoin(animalsA,', '), ...
        spec.groupLabels{gB}, strjoin(animalsB,', '));
    if cfg.requireBout && spec.applyBoutFilter
        fprintf('  MOST LIKELY CAUSE: the bout filter (>= %d contacts in a row within %.2f s).\n', ...
            cfg.boutMinLicks, cfg.boutWin_s);
        fprintf('  Early in learning animals may not bout at all, so those sessions lost every\n');
        fprintf('  trial. Set cfg.requireBout = false, or loosen cfg.boutMinLicks / cfg.boutWin_s,\n');
        logf('  and check the [bout] lines in the fit log for the day-1 sessions.\n');
    end
    fprintf('%s\n', repmat('=',1,84));
end
if pairedOK
CA = out.curves{gA,1};  CB = out.curves{gB,1};
DA = CA - CA(1,:);      DB = CB - CB(1,:);
mA = mean(DA(cfg.lateLicks,:),1,'omitnan')';
mB = mean(DB(cfg.lateLicks,:),1,'omitnan')';
gd = isfinite(mA) & isfinite(mB);
fprintf('\n%s\n%s vs %s | DELTA from contact 1, contacts %s | one-tailed paired, H1: %s LOWER\n%s\n', ...
    repmat('=',1,84), spec.groupLabels{gA}, spec.groupLabels{gB}, mat2str(cfg.lateLicks), ...
    spec.groupLabels{gB}, repmat('-',1,84));
fprintf('  animals: %s\n', strjoin(animalsA,', '));
if sum(gd) >= 2
    fprintf('  animal-level signrank p = %.4f (n=%d) | %d/%d animals lower | %s %+.3f, %s %+.3f\n', ...
        signrank(mB(gd), mA(gd), 'Tail','left'), sum(gd), sum(mB(gd) < mA(gd)), sum(gd), ...
        spec.groupLabels{gA}, mean(mA(gd)), spec.groupLabels{gB}, mean(mB(gd)));
    fprintf('  (0.0625 is the floor a one-tailed signrank can reach at n = 4)\n');
else
    fprintf('  insufficient paired data (n=%d)\n', sum(gd));
end
fprintf('%s\n', repmat('=',1,84));
end   % pairedOK


end   % repIdx -- REPEAT PASSES

%% LOCAL FUNCTIONS (shared)

function r = tongueRuns(v, cfg)
% Protrusions in one trial's tongue length: [onset offset] sample pairs, one row
% per protrusion. A protrusion is a contiguous run of VISIBLE tongue (length > 0;
% after the zero-fill "not visible" is exactly 0, so no threshold heuristic).
% TWO THINGS HAPPEN HERE, IN THIS ORDER, AND THE ORDER MATTERS:
%   1. runs separated by <= cfg.mergeGapSamples zeros are MERGED. A one- or
%      two-frame DLC dropout mid-protrusion is not the tongue going in and out
%      again; leaving it unmerged ends the fitting window at the dropout and
%      splits one lick into two on every downstream count.
%   2. only then is the minimum-duration rule applied. Applying it first would
%      delete the two halves of a dropout-split protrusion before they could be
%      rejoined.
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
% Which protrusion is contact 1. t = 0 IS the first post-go-cue port contact, and
% the tongue is necessarily OUT at that instant, so contact 1 is the protrusion
% that CONTAINS alignSample.
% The old rule -- the run whose ONSET is nearest t = 0 -- picks the wrong lick
% whenever a pre-cue protrusion happens to start closer to 0 than the real one
% does, and it silently reports that neighbour's retraction as the end of contact
% 1. Containment is the definition; nearest-onset is only the fallback for when
% tracking dropped out across t = 0 itself.
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
% Sample indices of ONE trial's fitting window, in one of two modes.
%           + cfg.postTongueRetract_s
%   The end is the RETRACTION, not the contact: contact 1 happens part way
%   still labelled before it went back in the mouth.
%   where the bout-end contact is the first one not followed by another within
%   cfg.winBoutGap_s. contactRel holds the post-go-cue contact times relative to
%   contact 1, so contactRel(1) is 0 by construction.
%   This window is defined ENTIRELY by contact times -- it never looks at the
%   tongue trace. That is deliberate: the whole point of the mode is to fit the
%   late contacts, and a kinematic rule would smuggle the same protrusion-shaped
%   assumption back in. A trial with fewer than two contacts, or whose bout ends
%   at contact 1, has no such window and is dropped.
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
% END OF LICK 2  ->  END OF BOUT.
% The start is the RETRACTION of the protrusion that carries contact
% 2 -- the last frame the tongue is still labelled before it goes
% back in the mouth -- plus the pad. NOT contact 2's time: the
% contact happens part way up the protrusion, exactly as it does for
% contact 1, so anchoring on it would put the window start inside
% lick 2 rather than after it.
% Which protrusion carries contact 2 is decided by the same
% containment rule the contact-1 window uses (pickFirstRun), just
% evaluated at contact 2's sample instead of at t = 0. So all three
% modes agree about what a lick is and where it ends.
            if numel(contactRel) < 2, return; end
            kEnd = boutEndContact(contactRel, cfg.winBoutGap_s);
            if kEnd < 2, return; end   % the bout is over before contact 2
            r = tongueRuns(tlTrial, cfg);
            if isempty(r), return; end
            [~, s2] = min(abs(traceTime - contactRel(2)));   % sample at contact 2
            i2 = pickFirstRun(r, s2);   % its protrusion
            tStart = traceTime(r(i2,2)) + cfg.winBoutStartPad_s;   % ITS RETRACTION
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
% whose gap to the NEXT contact exceeds gap_s. If no gap ever exceeds it, the
% bout runs to the last contact of the trial.
% Read as a stopping rule, not as "the last lick in the trial". Every trial's
% final lick is trivially not followed by another one, so taking that literally
% would make the rule vacuous and would sweep in a stray late lick seconds after
% the bout ended. What the gap identifies is the end of the CONTINUOUS train.
    d = diff(contactRel(:));
    k = find(d > gap_s, 1);
    if isempty(k), k = numel(contactRel); end
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
% COLUMN ORDER: lag-major, unit-minor; column (li-1)*nUnits + u is unit u at
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
% THE METHODS NAME THE IMPLEMENTATION, so this function honours that:
%   "To prevent overfitting, weights were estimated by ridge regression
%    (using the MATLAB function 'ridge')"
% algebraically identical closed form below, which is far faster when there are
% many predictors. Their agreement is checked at runtime -- see
% cfg.ridgeCheckEquivalence and checkRidgeEquivalence.
% ---- THE ARGUMENT THAT IS NOT OPTIONAL ----------------------------------
% The trailing 0 asks for coefficients ON THE SCALE OF THE ORIGINAL PREDICTORS,
% and the returned matrix is (p+1) x numel(lambdas) with the INTERCEPT IN ROW 1.
% that -- the two sets differ by whatever the predictor scales happen to be,
% several-fold on this data -- so the mistake is silent and large.
% ---- WHAT IS PENALIZED --------------------------------------------------
% ridge always penalizes the STANDARDIZED coefficients and always leaves the
% intercept unpenalized. That is the same estimator as the closed form below
% with standardizeCols = true, so standardizeCols = false has no counterpart in
% ridge and is rejected rather than quietly ignored.
% ---- EFFECTIVE DEGREES OF FREEDOM ---------------------------------------
% R.edf(k) = 1 + sum_i d_i/(d_i + lambda_k), over the eigenvalues d of the
% standardized Gram matrix; the 1 is the free intercept. It is the trace of the
% hat matrix. ridge does not return it, so it is computed from the singular
% values of the centred, scaled design -- ONE svd, whatever the number of
% lambdas. It is needed only by the 'gcv' lambda rule and by the fit printout,
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
% VERIFIED, not asserted: run head to head against ridge() on 8
% probe fits (35-session R1/R4 set, 650-1800 predictors), the two
% chose the same lambda every time and the held-out predictions
% agreed to a median of 1e-15 and a worst single sample of 1.8e-14
% on the 0-1 kinematic scale -- about 1e-12 of one camera pixel,
% and eleven orders below the third decimal the decoding index is
% reported to. Measured speedup ~30x per fit.
            mx = mean(X,1);  my = mean(y);
            Xc = X - mx;     yc = y - my;
            if standardizeCols
% replaces the standard deviation by 1 whenever
% zero. A unit firing in one or two window samples can land
% between those two tests, and then `== 0` divides that column
% by ~1e-10 while ridge() divides it by 1 -- a huge coefficient
% difference that has nothing to do with the algebra. Matching
% the threshold is what makes this a drop-in replacement.
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
% centred, scaled design ARE the Gram eigenvalues. The svd below
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
% at double precision. If this ever prints something larger, the fast path is
% wrong for that data and cfg.ridgeImpl should go back to 'matlab'.
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

function [kSel, gcv] = selectLambdaGCV(y, P, edf, cfg)
% Lambda by generalized cross-validation on the TRAINING FIT ONLY.
% No validation split, no folds, no held-out rows: RSS is the in-sample residual
% sum of squares and edf is the trace of the hat matrix, so the denominator is
% what stops the smallest lambda from always winning. Lambdas whose edf reaches n
% are unusable (the denominator collapses) and are excluded rather than allowed
% to produce an infinite score.
% cfg.lambdaFixed pins lambda instead; the nearest grid point is used.
    n   = numel(y);
    nK  = size(P, 2);
    gcv = inf(nK, 1);
    for kk = 1:nK
        den = n - edf(kk);
        if den <= 0, continue; end
        gcv(kk) = n * sum((y - P(:,kk)).^2) / den^2;
    end
    if isempty(cfg.lambdaFixed)
        [~, kSel] = min(gcv);
    else
        [~, kSel] = min(abs(log10(cfg.ridgeGrid(:)) - log10(cfg.lambdaFixed)));
    end
end

function [Wa, Wp, slotStart, traceLen] = warpedLickMatrix(Yact, Ypred, alignSample, cfg)
% Every protrusion is resampled onto cfg.warpedSegmentLength points and dropped
% into its own SLOT on a synthetic axis: lick k occupies
% with slotStart = warpedSlotFirst : warpedSlotSpacing : ... . The spacing is
% wider than the segment, so the gaps between slots stay NaN and the licks read
% as separate events instead of one connected line. Same three constants as
% before: segment 30, first slot 10, spacing 50.
% Resampling is interp1(..., 'pchip', NaN) over the FINITE samples of the segment,
% and a segment with a single finite sample is placed at the centre of its slot --
% both as in the original. Predictions are warped through the SAME sample indices
% as the actual trace, never re-detected, so the two curves are aligned segment
% WHY WARP AT ALL. Later licks in a bout are shorter as well as smaller. On a real
% time axis that difference smears the excursions together and flattens the mean
% what is left to compare is amplitude and shape.
% The one thing that is NOT carried over from the original is which run counts as
% lick 1: that now comes from the shared tongueRuns / pickFirstRun pair, so the
% warped figure, the fitting window and the decoding index all mean the same thing
% by "contact 1".
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
% One excursion resampled onto segL points -- the original's rule, with the
    wt = nan(segL,1);
    xv = v(isfinite(v));
    if numel(xv) == 1
        wt(ceil(segL/2)) = xv;
    elseif numel(xv) > 1
        wt = interp1(linspace(0,1,numel(xv))', xv(:), linspace(0,1,segL)', method, NaN);
    end
end

function drawWarpedSlotPanel(ax, Wa, Wp, slotStart, nShow, colP, nameStr, showX, cfg)
% ONE CONDITION ON ONE AXES: actual in black, predicted in that condition's
% colour, and the two labelled with TEXT ON THE AXES rather than a legend box --
% a key costs a corner of the plot and an eye movement, and with one condition
% per panel there is nothing to disambiguate.
% The labels sit upper-right because that is where the trace has room once the
% bout decays; if a session's licks stay large to the end they may land on data,
% and the two text() lines below are the place to move them.
    hold(ax, 'on')
    drawWarpedSet(ax, Wa, [0 0 0], slotStart, cfg, nShow);
    drawWarpedSet(ax, Wp, colP,    slotStart, cfg, nShow);
    if cfg.warpedTickCentre, tickAt = slotStart + cfg.warpedSegmentLength/2; else, tickAt = slotStart; end
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
% SEGMENT BY SEGMENT so no line is drawn across the gap between two licks.
% Normal 1.96 rather than a t quantile, matching the original figure.
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
% only ever targets a figure, so the position is computed the same way subplot
% does and the axes is parented explicitly.
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

function [decIdx, peakL, nTrialL, trialIdx] = contactDecodingIndex(Yact, Ypred, alignSample, cfg)
% 4th output trialIdx (nLA x nTrials): the same index for ONE trial,
%   1 - RMSE_decoded(trial, L) / RMSE_zero(session, L)
% i.e. the trial's own decoding error against the session-pooled null at that
% pooled across trials. RMSE_zero is the error a constant-zero predictor makes.
% Contact runs come from tongueRuns (dropout-merged, minimum duration applied)
% and contact 1 is the protrusion CONTAINING alignSample -- the same two rules
% the fitting window uses, so "contact 1" means one thing in this whole script.
    nLA = cfg.nLicksAnalyze;
    nTr = size(Yact, 2);
    sseDec = zeros(nLA,1);  sseZero = zeros(nLA,1);
    nTrialL = zeros(nLA,1); peaks = nan(nLA, nTr);
    nSampL  = zeros(nLA,1);
    sseDecT = nan(nLA, nTr);  nSampT = zeros(nLA, nTr);
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
            nSampL(L)    = nSampL(L) + sum(g);
            sseDecT(L,tt) = sum((a(g) - p(g)).^2);
            nSampT(L,tt)  = sum(g);
            peaks(L,tt) = max(a(g));
        end
    end
    decIdx = 1 - sqrt(sseDec ./ sseZero);   % n cancels inside the ratio
    decIdx(nTrialL == 0) = NaN;
    peakL = mean(peaks, 2, 'omitnan');
% per-trial index against the session-pooled null RMSE
    rmseZero = sqrt(sseZero ./ max(nSampL,1));   % nLA x 1
    rmseZero(nSampL == 0 | rmseZero <= 0) = NaN;
    trialIdx = 1 - sqrt(sseDecT ./ max(nSampT,1)) ./ rmseZero;   % implicit expansion
    trialIdx(nSampT == 0) = NaN;
end

function d = trialDrop(M, licks, minN)
% one value per trial: mean(index over licks) - index at lick 1, for
% trials with a lick-1 value and at least minN finite values among licks.
    d = zeros(0,1);
    if isempty(M), return; end
    licks = licks(licks >= 2 & licks <= size(M,1));
    if isempty(licks), return; end
    base = M(1,:);
    late = M(licks,:);
    n    = sum(isfinite(late), 1);
    v    = mean(late, 1, 'omitnan') - base;
    d    = v(isfinite(base) & n >= minN);
    d    = d(:);
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

function col = pickColour(spec, g, ts, grpCols, tsCols)
% Colour by GROUP when groups share an axes (the learning script), otherwise by
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
% is the fixed range every decoding-index axis uses.
% version-fragile, line() works everywhere.
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
% One asterisk per contact that survives Benjamini-Hochberg at q < alpha.
% the signed rank -- both at UNCORRECTED p < 0.05, while the BH q values went
% only to the console. So the published panel and the stated correction did not
% agree, and half the stars came from a test the Methods do not use. The third
% argument is ignored and kept only so old call sites still run.
    yl  = ax.YLim;
    pad = 0.04 * range(yl);
    for L = contacts(:)'
        if L > numel(q), continue; end
        if isfinite(q(L)) && q(L) < alpha
            text(ax, L, yl(2) - 1.0*pad, '*', 'Color','r', 'FontSize',16, ...
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
%   S.time, S.trialdat, S.nUnitsTotal, S.tongueRaw, S.cluid, S.trialid,
%   S.anm, S.date, S.Ntrials, S.lickL, S.goCue, S.trialTypes, S.hit
% WHY THIS EXISTS. loadSessionData + getKinematics + loadMotionEnergy is where
% essentially all the wall-clock time goes -- the ridge fit itself takes well
% under a second. Those three are deterministic given (loader, date, params), and
% a lag or a penalty and re-run, the old code paid the full load again for
% nothing. With the cache warm, a re-run skips straight to the fit.
% CACHE KEY. The full set of params fields that can change what gets loaded is
% written into the .mat and compared on read, so a stale cache can never be used
% silently: change params.smooth or params.quality or the condition list and the
% whole point.
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

function [kSel, sse] = selectLambdaCV5(X, y, rowTrial, cfg)
% Lambda by k-fold cross-validation WITHIN THE TRAINING SET, which is what the
% Methods describe:
%    predict the held-out fold, pooling squared errors across folds. The decoder
% Squared errors are POOLED across folds -- summed, not averaged per fold --
% so folds that happen to contain longer trials carry proportionally more
% weight, which is what "pooling squared errors" means and what makes the
% criterion a plain held-out SSE.
% FOLDS ARE OVER WHOLE TRIALS, NOT ROWS. Each row of X is one time sample, and
% consecutive samples from the same trial are strongly autocorrelated. Splitting
% rows at random would leave a held-out sample's immediate neighbours sitting in
% the training set, so the held-out error would measure interpolation rather
% than generalisation and would keep choosing a lambda that is too small. That
% is exactly the failure mode of the GCV rule this replaces. rowTrial comes from
% laggedDesign and says which trial each row belongs to.
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

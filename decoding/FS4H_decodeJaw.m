%% FS4H_decodeJaw.m
%  Ridge decoder of jaw position from population spike rates, VTA Reward Task.
%  Spike rates are z-scored and lagged, then fit per session with a ridge
%  penalty chosen by cross-validation over held-out training trials. The
%  decoding index is 1 - RMSE(decoded) / RMSE(zero baseline), computed for
%  each lick contact and averaged across sessions.
%  Produces Fig. S4H.
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
RUN.cvFolds  = 5;   % Methods: five-fold cross-validation (was 3), as the tongue scripts
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

% FITTING WINDOW -- per trial, variable length.
%           + postJawClose_s
% THE JAW ANALOGUE OF THE TONGUE RULE, and it has to be, because jaw
% displacement has no zero. For tongue the window closed when the tongue went
% back in the mouth -- the last visible sample of the protrusion. The jaw never
% disappears; what corresponds to "back in the mouth" is the jaw CLOSING, and
% the closing trough of contact 1's cycle is where that happens. Both rules say
% the same thing: run the window to the end of the first movement, not to the
% contact, which happens part way through it.
% The trough comes from localJawSegments (the persistence-based cycle finder
% index and the warped figure all inherit ONE definition of "cycle 1".
cfg.preGoCue_s     = 0.20;   % window opens this long BEFORE the go cue
cfg.postJawClose_s = 0.005;   % window closes this long AFTER contact 1's cycle does
cfg.minWindowSamples = 10;   % trials with a shorter window are dropped

%% WHICH FITTING WINDOW (three options)
% 'contact1'        THE DEFAULT, and what the block above describes:
%                   go cue - preGoCue_s  ->  contact 1's jaw cycle CLOSES,
%                   + postJawClose_s.
% 'bout'            From winBoutStartPad_s AFTER CONTACT 2, to winBoutEndPad_s
%                   after the contact that ENDS the bout -- the first one not
%                   followed by another within winBoutGap_s. Defined ENTIRELY
%                   from contact times; it never looks at the jaw trace, so it
%                   is identical to the tongue script's 'bout' window and the
%                   two variables are fitted over exactly the same samples.
% 'lick2ToBoutEnd'  END OF LICK 2 -> END OF BOUT. Same end as 'bout', but the
%                   start is the CLOSING TROUGH of cycle 2 + winBoutStartPad_s,
%                   not contact 2's time. The jaw analogue of "after lick 2 has
%                   finished": contact 2 happens part way up the opening, so
%                   anchoring on it would put the start inside lick 2.
% the FIRST cycle, so contacts 2-8 are extrapolation -- a decline across the bout
% could be the model leaving its training distribution rather than cortex
% disengaging, and the two explanations predict the same curve. Move the window
% to the late contacts and they diverge: if the decline survives it is not a
cfg.fitWindowMode     = 'contact1';   % EARLY model = Fig. S4H (was 'contact1'). Options: 'contact1' | 'bout' | 'lick2ToBoutEnd'
cfg.winBoutStartPad_s = 0.02;   % opens this long after contact 2 / cycle 2 closing
cfg.winBoutEndPad_s   = 0.20;   % closes this long AFTER the bout's last contact
cfg.winBoutGap_s      = 0.75;   % a gap longer than this is what ENDS the bout

% NEURAL Z-SCORING WINDOW
cfg.zscoreWin_s = [-2.0 2.0];   % Methods: -2 s to +2 s around the first port contact
cfg.minBaselineStd = 1e-6;   % units below this are DROPPED, not divided by

% NEURAL LAG CONTEXT (lag set is -preBins : +postBins, inclusive)
cfg.lagPre_s  = 0.06;
cfg.lagPost_s = 0.04;

% ANALYSIS WINDOW for the per-contact decoding index
cfg.analysisWin_s = [-0.09 3.0];

%% TARGET: jaw displacement
% THE ONE PLACE WHERE JAW AND TONGUE GENUINELY DIVERGE. Read it before changing
% anything here.
% For tongue_length a NaN means "tongue not visible", which physically means
% length ZERO. It is a real observation, so the tongue scripts zero-fill it and
% score it like any other sample.
% jaw_ydisp_view1 IS NOT LIKE THAT. A NaN means the tracker lost the jaw marker.
% is a fabricated observation, and it breaks three things -- all three are the
%   1. TRAINING. Rows with a fabricated 0 target teach the decoder to output ~0
%      whenever tracking drops out.
%   2. LAMBDA. If the prediction is NaN-masked at a dropout but the target is a
%      real 0 there, SSres drops the sample and SStot keeps it -- two different
%      sample sets, so R^2 rises with the dropout rate.
%   3. THE DECODING INDEX. The same asymmetry, so more tracking failure looks
%      like better decoding.
% Here NaN stays NaN: training rows with a NaN target are dropped by the same
% okTr mask that drops NaN predictors, and every score -- the ridge R^2, the
% decoded RMSE and the baseline RMSE -- is taken over ONE common finite mask
% across actual AND predicted, so numerator and denominator always cover
% identical samples. Predictions are never masked: the neural data is fine, so
% the model can predict everywhere, which is also what keeps the full-trial
% trace continuous. Segmentation runs on a gap-filled COPY (findpeaks cannot
% take NaN) and the indices it returns are applied to the real NaN-bearing
% arrays, so nothing interpolated is ever scored.
cfg.jawFeature  = 'jaw_ydisp_view1';
cfg.normPctLow  = 2;   % Methods: 2nd percentile, as the tongue scripts
cfg.normPctHigh = 99;   % Methods: 99th percentile

% the biggest jaw openings, which are exactly the peaks the decoder is judged on.
% Percentile SCALING is kept, so the units are unchanged; values just run
% slightly outside [0,1]. true reproduces the old behaviour.
cfg.clipNormalizedRange = false;
cfg.reportDropout       = true;   % print the % of untracked jaw samples per session

%% RIDGE
cfg.ridgeGrid       = logspace(-2, 6, 50);
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
cfg.cvSeed      = 7;   % base seed; each pass offsets it by repIdx
RUN.cvSeed0     = cfg.cvSeed;   % remembered so the offset is from a fixed base            % [] = folds in trial order, no shuffle
cfg.lambdaFixed = [];   % [] = GCV. Set a number to pin lambda instead.

%% DECODING INDEX: what it is measured against
% THE BASELINE IS ZERO. The null model is "predict a scaled jaw of 0 and never
% move", the same null the tongue scripts use, so a jaw index and a tongue index
% from the matching session are the same quantity and can be put side by side.
% There is no rest-baseline option in this script and no per-trial reference to
% compute -- the denominator is fixed by construction.
% WHAT 0 IS HERE, so the number is read correctly. The target is scaled as
% (jaw - p3)/(p99 - p3) with the percentiles taken over every tracked sample in
% the session (see the TARGET block). Scaled 0 is therefore the 3rd percentile of
% that distribution. Most of a session is the animal not licking, so the bulk of
% the distribution is resting jaw and the 3rd percentile sits somewhat BELOW the
% resting position, in the low tail of resting jitter and tracker noise.
% The consequence is worth stating plainly rather than hiding in a comment: the
% resting jaw is at some scaled value r > 0, so predicting 0 is already wrong by
% about r even while the animal is still. RMSE_0 is therefore LARGER than the
% error of predicting rest, the ratio is smaller, and the index runs HIGHER than
% a rest-referenced one would. How much higher depends on where the tracker's
% reference happened to sit in that session, which is a rig property, not a
% neural one -- so treat absolute jaw index values as calibrated only within a
% session, and lean on the WITHIN-session comparison across contacts (contact L
% vs contact 1), which is what the delta figures show and what the disengagement
% question actually asks.
% The one confound that survives either baseline: a constant baseline shrinks
% when the excursion shrinks, so a purely behavioural decline in jaw opening
% across the bout lowers the index on its own. peakL, the per-contact excursion,
% is carried alongside the index for exactly this reason -- check it before
% reading a falling index as cortical.

cfg.nLicksAnalyze = 6;
% ANALYSIS WINDOW the per-contact index is computed over (declared with the other
% windows above; repeated here because the segmentation runs on exactly this
% span, and the persistence threshold is a percentile of the trace INSIDE it).

%% JAW CYCLE SEGMENTATION
% this file is the same 937 lines, comments included. It is the jaw analogue of
% the tongue scripts' above-zero run detection, and it has to be a different
% algorithm because jaw displacement has no resting floor to threshold.
% In one paragraph: build the complete alternating max/min sequence of the
% smoothed trace, then repeatedly delete the least prominent adjacent extremum
% PAIR until everything left clears cyclePersistFrac x the trial's robust jaw
% pair separated by a tiny amplitude, so it is merged away and can never become
% a cycle boundary; a real lick cycle has persistence on the order of the jaw's
% working range and always survives. Surviving minima ARE the boundaries, and
% between consecutive minima there is exactly one peak by construction. Cycles
% are then matched to lick contacts CONTAINMENT FIRST, with a nearest-peak pass
% genuinely smaller than its neighbours.
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
% Only read by the legacy 'peakDetect' anchor, kept so it still runs:
cfg.jawLickIndexRule   = 'contactMatched';   % 'contactMatched' | 'peakOrder'
cfg.peakMinHeight      = 0.08;
cfg.peakProminence     = 0.05;
cfg.peakMinDist        = 5;
cfg.troughProminence   = 0.0001;
cfg.troughMinDist      = 7;
cfg.troughMaxValue     = 0.25;
cfg.segHalfWin         = 5;
cfg.requireTroughBelowMax = false;
% localJawSegments reads this name; it is the same number as cfg.nLicksAnalyze
% and is set from it so the two can never drift apart.
cfg.numLicksToAnalyze = cfg.nLicksAnalyze;
cfg.reportSegmentation = true;   % print how many contacts got a cycle, per session

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
cfg.requireBout   = false;   % <-- the toggle
cfg.boutMinLicks  = 6;
cfg.boutWin_s     = 1.50;
cfg.boutMaxILI_s  = 0.25;

%% STATISTICS
cfg.earlyLicks = 1:2;
cfg.lateLicks  = 3:6;
cfg.alpha      = 0.05;

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
cfg.nDiagSessions    = 0;

% The WARPED figure can be decoupled from the six-panel one. It is the panel that
% shows whether the CYCLE SEGMENTATION actually found the licks, which is worth
% having for every session, whereas 26 six-panel tabs is a lot to scroll through.
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

% HOW THE WARPED FIGURE IS BUILT.
%                 map per inter-lick interval takes this trial's contact times
%                 onto the across-trial median template (uniform winBoutGap-free
%                 spacing of warpDelta_s), and the trace is resampled through the
%                 INVERSE of that map. The output is ONE CONTINUOUS TRACE on a
%                 warped time axis, with a dotted line at each template contact
%                 -- no gaps, because nothing is cut into slots.
%   'slots'       the tongue scripts' layout: each cycle resampled into its own
%                 fixed-width slot with NaN between them. Kept so the two can be
%                 compared, but it is not what the jaw figures looked like.
% Continuous is right for jaw: jaw displacement never leaves the trace, so
% cutting it into per-cycle slots throws away the inter-lick trace, which is
% exactly where the jaw sits between licks and is worth seeing.
cfg.warpStyle    = 'continuous';   % 'continuous' | 'slots'
cfg.warpDelta_s  = 0.15;   % fixed inter-lick spacing of the template

% HOW FAR THE WARP REACHES vs HOW MUCH IS SHOWN -- two different numbers.
%   cfg.warpNumLicks   how many contacts the TEMPLATE spans, i.e. how much of the
%                      bout actually gets warped. Set it generously: any part of
%                      the trace past the last mapped interval is handled by
%                      cfg.warpEdgeMode rather than by a real landmark, so a
%                      short template leaves a long extrapolated tail. 12 covers
%                      essentially every bout in this data.
%   spec.warpShowLicks how many contacts are LABELLED and kept inside the x
%                      limits. Per study, because the tasks differ: 8 for the
%                      R-series and learning, 6 for VTA, whose test trials are
%                      selected on having >= 4 contacts and often have few more.
cfg.warpNumLicks = 12;

% What happens OUTSIDE the mapped intervals -- before the first and after the
% original) samples those regions on the RAW clock while everything between
% landmarks is on the WARPED clock, and the seam between them is the sharp
% step that appears just before contact 1 in both the actual and the predicted
% trace. 'extend' carries the nearest affine map outwards so the trace is
cfg.warpEdgeMode = 'extend';   % 'extend' | 'nan' | 'identity'

% With TWO test sets (R1 vs R4, R1 vs R6), draw both on the SAME warped axes --
% solid actual, dashed predicted, each in that test set's colour -- instead of
% pooling them into one pair of traces. Ignored where there is only one test set.
cfg.warpSplitTestSets = true;
cfg.warpedTickCentre = false;

% THE DENSER LAYOUT (what the first tabbed version drew, before the port): a 25-
% sample segment with a 1-sample gap, linear interpolation, centred labels. Paste
% these five lines over the ones above to get it back.

%% PARAMS (field names read by the loaders)
params.alignEvent = 'firstLick';
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 8;
params.lowFR  = 0.01;   % minimum mean firing rate, Hz
params.quality    = {'good'};   % good units only (findClusters trims blanks and ignores case)
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
%% STUDY: VTA cohort
spec.name = 'VTA';
spec.sessionDates = { ...
    '2025-02-15','2025-02-17','2025-02-18','2025-02-19', ...   % TDv1
    '2025-02-25','2025-02-26','2025-02-27','2025-02-28', ...   % TDv4
    '2025-08-21','2025-08-22','2025-08-23','2025-08-24','2025-08-25','2025-08-26', ...   % TDv6
    '2025-08-25','2025-08-26','2025-08-28','2025-08-30','2025-08-31' };   % TDv5
spec.sessionLoaders = { ...
    @loadTDv1_neur22, @loadTDv1_neur, @loadTDv1_neur, @loadTDv1_neur, ...
    @loadTDv4_neur,   @loadTDv4_neur, @loadTDv4_neur, @loadTDv4_neur, ...
    @loadTDv6_neur,   @loadTDv6_neur, @loadTDv6_neur, @loadTDv6_neur, @loadTDv6_neur, @loadTDv6_neur, ...
    @loadTDv5_neur,   @loadTDv5_neur, @loadTDv5_neur, @loadTDv5_neur, @loadTDv5_neur };
spec.groupMaps   = { [1 1 1 1  1 1 1 1  2 2 2 1 2 2  2 2 0 1 2] };   % M1
spec.groupLabels = {'M1'};
%         (fewer than lickCountMin contacts in the first lickCountWin_s), topped
%         up with long-bout trials only when the short ones run out.
spec.poolCondIdx    = 1;
% first and topped up with long-bout trials only if they do not fill the quota.
% 'lickCount'     = the old hard partition (all short train, all long test), which
%                   is trainFrac-free and completely deterministic.
spec.splitRule      = 'lickCountFrac';   % 'lickCountFrac' | 'lickCount'
spec.trainFrac      = 0.70;   % Methods: ~70% for the VTA Reward task
spec.lickCountMin   = 4;
spec.lickCountWin_s = 1.25;
spec.testSets       = struct('name','manyLick','condIdx',{[]});
spec.trialCaps           = { 'TDv1','2025-02-15', 211 };
spec.singleProbeSessions = { 'TDv1','2025-02-15' };

spec.minTestTrials    = 10;   % 'lickCountFrac': never leave fewer long-bout
% trials than this to score on
spec.testFromManyOnly = true;   % 'lickCountFrac': test on long-bout trials only,
% as the plain 'lickCount' rule does

spec.warpShowLicks    = 6;   % contacts LABELLED on the warped figure

%% bout filter + diagnostics
% BOUT FILTER OFF HERE, and it has to be. This script's split rule IS a lick
% count -- fewer than 4 contacts trains, 4 or more tests -- so demanding 8
% contacts in a row would delete the training set outright. cfg.requireBout is
% ignored in this script.
spec.applyBoutFilter = false;
spec.diagRule        = 'random';   % 'random' = cfg.nDiagSessions sessions at random
spec.diagGroups      = [];   % used only by diagRule 'groups'

%% colours
% One test set, GREEN throughout.
spec.tsColours  = [0.10 0.60 0.25];   % manyLick
spec.grpColours = [];
spec.predColour = [0.10 0.60 0.25];   % prediction on the diagnostics

%% figure behaviour
spec.overlayGroups = false;   % false = never draw two groups on one axes
spec.showFigCI     = true ;   % Figure 1: decoding index with error bars
spec.showFigDelta  = false;   % Figure 3: delta from contact 1
spec.deltaGroups   = [];   % groups on the delta figures ([] = all)
spec.figRef       = 'Fig. S4H';   % the manuscript panel this file produces
spec.vsC1Tail     = 'left';   % 'both' | 'left' | 'right', from the legend
spec.vsC1Test     = true ;   % per-contact test vs contact 1
spec.betweenTest  = false;   % per-contact test BETWEEN the two trial types
spec.groupTest    = 'none';   % 'none' | 'ranksum' | 'crossTask'
spec.groupSummaryLicks = 3:6;   % contacts averaged into the single summary value
spec.summaryFile  = fullfile(tempdir, 'decodeSummary_jaw_vta_clean.mat');
spec.crossTaskFile = fullfile(tempdir, 'decodeSummary_jaw_r1_clean.mat');
spec.compareGroups = [];   % [a b] = per-contact paired test between two groups
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
fprintf('  decoding index baseline: ZERO (1 - RMSE_dec/RMSE_0, same null as the tongue scripts)\n');
fprintf('  jaw cycles: %s segmentation, ported unchanged from the jaw scripts\n', cfg.segAnchor);
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
assert(size(S.jawRaw,1) == nT, 'sess %d: kinematics %d time samples, time %d.', ...
    sessionIdx, size(S.jawRaw,1), nT);

% ---------------- TARGET: JAW DISPLACEMENT ----------------
% Percentile-scaled, NOT clipped, and NaN LEFT AS NaN. A NaN here is the tracker
% losing the jaw marker, not the jaw being anywhere in particular -- see the long
% note in SETTINGS. Nothing downstream fills it, and every score uses a finite
% mask shared by actual and predicted.
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
% and its output is then used for three different things: the per-trial fitting
% warped figure. Running it once is not just a speed choice -- it is what makes
% "contact 1" mean exactly one thing in this script.
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
        fprintf('          %d trial(s) where cycle order and contact order disagree (see localJawSegments)\n', nSegDisagree);
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
% The ORIGINAL rule, kept: a hard partition on contact count. Every
% short-bout trial trains, every long-bout trial tests, nothing is
% sampled -- which is why this rule gives the same answer every run
% whatever cfg.rngSeed is set to.
        testTrials  = poolTrials(lickCount(poolTrials) >= spec.lickCountMin);
        trainTrials = poolTrials(lickCount(poolTrials) <  spec.lickCountMin);

    case 'lickCountFrac'
% pool and fills it FEW-LICK TRIALS FIRST: all the short-bout trials go
% in, and only if they do not fill the quota is it topped up with
% long-bout trials drawn at random. 100 pool trials at trainFrac 0.7 with
% 40 short-bout trials -> all 40 short + 30 long train, 30 long test.
% TWO GUARDS, both of which matter:
%   spec.minTestTrials    the top-up is capped so at least this many
%                         long-bout trials are LEFT to score. Without it
%                         a high trainFrac on a session with few short
%                         trials eats the test set entirely.
%   spec.testFromManyOnly true keeps the test set to long-bout trials, as
%                         the original rule did, so the per-contact n at
%                         contacts 4-8 is not diluted by trials that
%                         never had a contact 4. Leftover short-bout
%                         trials are simply unused. false puts every
%                         non-training trial in the test set instead.
% READ THIS BEFORE COMPARING TO THE 'lickCount' NUMBERS. That rule is a
% pure DOMAIN SHIFT: the model never sees a long bout. As soon as the
% top-up is non-zero this becomes a blend, and at a trainFrac high enough
% to exhaust the short trials it is mostly a within-distribution split.
% The decoding index will rise for that reason alone. The split is also
% RANDOM now, so runs differ unless cfg.rngSeed is set.
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
% resampled to cfg.warpedSegmentLength points and dropped into its own slot, mean
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
% INVERSE of that map. The result is ONE CONTINUOUS TRACE on a warped time
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
% condition: actual black, predicted in that condition's colour, labelled
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
    if cfg.warpedTickCentre, tickAt = slotStart + cfg.warpedSegmentLength/2; else, tickAt = slotStart; end
    xticks(ax, tickAt(1:cfg.maxLickShow));
    xticklabels(ax, arrayfun(@(L) sprintf('contact %d', L), 1:cfg.maxLickShow, 'UniformOutput', false));
    xtickangle(ax, 45);
    xlim(ax, [slotStart(1)-5, slotStart(cfg.maxLickShow)+cfg.warpedSegmentLength+5]);
    legend(ax, [hA hP], {'actual','predicted'}, 'Location','northeast','FontSize',8);
    title(ax, {hdr, sprintf('warped jaw cycles 1-%d, TEST | mean \\pm 1.96 SEM | n = %d trials', ...
        cfg.maxLickShow, numel(testTrials))}, 'FontSize',10,'FontWeight','normal');
    box(ax,'off')
end

end   % wantWarp

clear spikes Xtrain Xtest Ptr Pte predFull Mact Mpred Wa Wp
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
  if ~spec.vsC1Test, break; end
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
%    between the Simple and Double Reward Tasks, and over the final two licks of
%    the analysis window for the comparison between the first and fifth day of
%    training."
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
% only the samples inside ITS OWN window (variable length), so the trace is
% averaged over however many trials cover each sample; the count is shown on
% the right axis and thins out at the edges.
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
            ylabel(ax, 'jaw displacement (norm)', 'FontSize', 11);   % NO [0 1] clamp: jaw is
% percentile-scaled, not bounded
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
    finishAxes(ax, nLA, '1 - RMSE_{dec}/RMSE_{0}', h, lbl, [0 1]);
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
fprintf('\n%s\n', repmat('=',1,110));
fprintf('%-5s %-8s %-12s %5s %-8s %6s %6s %6s %10s %7s %9s %9s %11s\n', ...
    'sess','animal','date','probe','region','units','nTrn','nTst','lambda','edf', ...
    'R2 train','R2 test','win samp');
fprintf('%s\n', repmat('-',1,110));
for r = res
    wl = r.winLenSamples;
    fprintf('%-5d %-8s %-12s %5d %-8s %6d %6d %6d %10.3g %7.1f %9.3f %9.3f %5d-%-5d\n', ...
        r.session, r.anm, r.date, r.probe, r.region, r.nUnits, r.nTrain, r.nTest, ...
        r.lambda, r.edf, r.r2TrSel, r.r2TeSel, min(wl), max(wl));
end
fprintf('%s\n', repmat('=',1,110));

end   % repIdx -- REPEAT PASSES

%% LOCAL FUNCTIONS (shared)

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
%                   winBoutEndPad_s. Contact times only -- it never looks at the
%                   jaw trace, which is deliberate: the point of the mode is to
%                   fit the late contacts, and a kinematic rule would smuggle the
%                   same cycle-shaped assumption back in. It also makes this
%                   window identical to the tongue script's.
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
% Read as a stopping rule, not as "the last lick in the trial". Every trial's
% final lick is trivially not followed by another, so the literal reading is
% vacuous and would sweep in a stray lick seconds after the bout ended. What the
% gap identifies is the end of the CONTINUOUS train.
    d = diff(contactRel(:));
    k = find(d > gap_s, 1);
    if isempty(k), k = numel(contactRel); end
end

function [decIdx, peakL, nTrialL] = jawDecodingIndex(Yact, Ypred, segs, cfg)
% samples, pooled across trials.
% RMSE_0 is the error made by predicting a scaled jaw of 0 and never moving --
% the same null the tongue scripts use. It is a CONSTANT, identical for every
% trial and every session, so nothing about the denominator varies with the
% animal's posture or with how the bout ended.
% It is still true that a fixed baseline shrinks when the excursion shrinks: a
% purely behavioural decline in jaw opening across the bout lowers the index on
% its own, with no change in how well the neural data predicts. peakL, the
% per-contact peak above 0, is returned alongside so that can be checked.
% ONE FINITE MASK, shared by actual and predicted, for BOTH RMSEs. This is the
% asymmetry that has to be avoided with a NaN-bearing target: mask the numerator
% at dropout samples but not the denominator and the index rises with the
% tracking failure rate.
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

function h = drawWarpedSet(ax, W, col, slotStart, cfg)
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
% colour, and the two labelled with TEXT ON THE AXES rather than a legend box --
% a key costs a corner of the plot and an eye movement, and with one condition
% per panel there is nothing to disambiguate.
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
%   S.time, S.trialdat, S.nUnitsTotal, S.jawRaw, S.cluid, S.trialid,
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
% Jaw lick segmentation.
% THE BUG THIS ROUND FIXES: TROUGHS LANDING ON PLATEAU WIGGLES
% When the jaw reaches maximum opening it frequently hangs there and jitters
% before closing. Every one of those jitters is a legitimate local minimum, and
% it sits CLOSER to the peak than the real inter-lick trough does. The previous
% proximity won and the segment boundary landed on the wiggle at the top of the
% excursion instead of at the bottom of the cycle. Licks whose excursion
% happened to have no plateau were unaffected, which is why some licks in a
% trial looked right and others did not.
% Two things are wrong with the old approach and both are fixed by inverting the
% anchor. Peak detection is a LOCAL criterion applied to a noisy trace, and lick
% contacts are ground truth we already have. So:
%                                    see the persistence branch below)
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
%      climb into the next cycle always does. Crucially, NOTHING DEPENDS ON A
%      WINDOW BEING THE RIGHT WIDTH -- see the long comment at the walk itself
%      licks and on licks with a wider-than-median gap.
%   3. BROAD BASINS. If the bottom is a flat stretch rather than a sharp V, the
%      boundary moves to the end of the CONTIGUOUS run of samples within
%      cfg.troughFlatTol of the trough -- the moment opening begins on the way
%      in, bottom reached on the way out. Contiguity matters: the earlier rule
%      accepted any sample in the window within tolerance, so a shoulder further
%      up the flank could cut the segment short. Default tolerance is now 0,
%      i.e. take the actual trough; wiggle protection comes from the walk.
%   4. NUMBERING IS EXACT BY CONSTRUCTION. Segment k IS contact k. There is no
%      peak-to-contact matching step left, so nothing can skip a lick or slide a
%      label -- the earlier "band 2 over contact 3, then band 4" failure is
%      structurally impossible now rather than merely fixed.
%   5. AMPLITUDE GATE. peak minus the deeper of its two troughs must reach
%      cfg.minCycleAmp, so a contact with no real jaw excursion (dropout, or a
%      lick made with almost no opening) does not get a window drawn around
%      noise. Set cfg.minCycleAmp = 0 to segment every contact regardless.
% findpeaks still runs, but under 'contactCycle' it only supplies a count for
% the audit; it decides no boundary. Set cfg.segAnchor = 'peakDetect' to get the
% honours cfg.jawLickIndexRule, cfg.requireTroughBelowMax and cfg.segMode).
% ALSO STILL TRUE FROM THE PREVIOUS ROUND: segments are one full cycle
% (trough -> peak -> trough), so their LENGTH VARIES between licks and trials.
% NaN HANDLING is unchanged: findpeaks and min/max cannot be trusted across
% NaN, so detection runs on a gap-filled COPY and the returned indices are
% applied to the real NaN-bearing trace by the caller. No interpolated sample is
% ever scored.
% OUTPUTS
%   nDisagree  trials where contact-anchored windows differ from what plain
%              peak-order would have produced (a scale for how much the anchor
%              matters on this data)
%   segDiag    nTrials x 4: [nContacts nSegmented nPeaksFindpeaks nFailedAmp]

% ---- defaults, so this drops in without touching the cfg block ----
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
% Every previous version tried to find each boundary by a LOCAL test --
% "is the nearest labelled trough here", "has the trace risen enough
% above the running minimum", "has it descended far enough first". Each
% trace, because a local test cannot tell a plateau wiggle from a real
% trough: the two have the same SHAPE. They differ only in SCALE, and
% scale is a global property of the waveform, not something visible in a
% three-sample neighbourhood.
% So this version stops asking local questions. It builds the full
% alternating max/min sequence of the trace and then repeatedly deletes
% the least prominent adjacent extremum PAIR until every remaining
% extremum has prominence at least persistThresh. This is topological
% persistence simplification, and it is exactly the right tool here:
%   - a plateau wiggle is a max/min pair separated by a tiny amplitude.
%     It has near-zero persistence and is merged away early, so it can
%     NEVER become a boundary. Not "usually not" -- never, for any
%     wiggle below threshold, regardless of where it sits on the cycle
%     or how it is shaped.
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
% This is the lick-1 fix. Persistence normally takes the SMALLER
% of the two drops on either side of an extremum, which is right
% when both neighbours are real troughs. The first element of E
% is not a real trough -- it is wherever the analysis window
% happened to be cut (t = -0.09). When lick 1 has no trough
% before it inside the window, the jaw is already part-way open
% at that cut, so the "drop" on the left is small, and lick 1's
% arbitrary number.
% So at a window edge, persistence is measured on the side that
% has a real trough, and only that side.
                if E(1) == 1 && ~isempty(pr)
                    pr(1) = d2(1);   % left neighbour is the window start
                end
                if E(end) == nSamp && ~isempty(pr)
                    pr(end) = d1(end);   % right neighbour is the window end
                end

                [mn, jj] = min(pr);
                if mn >= persistThresh, break; end
                j = jj + 1;   % index into E

% Remove the low-persistence extremum together with the
% neighbour it merges into. Removing a PAIR is what keeps the
% sequence alternating; removing one alone would leave two
% maxima adjacent and corrupt every later step.
                takeLeft = (d1(jj) <= d2(jj));
% ...but never delete a window edge. It is the boundary that
% lets an incomplete first or last cycle still have a window,
% and merging it away is how lick 1 disappeared.
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
% sitting exactly ON the left window edge matched bnd = 1 on
% BOTH sides; the "peak is on an edge" branch then pushed the
% closing boundary out to nSamp, and lick 1's cycle became the
% entire 5-second window. The edge fallback has to know WHICH
% edge the peak is on, and strict comparisons plus a per-side
% default give that for free.
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
% and the last cycle in trial 149.
% A lick cannot last much longer than the animal's own lick
% period, so that is the bound: at most
% cfg.cycleMaxHalfILIFrac x the trial's median inter-lick
% interval on each side of the peak. When the real trough is
% further out than that, the boundary becomes the LOWEST point
% within the allowance rather than a flat cutoff, so on a slow
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

% amplitude check, ignoring any side that is a window edge rather
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
% THE PROBLEM THIS SOLVES. cfg.cyclePersistFrac is a GLOBAL,
% per-trial threshold: an excursion has to clear a fraction of the
% WHOLE trial's jaw range to survive simplification. That is what
% makes it robust against plateau wiggles, and it is also its one
% blind spot -- a lick whose opening is genuinely smaller than its
% neighbours in the same trial falls under the bar and is merged
% 5 missing while 2 and 4 are correct, and lick 1 missing because
% "the peak isn't as high as the rest".
% Raising or lowering the global threshold cannot fix that. Lower it
% and the wiggles come back everywhere; raise it and more small licks
% vanish. The information that resolves it is not in the amplitude at
% has no cycle but contacts k-1 and k+1 both do, then the stretch of
% trace between k-1's closing trough and k+1's opening trough is,
% by construction, exactly one inter-lick interval with a contact
% sitting inside it. Whatever excursion is in there IS lick k's,
% however small, because there is nowhere else for it to be.
% So the span is bounded by two boundaries already trusted, and the
% rescue only has to find the peak inside it. No global threshold is
% involved, which is why this recovers a small lick without
% loosening anything for the rest of the trial.
% GUARDS, so this cannot invent cycles out of flat trace:
%   - the span must be no wider than cfg.rescueMaxSpanILIFrac x the
%     trial's own inter-lick interval. A wide span means a contact is
%     missing on one side too and the bracket is not one cycle.
%   - the contact itself must fall inside the span.
%   - the peak must clear BOTH ends by cfg.rescueMinAmpFrac of the
%     trial's jaw range. A flat stretch has no peak and is refused.
%   - the peak must be strictly interior; a maximum sitting on a span
%     edge means the excursion belongs to the neighbour, not here.
% EDGE CONTACTS (cfg.rescueEdgeContacts). Lick 1 has no left
% neighbour and the last lick has no right one, so the window edge
% stands in for the missing boundary. That is the image-3 case: lick
% 1's opening is off-screen or shallow, and the span from the window
% start to lick 2's opening trough contains exactly one excursion.
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

% Usually that lands on the span edge itself, which is the
% right answer: the previous cycle's close IS this cycle's
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
% contact, the last lick has no next one) -- the trough search
% backstop is now taken from the LOCAL gap, not from this median.
% Using the median as a cap everywhere is what truncated the search
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

% half-distance to each neighbouring contact bounds the window;
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
% Replaces "argmin of a capped window", which failed in two places
% notion of "have I actually left the trough yet", so it just returns
% whatever the lowest sample in the window happens to be:
%   EDGE LICKS (lick 1's opening trough, the last lick's closing
%   trough). There is no neighbouring contact peak to bound the
%   search, so the window ran out to a fixed cap or the trace edge
%   and argmin returned a point on the pre-lick baseline or the
%   post-lick settle, not the trough of that cycle.
%   where c4->c5 and c6->c7 are far wider than the rest). The cap
%   THE TROUGH and argmin returned the window edge -- a boundary
%   sitting on the flank, mid-descent.
% The walk instead asks the right question directly: step outward
% from the peak tracking the running minimum, and stop as soon as
% the trace has climbed back up by more than a threshold above that
% minimum -- meaning we are demonstrably out of the trough and into
% the next rise. The running minimum at that moment IS the trough.
% The threshold is cfg.troughRiseFrac of the excursion depth so far,
% with cfg.troughRiseFloor as an absolute noise floor, so:
%   - a plateau wiggle near the peak cannot stop the walk (depth is
%     small there, so the floor governs and the wiggle is under it),
%   - a real climb out of the trough always does (depth is large by
%     then, so a genuine rise clears the fraction),
%   - and nothing depends on a window being the right width, which
%     is what broke at the edges and on the long gaps.
            segStart = nan(nC,1); segEnd = nan(nC,1);
            for k = 1:nC
                pk = pkIdx(k);
                if ~isfinite(pk), continue; end

% hard bounds: the neighbouring CONTACT peaks, so a cycle can
% never cross into an adjacent lick's excursion
                prevPk = 1;
                kk = find(isfinite(pkIdx(1:k-1)), 1, 'last');
                if ~isempty(kk), prevPk = pkIdx(kk); end
                nextPk = nSamp;
                kk = find(isfinite(pkIdx(k+1:end)), 1, 'first');
                if ~isempty(kk), nextPk = pkIdx(k+kk); end

% Backstop, now from the LOCAL gap rather than the trial median.
% This is the second half of the long-gap fix: a lick with a wide
% gap on one side gets a correspondingly wide allowance on that
% side instead of being held to the typical interval. It is only
% a backstop -- the walk normally stops well inside it.
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
% maximum opening it can hang high and jitter for tens of ms
% before closing. A jitter is a rise above the running minimum
% terminating up on the plateau and calling a HIGH point the
% trough -- "mis-identifying the times when the peak stays high
% and wiggles as the trough".
% Rises alone cannot distinguish the two, because the shapes are
% the same; what separates them is HOW FAR DOWN the trace has
% come. So: look ahead over the whole permitted interval, see how
% deep this side can possibly go, and refuse to let the walk stop
% until it has covered cfg.troughMinDepthFrac of that. A plateau
% wiggle sits at a small fraction of the full depth and is now
% simply ignored; a rise near the bottom still stops the walk
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
% That means the trough is off-screen or past the neighbouring
% peak, not that we found one.
                truncBack = (lo <= lo0);
% contiguous flat basin only: if the bottom is a flat stretch
% rather than a sharp V, start at the end of THAT stretch, the
% moment opening begins. Contiguity matters -- the old rule took
% any sample in the window within tolerance of the minimum,
% which let an unrelated shoulder further up the flank win.
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

% --- amplitude gate, now SIDE-AWARE ---
% This is the other half of "lick 1 is missing". The gate used
% the SHALLOWER of the two troughs (max of the two trough
% values). For lick 1 the opening trough is frequently off the
% left edge of the analysis window -- the jaw is already open at
% that reads "1 below amplitude gate" is missing lick 1 for this
% So a boundary that hit the window edge without turning is not
% evidence about amplitude, and is excluded from the test. Only
% unmeasurable, and then the shallower side stands.
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

% ---- optional legacy window shape, same anchors ----
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
% with the NEAREST trough on either side, then match peaks to contacts.
% This is the path whose nearest-trough rule put boundaries on plateau
% wiggles; it is here so the two can be plotted side by side, not
% because it should be used.
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
% ADAPTED IN ONE PLACE ONLY: this loader hands the decoder a slim struct, so the
% goCue and lickL arrays are passed in directly instead of sessionObj. Every
% REAL lick-contact times for this trial (not jaw-oscillation peaks).
% sessionObj.bp.ev.goCue/lickL are on each trial's own ABSOLUTE clock, while
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
% timeWarp2's findMedianJawTimes, applied to REAL contact times: per lick index,
% the median across trials that have that many licks, then timeWarp2's fixed
% uniform-spacing override. licksRel(1) is always exactly 0, so the template
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
% timeWarp2's trialWarpFits_Jaw on real contact times: one affine map per
% inter-lick interval (goCue->lick1, lick1->lick2, ..., last lick->end), taking
% this trial's boundaries onto the across-trial median template. mode(goCue) in
% the original is replaced with nanmedian(goCue) -- goCue times are near-
% Interior boundaries use the EXACT shared anchor (no "-dt" back-off): the
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
% violate. Grid points outside every mapped interval keep the unwarped
% (identity) value -- there is no landmark to warp against there.
    ls = contacts;
    mls = medStart;
    nL = numel(ls);
% it fixes the sharp step just before contact 1.
%              UNWARPED value. The first interval starts at the median go cue, so
%              everything before that is sampled on the raw clock while
%              everything after is sampled on the warped clock -- and the two
%              disagree, which puts a discontinuity exactly at the go cue. That
%              is the small bump that rises and then drops off a cliff. It is an
%              artifact of the seam, not jaw movement; it appears in the actual
%              AND the predicted trace at the same x, which is the giveaway,
%              since no physiological event could do that to both.
%   'extend'   the FIRST interval's affine map is extrapolated
%              backwards to cover everything before it, and the LAST interval's
%              forwards. Same map on both sides of the seam, so the trace is
%              continuous by construction and the pre-cue baseline is still real
%              data, just warped by the nearest available map.
%   'nan'      leave the unmapped edges empty. The most conservative option:
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
% Mean +/- 95% CI band. Replaces the shadedErrorBar dependency the original
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

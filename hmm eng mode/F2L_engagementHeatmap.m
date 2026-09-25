%% F2L_engagementHeatmap.m
%  Engaged and disengaged motor cortical states, Delayed Reward Task.
%  Per-trial disengagement times come from the HMM-GLM fits in the
%  Disengagement Times folder. Neural activity is aligned to those times to
%  define the engagement mode and measure how fast the state switches.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%    Disengagement Times\<session>\   HMM-GLM state fits
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/100
%    params.smooth      15
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


%% CONFIG (edit these)
% ---- the unit test ----
% 'signrank' = paired two-tailed Wilcoxon signed-rank   <- Methods, and current
% 'ttest1t'  = one-tailed t-test in BOTH directions     <- what the Fig. 2L legend says
cfg.unitTest  = 'signrank';
cfg.unitAlpha = 0.005;
cfg.unitFDR   = false;   % true = Benjamini-Hochberg across units instead of a
% fixed alpha. The run already prints the expected
% false-positive count; this is the defensible version.

% ---- engagement mode ----
% Half-width of the pre and post windows around the transition. The windows are
% [-W, 0) and [0, +W), derived from the real time axis, equal in length and
% non-overlapping. Methods say 0.400; every one of these scripts uses 0.300.
cfg.engModeWin_s = 0.400;   % +/-400 ms around the transition

% Methods: p(t) = sum_i w_i r_i(t) with sum|w_i| = 1. 'mean' divides that by the
% unit count, which rescales every session by its own N before sessions are
% averaged together.
cfg.projMode = 'sum';   % 'sum' (Methods) | 'mean' (previous behaviour)

% ---- figures ----
cfg.xlimPlot    = [-1.25 1];   % time axis on every heatmap
cfg.sortWin_s   = 0.5;   % peak-rate sort window either side of t = 0
cfg.climFixed   = [0 150];   % fixed colour scale (Hz)
cfg.climPctile  = 99.7;   % percentile scale, for the first figure

% ---- paths ----
% HMM results now live in <MATLAB Codes _ v2>\Data\HMM: cfg.hmmRoot  = 'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final';

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
cfg.hmmRoot = fullfile(v2Root, 'Disengagement Times');   % HMM-GLM disengagement times, one folder per session and probe
assert(exist(cfg.hmmRoot, 'dir') == 7, 'No disengagement-time folder: %s', cfg.hmmRoot);

datapth = '';   % raw data folder not used (was: datapth = fullfile(cfg.codeRoot,'data');)

%% PARAMS
params.alignEvent = 'goCue';
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz
params.tmin       = -2;
params.tmax       = 5;
params.dt         = 1/100;
params.smooth     = 15;
params.quality    = {'good'};   % good units only (findClusters trims blanks and ignores case)

params.condition = { ...
    'hit==1 | hit==0', ...
    'hit==1 & trialTypes == 1 & rewardedLick == 1', ...
    'hit==1 & trialTypes == 2 & rewardedLick == 1', ...
    'hit==1 & trialTypes == 3 & rewardedLick == 1', ...
    'hit==1 & trialTypes == 1 & rewardedLick == 4', ...
    'hit==1 & trialTypes == 2 & rewardedLick == 4', ...
    'hit==1 & trialTypes == 3 & rewardedLick == 4', ...
    'hit==1 & rewardedLick == 1', ...
    'hit==1 & rewardedLick == 4', ...
    'hit==1' };

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

%% SESSIONS
% {loader, date, HMM folder for probe 1, HMM folder for probe 2}
% Where the two folders are the same, only one HMM result exists for that
% session and both regions are aligned with it -- see the consistency check.
sessionTable = { ...
 @loadTD1_neural , '2023-02-21', 'TD1d_2023_02_21_P1' , 'TD1d_2023_02_21_P2'  ; ...
 @loadTD1_neural , '2023-02-22', 'TD1d_2023_02_22_P1' , 'TD1d_2023_02_22_P2'  ; ...
 @loadTD1_neural , '2023-02-23', 'TD1d_2023_02_23_P1' , 'TD1d_2023_02_23_P2'  ; ...
 @loadTD1_neural , '2023-02-24', 'TD1d_2023_02_24_P1' , 'TD1d_2023_02_24_P2'  ; ...
 @loadTD4_neural , '2023-02-21', 'TD4d_2023_02_21_P2' , 'TD4d_2023_02_21_P2'  ; ...
 @loadTD4_neural , '2023-02-24', 'TD4d_2023_02_24_P2' , 'TD4d_2023_02_24_P2'  ; ...
 @loadTD4_neural , '2023-02-25', 'TD4d_2023_02_25_P2' , 'TD4d_2023_02_25_P2'  ; ...
 @loadTD4_neural , '2023-03-19', 'TD4d_2023_03_19_P1' , 'TD4d_2023_03_19_P2'  ; ...
 @loadTD13_neural, '2024-11-12', 'TD13d_2024_11_12_P1', 'TD13d_2024_11_12_P1' ; ...
 @loadTD13_neural, '2024-11-13', 'TD13d_2024_11_13_P2', 'TD13d_2024_11_13_P2' ; ...
 @loadTD13_neural, '2024-11-21', 'TD13d_2024_11_21_P2', 'TD13d_2024_11_21_P2' ; ...
 @loadTD15_neural, '2024-11-24', 'TD15d_2024_11_24_P2', 'TD15d_2024_11_24_P2' ; ...
 @loadTD15_neural, '2024-11-25', 'TD15d_2024_11_25_P1', 'TD15d_2024_11_25_P1' ; ...
 @loadTD15_neural, '2024-11-26', 'TD15d_2024_11_26_P1', 'TD15d_2024_11_26_P1' ; ...
 @loadTD15_neural, '2024-11-27', 'TD15d_2024_11_27_P2', 'TD15d_2024_11_27_P2' ; ...
 @loadTD22_neural, '2025-06-17', 'TD22d_2025_06_17_P2', 'TD22d_2025_06_17_P2' ; ...
 @loadTD22_neural, '2025-06-18', 'TD22d_2025_06_18_P1', 'TD22d_2025_06_18_P2' ; ...
 @loadTD22_neural, '2025-06-19', 'TD22d_2025_06_19_P1', 'TD22d_2025_06_19_P2' ; ...
 @loadTD22_neural, '2025-06-20', 'TD22d_2025_06_20_P1', 'TD22d_2025_06_20_P2' ; ...
 @loadTD22_neural, '2025-06-21', 'TD22d_2025_06_21_P1', 'TD22d_2025_06_21_P1' ; ...
 @loadTD23_neural, '2025-06-17', 'TD23d_2025_06_17_P1', 'TD23d_2025_06_17_P2' ; ...
 @loadTD23_neural, '2025-06-18', 'TD23d_2025_06_18_P1', 'TD23d_2025_06_18_P2' ; ...
 @loadTD23_neural, '2025-06-19', 'TD23d_2025_06_19_P1', 'TD23d_2025_06_19_P2' ; ...
 @loadTD23_neural, '2025-06-20', 'TD23d_2025_06_20_P1', 'TD23d_2025_06_20_P2' ; ...
 @loadTD23_neural, '2025-06-21', 'TD23d_2025_06_21_P1', 'TD23d_2025_06_21_P2' };

nSessions = size(sessionTable, 1);

all_meta = [];
for s = 1:nSessions
    all_meta = [all_meta; slimMeta(sessionTable{s,1}, sessionTable{s,2})];   %#ok<AGROW>
end

dataDirs = cell(nSessions, 2);
for s = 1:nSessions
    dataDirs{s,1} = fullfile(cfg.hmmRoot, sessionTable{s,3});
    dataDirs{s,2} = fullfile(cfg.hmmRoot, sessionTable{s,4});
end

%% PROBE MAP (defined ONCE)
% recorded in that session. Derived from spec.groupMaps in the tongue-length
% decoding script, matched by animal and date.
% Position i refers to session i, so a length mismatch silently reassigns
% probes -- hence the assert.
% Probe maps cross-checked against tongue_r14.m (spec.groupMaps, matched by
% animal and date): all 25 sessions and every M1/ALM probe agree. The one
% change is session 25 (TD23d 2025-06-21) ALM: tongue_r14 puts ALM on probe 2
% (probe 1 carries no region that day), but the only HMM-GLM fit for this
% session is _P1, so ALM spikes were being aligned to transition times from
% Session 25 ALM now uses the TD23d_2025_06_21_P2 fit. ALM: 15 sessions, 4 animals.
m1  = [1 1 0 1  2 0 2 1  1 2 2  2 1 1 2  0 1 2 2 1  2 2 2 1 0];
alm = [2 2 2 2  0 2 0 2  0 0 0  0 0 0 0  2 2 1 1 0  1 1 1 2 2];   % session 25 ALM restored to probe 2 (a _P2 fit now exists)

assert(numel(m1) == nSessions && numel(alm) == nSessions, ...
    ['probe map has %d (m1) / %d (alm) entries but there are %d sessions. ' ...
     'Position i refers to session i, so a mismatch reassigns probes silently.'], ...
    numel(m1), numel(alm), nSessions);

regions = struct('name', {'M1','ALM'}, 'map', {m1, alm});

%% SETTINGS PRINTOUT
fprintf('\n%s\n', repmat('=',1,72));
logf(' ANALYSIS SETTINGS\n');
fprintf('%s\n', repmat('-',1,72));
logf('  params.dt      = %g s  (%.0f Hz)\n', params.dt, 1/params.dt);
logf('  params.lowFR   = %g Hz\n', params.lowFR);
logf('  params.smooth  = %g\n', params.smooth);
logf('  params.quality = %s\n', strjoin(strtrim(params.quality), ', '));
logf('  engagement-mode window = +/- %g ms   (Methods say 400)\n', 1000*cfg.engModeWin_s);
logf('  projection     = %s\n', cfg.projMode);
fprintf('  unit test      = %s, alpha = %.4g%s\n', cfg.unitTest, cfg.unitAlpha, ...
    ternaryStr(cfg.unitFDR, ' (BH-FDR)', ' (uncorrected)'));
for r = 1:numel(regions)
    fprintf('  %-4s sessions used = %d of %d\n', regions(r).name, ...
        sum(regions(r).map ~= 0), nSessions);
end
fprintf('%s\n\n', repmat('=',1,72));

%% dataDir / PROBE CONSISTENCY CHECK
% The HMM folder names end in _P1 or _P2. For every (session, region) the plots
% actually use, check that the suffix matches the probe the map asks for. A
% mismatch means that region's units are being aligned with the OTHER probe's
% transition times, which is silent and wrong.
fprintf('%s\n', repmat('-',1,72));
logf(' dataDir / probe consistency\n');
nMismatch = 0;
for r = 1:numel(regions)
    for s = 1:nSessions
        pc = regions(r).map(s);
        if pc == 0, continue; end
        leaf = sessionTable{s, 2 + pc};
        tok  = regexp(leaf, '_P(\d)$', 'tokens', 'once');
        if isempty(tok), continue; end
        if str2double(tok{1}) ~= pc
            nMismatch = nMismatch + 1;
            logf('  MISMATCH sess %2d %-5s: map says probe %d, folder is %s\n', ...
                s, regions(r).name, pc, leaf);
        end
    end
end
if nMismatch == 0
    logf('  all used (session, region) pairs match their folder suffix\n');
else
    fprintf('  %d mismatch(es) above -- those regions are aligned with the\n', nMismatch);
    fprintf('  other probe''s HMM transitions. Check before trusting them.\n');
end
fprintf('%s\n\n', repmat('-',1,72));

%% PREALLOCATE
r1top = cell(nSessions, 2);   r1bot = cell(nSessions, 2);
r4top = cell(nSessions, 2);   r4bot = cell(nSessions, 2);
All_R1Modes_trials = cell(nSessions, 2);
All_R4Modes_trials = cell(nSessions, 2);

n_tests      = zeros(nSessions, 2);
n_sig        = zeros(nSessions, 2);
n_considered = zeros(nSessions, 2);
n_trials_R1  = zeros(nSessions, 2);

%% MAIN LOOP
for sessnum = 1:nSessions

    meta = all_meta(sessnum, 1);
    params.probe = {meta.probe};
    params.cluid = {};
    [obj, params] = slimToLegacy(meta, params);

    [reg1, reg2] = regionSplit(obj(1), params(1));
    allRegions = {reg1, reg2};

    for aa = 1:2
% Only probes the maps use (same rule as tongue_r14's probesOf). Unmapped
% probes -- e.g. tjS1 on probe 1 in sessions 10 and 15 -- were being
% aligned to the other probe's HMM folder and then discarded.
        if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end
        brainRegion = allRegions{aa};
        if isempty(brainRegion)
            logf('[skip] sess %2d region %d: no units (single-probe session)\n', sessnum, aa);
            continue
        end

% ---------------- HMM tables ----------------
        HMM = readHMM(dataDirs{sessnum, aa}, ...
                      {'R1_Trial_Track','R4_Trial_Track','R4_dt','R1_dt'});
        trials1 = table2array(HMM.R1_Trial_Track);  trials1 = trials1(:);
        trials4 = table2array(HMM.R4_Trial_Track);  trials4 = trials4(:);
%% HMM dt files are in milliseconds. Convert to bins at THIS script's own
%% *100, which are only correct at params.dt = 1/100.
        dtBins1 = round(table2array(HMM.R1_dt) * 0.001 / params.dt);  dtBins1 = dtBins1(:);
        dtBins4 = round(table2array(HMM.R4_dt) * 0.001 / params.dt);  dtBins4 = dtBins4(:);
% Sanity check on the HMM dt units (ms from the go cue). Compare these
% against the Fig. 2K / S6B disengagement times: if they come out ~10x
% too small the files are in 10-ms HMM bins, and if they are ~140 ms off
% the export is relative to the HMM window start (GC - 140 ms), not the GC.
        logf('  [hmm] sess %2d probe %d | median transition re: GC  R1 %.3f s  R4 %.3f s  (%d / %d trials)\n', ...
            sessnum, aa, median(dtBins1,'omitnan')*params.dt, median(dtBins4,'omitnan')*params.dt, ...
            numel(trials1), numel(trials4));

% ---------------- drop NaN-shift trials ----------------
        bad1 = isnan(dtBins1);
        if any(bad1)
            warning('Sess %d reg %d: %d R1 trials with NaN shift removed.', sessnum, aa, sum(bad1));
            trials1(bad1) = [];  dtBins1(bad1) = [];
        end
        bad4 = isnan(dtBins4);
        if any(bad4)
            warning('Sess %d reg %d: %d R4 trials with NaN shift removed.', sessnum, aa, sum(bad4));
            trials4(bad4) = [];  dtBins4(bad4) = [];
        end

% ---------------- align to the transition ----------------
        allDat = obj.trialdat(:, brainRegion, :);
        T = size(allDat, 1);
        N = size(allDat, 2);
        n_considered(sessnum, aa) = N;
        n_trials_R1(sessnum, aa)  = numel(trials1);

        aligned1 = doAlign(allDat, trials1, dtBins1);
        aligned4 = doAlign(allDat, trials4, dtBins4);

        psth1 = mean(aligned1(:,:,trials1), 3, 'omitnan');   % T x N
        psth4 = mean(aligned4(:,:,trials4), 3, 'omitnan');

% ---------------- the two windows, FROM THE TIME AXIS ----------------
% [-W, 0) and [0, +W): equal length, no overlap, split exactly at the
% transition. The old hard-coded 170:200 / 200:230 were off by one and
% shared bin 200.
        tAxis = params.tmin + (0:T-1)*params.dt;
        win1  = find(tAxis >= -cfg.engModeWin_s & tAxis <  0);
        win2  = find(tAxis >=  0 & tAxis <  cfg.engModeWin_s);
        if sessnum == 1 && aa == 1
            logf('  [windows] pre = bins %d:%d (%.3f to %.3f s, %d bins)\n', ...
                win1(1), win1(end), tAxis(win1(1)), tAxis(win1(end)), numel(win1));
            logf('  [windows] post = bins %d:%d (%.3f to %.3f s, %d bins)\n\n', ...
                win2(1), win2(end), tAxis(win2(1)), tAxis(win2(end)), numel(win2));
        end

% ---------------- engagement mode from R4 ----------------
        w = mean(psth4(win1,:), 1) - mean(psth4(win2,:), 1);   % 1 x N
        w = w / sum(abs(w));   % sum|w| = 1

        mod1 = projectMode(psth1, w, cfg.projMode);   % T x 1
        mod4 = projectMode(psth4, w, cfg.projMode);

% ---------------- per-unit test on R1 trials ----------------
        nT1   = numel(trials1);
        pvals = nan(N, 1);
        for neuron = 1:N
            trialData   = squeeze(aligned1(:, neuron, trials1));   % T x nT1
            pre_trials  = mean(trialData(win1,:), 1);
            post_trials = mean(trialData(win2,:), 1);
            pvals(neuron) = pairedTest(pre_trials, post_trials, cfg.unitTest);
        end

        n_tests(sessnum, aa) = N;
        if cfg.unitFDR
            sigNeurons = find(bhFDRlocal(pvals) < cfg.unitAlpha);
        else
            sigNeurons = find(pvals < cfg.unitAlpha);
        end
        n_sig(sessnum, aa) = numel(sigNeurons);

% ---------------- single-trial projections (always stored) ----------
        proj1 = nan(T, nT1);
        for k = 1:nT1
            proj1(:,k) = projectMode(aligned1(:,:,trials1(k)), w, cfg.projMode);
        end
        proj4 = nan(T, numel(trials4));
        for k = 1:numel(trials4)
            proj4(:,k) = projectMode(aligned4(:,:,trials4(k)), w, cfg.projMode);
        end

        All_R1Modes_trials{sessnum,aa} = proj1;
        All_R4Modes_trials{sessnum,aa} = proj4;

        if isempty(sigNeurons)
            warning('Sess %d reg %d: no units below alpha = %.4g.', sessnum, aa, cfg.unitAlpha);
            continue
        end

% ---------------- split by weight sign, sort by |w| ----------------
        topIdx = sigNeurons(w(sigNeurons) >  0);   % decrease at the transition
        botIdx = sigNeurons(w(sigNeurons) <= 0);   % increase
        [~, o] = sort(w(topIdx), 'descend');  topIdx = topIdx(o);
        [~, o] = sort(w(botIdx), 'ascend');   botIdx = botIdx(o);

        r1top{sessnum,aa} = psth1(:, topIdx).';   % units x time
        r4top{sessnum,aa} = psth4(:, topIdx).';
        r1bot{sessnum,aa} = psth1(:, botIdx).';
        r4bot{sessnum,aa} = psth4(:, botIdx).';
    end
end

%% SUMMARY
n_considered_plot = 0;  n_tests_plot = 0;  n_sig_plot = 0;  minTrialsUsed = Inf;
for r = 1:numel(regions)
    for s = 1:nSessions
        pc = regions(r).map(s);
        if pc == 0, continue; end
        n_considered_plot = n_considered_plot + n_considered(s, pc);
        n_tests_plot      = n_tests_plot      + n_tests(s, pc);
        n_sig_plot        = n_sig_plot        + n_sig(s, pc);
        if n_trials_R1(s, pc) > 0
            minTrialsUsed = min(minTrialsUsed, n_trials_R1(s, pc));
        end
    end
end

fprintf('\n%s\n', repmat('=',1,72));
fprintf(' SUMMARY (sessions/probes that reach the plots)\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('  Test               : %s%s\n', testLabel(cfg.unitTest), ...
    ternaryStr(cfg.unitFDR, ', BH-FDR across units', ', uncorrected'));
fprintf('  Comparison         : R1 trials, pre-window vs post-window at t = 0\n');
fprintf('  Window             : +/- %g ms\n', 1000*cfg.engModeWin_s);
fprintf('  Alpha              : %.4g\n', cfg.unitAlpha);
fprintf('  Neurons considered : %d\n', n_considered_plot);
fprintf('  Tests performed    : %d\n', n_tests_plot);
if n_tests_plot > 0
    fprintf('  Neurons passing    : %d  (%.1f%%)\n', n_sig_plot, 100*n_sig_plot/n_tests_plot);
end
if ~cfg.unitFDR
    fprintf('  Expected FP        : %.1f\n', n_tests_plot * cfg.unitAlpha);
else
    fprintf('  Expected FP        : n/a (FDR controls the proportion, not a count)\n');
end

% THE ALPHA FLOOR. The exact signed-rank test cannot produce a two-tailed p
% below 2^(1-n) with n paired observations, so on a session with too few R1
% trials NO unit can reach alpha however large its effect -- silently, and only
if strcmpi(cfg.unitTest, 'signrank')
    nNeeded = ceil(1 - log2(cfg.unitAlpha));
    fprintf('%s\n', repmat('-',1,72));
    fprintf('  signed-rank alpha floor: needs n >= %d R1 trials for p < %.4g\n', ...
        nNeeded, cfg.unitAlpha);
    if isfinite(minTrialsUsed)
        logf('  smallest R1 trial count in use: %d\n', minTrialsUsed);
    else
        logf('  smallest R1 trial count in use: n/a (no region produced trials)\n');
    end
    if isfinite(minTrialsUsed) && minTrialsUsed < nNeeded
        fprintf('  *** at least one session cannot produce a significant unit ***\n');
        for r = 1:numel(regions)
            for s = 1:nSessions
                pc = regions(r).map(s);
                if pc == 0 || n_trials_R1(s,pc) == 0, continue; end
                if n_trials_R1(s,pc) < nNeeded
                    fprintf('      sess %2d %-4s: %d R1 trials\n', s, regions(r).name, n_trials_R1(s,pc));
                end
            end
        end
    end
end
fprintf('%s\n\n', repmat('=',1,72));

%% ASSEMBLE PLOT DATA (once)
% Both regions pooled: every (session, region) pair the maps select contributes
% its units. use_group is the 2*nSessions key the old code built by hand.
use_group = [alm(:).', m1(:).'];
stkR1t = [r1top; r1top];   % row s = ALM pass, row s+nSessions = M1 pass
stkR1b = [r1bot; r1bot];
stkR4t = [r4top; r4top];
stkR4b = [r4bot; r4bot];

widths = [];
for s = 1:numel(use_group)
    pc = use_group(s);
    if pc == 0, continue; end
    for c = {stkR1t{s,pc}, stkR1b{s,pc}}
        if ~isempty(c{1}), widths(end+1) = size(c{1}, 2); end   %#ok<AGROW>
    end
end
assert(~isempty(widths), 'No sessions produced any significant units.');
nT = min(widths);

t = params.tmin + (0:nT-1)*params.dt;
[~, t0bin] = min(abs(t));
logf('Plot time axis: %.3f to %.3f s | nT = %d | t = 0 at bin %d\n', t(1), t(end), nT, t0bin);

r1_top_all = [];  r4_top_all = [];
r1_bot_all = [];  r4_bot_all = [];
for s = 1:numel(use_group)
    pc = use_group(s);
    if pc == 0, continue; end
    if ~isempty(stkR1t{s,pc}) && ~isempty(stkR4t{s,pc})
        r1_top_all = [r1_top_all; stkR1t{s,pc}(:, 1:nT)];   %#ok<AGROW>
        r4_top_all = [r4_top_all; stkR4t{s,pc}(:, 1:nT)];   %#ok<AGROW>
    end
    if ~isempty(stkR1b{s,pc}) && ~isempty(stkR4b{s,pc})
        r1_bot_all = [r1_bot_all; stkR1b{s,pc}(:, 1:nT)];   %#ok<AGROW>
        r4_bot_all = [r4_bot_all; stkR4b{s,pc}(:, 1:nT)];   %#ok<AGROW>
    end
end

nTop = size(r1_top_all, 1);
nBot = size(r1_bot_all, 1);
nNeurons = nTop + nBot;
fprintf('Decreasing units (positive w): %d\n', nTop);
fprintf('Increasing units (negative w): %d\n', nBot);
fprintf('Total units in the plots: %d\n\n', nNeurons);

% Sort: decreasing group by peak BEFORE t = 0, increasing group by peak AFTER.
pre_win  = t >= -cfg.sortWin_s & t <  0;
post_win = t >   0             & t <= cfg.sortWin_s;
[~, oTop] = sort(max(r1_top_all(:, pre_win),  [], 2), 'descend');
[~, oBot] = sort(max(r1_bot_all(:, post_win), [], 2), 'ascend');

data_r1 = [r1_top_all(oTop,:); r1_bot_all(oBot,:)];
data_r4 = [r4_top_all(oTop,:); r4_bot_all(oBot,:)];

ttlSuffix = sprintf('%s, alpha = %.3g%s', cfg.unitTest, cfg.unitAlpha, ...
    ternaryStr(cfg.unitFDR, ' (BH-FDR)', ' (uncorrected)'));

%% FIGURE 1: percentile colour scale
clim1 = [0, prctile([data_r1(:); data_r4(:)], cfg.climPctile)];
drawPairedHeatmap(t, data_r1, data_r4, nTop, nBot, clim1, cfg.xlimPlot, ...
    sprintf('M1 + ALM | n = %d | %s | colour to %gth pct', nNeurons, ttlSuffix, cfg.climPctile));

% Only the percentile-scaled heatmap is produced. The fixed-colour-scale
% version and the per-region engagement-mode trace figures were removed,
% along with the modeTrace assembly that only they used.
%% LOCAL FUNCTIONS
function p = pairedTest(a, b, which)
% One paired test, selected by name. a and b are matched samples.
    switch lower(which)
        case 'signrank'
            p = signrank(a, b);   % two-tailed by default
        case 'ttest2t'
            [~, p] = ttest(a, b);
        case 'ttest1t'
% The Fig. 2L legend's version: one-tailed in BOTH directions, so a
% unit counts if it changes significantly either way.
            [~, pUp] = ttest(a, b, 'Tail', 'left');
            [~, pDn] = ttest(a, b, 'Tail', 'right');
            p = min(pUp, pDn);
        otherwise
            error('cfg.unitTest must be ''signrank'', ''ttest2t'' or ''ttest1t''; got ''%s''.', which);
    end
end

function s = testLabel(which)
    switch lower(which)
        case 'signrank', s = 'Paired two-tailed Wilcoxon signed-rank across trials';
        case 'ttest2t',  s = 'Paired two-tailed t-test across trials';
        case 'ttest1t',  s = 'Paired one-tailed t-test, both directions';
        otherwise,       s = which;
    end
end

function aligned = doAlign(data, trials, dtBins)
% Shift each listed trial by its own dtBins so t = 0 becomes the transition.
% Samples shifted in from outside the trial are left NaN.
    [T, N, TT] = size(data);
    aligned = nan(T, N, TT);
    for k = 1:numel(trials)
        tr = trials(k);
        sh = dtBins(k);
        X  = data(:,:,tr);
        if sh >= 0
            aligned(1:T-sh, :, tr) = X(1+sh:end, :);
        else
            aligned(1-sh:T, :, tr) = X(1:end+sh, :);
        end
    end
end

function [reg1, reg2] = regionSplit(obj, params)
% Cluster indices for the two probes, with BOTH always assigned. A single-probe
% session returns an empty second region rather than whatever the previous
% session left behind.
    Ncells = size(obj.trialdat, 2);
    if iscell(params.cluid) && ~isempty(params.cluid)
        n1 = size(params.cluid{1,1}, 1);
    else
        n1 = Ncells;
    end
    n1   = min(n1, Ncells);
    reg1 = 1:n1;
    if n1 >= Ncells
        reg2 = [];
    else
        reg2 = (n1+1):Ncells;
    end
end

function HMM = readHMM(dataDir, fileBases)
% Read the HMM-GLM tables into a STRUCT, rebuilt per session and per region.
% name, which silently reuses the previous session's table if a read is skipped.
    HMM = struct();
    for k = 1:numel(fileBases)
        fn = fullfile(dataDir, [fileBases{k} '.csv']);
        if ~isfile(fn), error('Missing HMM result file: %s', fn); end
        HMM.(fileBases{k}) = readtable(fn);
    end
end

function p = projectMode(X, w, mode)
% Project activity onto the engagement mode.
%   Methods:  p(t) = sum_i w_i * r_i(t),  with sum_i |w_i| = 1
% 'mean' divides that by the unit count. Because the weights are already
% normalised, that extra division rescales each session by its own N, which
% matters as soon as sessions are averaged together.
    N = size(X, 2);
    switch lower(mode)
        case 'sum',  p = sum(X .* reshape(w,1,N), 2);
        case 'mean', p = mean(X .* reshape(w,1,N), 2);
        otherwise,   error('cfg.projMode must be ''sum'' or ''mean''; got ''%s''.', mode);
    end
end

function drawPairedHeatmap(t, dR1, dR4, nTop, nBot, climVal, xl, ttl)
% R1 over R4, same colour scale, ONE colorbar spanning the full height of both
% panels, and a dashed line between the decreasing and increasing blocks.
% The panels are positioned MANUALLY rather than with subplot(), because a
% colorbar has to be told where to go to span two axes: colorbar(ax) attaches
% itself to that one axes and sizes itself to it, which is why the bar used to
% sit beside the lower panel only. Creating it also SHRINKS the axes it is
% attached to, so both positions are reset afterwards and the bar is then
% stretched from the bottom of the lower panel to the top of the upper one.
    n = nTop + nBot;

% , normalised figure units
    L     = 0.11;   % left edge, both panels
    W     = 0.70;   % width, both panels
    Hp    = 0.36;   % height of each panel
    Bl    = 0.085;   % bottom of the LOWER panel
    Bu    = 0.525;   % bottom of the UPPER panel
    cbGap = 0.025;   % gap between the panels and the colorbar
    cbW   = 0.030;   % colorbar width

    figure('Color','w', 'Name',ttl, 'Units','normalized', ...
           'Position',[0.30 0.06 0.42 0.84]);

    pos = {[L Bu W Hp], [L Bl W Hp]};
    ax  = gobjects(1,2);
    for k = 1:2
        ax(k) = axes('Position', pos{k});   %#ok<LAXES>
        if k == 1
            D = dR1;  lbl = 'R1';
        else
            D = dR4;  lbl = 'R4';
        end
        imagesc(ax(k), t, 1:n, D);
        set(ax(k), 'YDir','reverse', 'Color','w', 'FontSize',11);
        caxis(ax(k), climVal);
        colormap(ax(k), safeColormap());
        xlim(ax(k), xl);
        ylabel(ax(k), 'Neurons');
        title(ax(k), sprintf('%s trials', lbl), 'FontWeight','bold');
        hold(ax(k), 'on')
% of those is version-fragile, line() works everywhere.
        line(ax(k), xl, [nTop nTop] + 0.5, 'Color','w', 'LineStyle','--', ...
            'LineWidth', 1.5, 'HandleVisibility','off');
        line(ax(k), [0 0], [0.5 n+0.5], 'Color','w', ...
            'LineWidth', 2, 'HandleVisibility','off');
        text(ax(k), xl(1)+0.03, max(nTop/2,1), [lbl char(8595)], ...
            'Color','w', 'FontSize',12, 'FontWeight','bold');
        text(ax(k), xl(1)+0.03, nTop + max(nBot/2,1), [lbl char(8593)], ...
            'Color','w', 'FontSize',12, 'FontWeight','bold');
        hold(ax(k), 'off')
        box(ax(k), 'off');
    end

% Only the lower panel carries the time axis; the two are linked, so the
% upper one's tick labels would just repeat it.
    set(ax(1), 'XTickLabel', {});
    xlabel(ax(2), 'Time from disengagement (s)');

% ONE colorbar, spanning both panels.
    cb = colorbar(ax(2));
    set(ax(1), 'Position', pos{1});   % undo the shrink colorbar() applies
    set(ax(2), 'Position', pos{2});
    cb.Position     = [L + W + cbGap, Bl, cbW, (Bu + Hp) - Bl];
    cb.Label.String = 'Firing rate (Hz)';
    cb.FontSize     = 11;

    linkaxes(ax, 'xy');
    sgtitle(ttl, 'FontSize', 10);
end

function q = bhFDRlocal(p)
% Benjamini-Hochberg adjusted p-values across units; NaNs preserved in place.
    p  = p(:);
    ok = isfinite(p);
    q  = nan(size(p));
    pv = p(ok);
    m  = numel(pv);
    if m == 0, return; end
    [ps, ix] = sort(pv);
    qa = nan(m,1);
    for i = 1:m
        qa(ix(i)) = min(ps(i:end) .* m ./ (i:m)');
    end
    q(ok) = min(qa, 1);
end

function cm = safeColormap()
% linspecer if it is on the path, parula otherwise, so the figures still draw
% on a machine that does not have it.
    if exist('linspecer', 'file') == 2
        cm = linspecer;
    else
        cm = parula(256);
    end
end

function s = ternaryStr(c, a, b)
    if c, s = a; else, s = b; end
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

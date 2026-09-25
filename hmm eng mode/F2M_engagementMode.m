%% F2M_engagementMode.m
%  Engaged and disengaged motor cortical states, Delayed Reward Task.
%  Per-trial disengagement times come from the HMM-GLM fits in the
%  Disengagement Times folder. Neural activity is aligned to those times to
%  define the engagement mode and measure how fast the state switches.
%  Produces Fig. 2M.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%    Disengagement Times\<session>\   HMM-GLM state fits
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/300
%    params.smooth      15
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


sz = 26;

params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping only operates on neural data for now.
params.behav_only          = 0;

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

% ---- THE PARAMETERS BLOCK, RESTORED ---------------------------------------
% A previous pass deleted the "duplicate" parameters block here on the grounds
% that params.condition(1) = {...} overwrote element 1 without clearing the
% array, appending nine more entries onto ten that already existed and leaving
% lowFR, smooth, quality, tmin, tmax, condition and traj_features were never
% set at all, so this script errored on the first reference to params.dt and
% could not run.
% The block is restored below, written as a single cell literal so the
% overwrite-versus-append problem cannot come back. Values are the ones this
% file's header documents for engMode: dt = 1/300, lowFR = 0.2, smooth = 10.
params.timeWarp = 0;
params.nLicks   = 20;
params.lowFR    = 0.01;   % minimum mean firing rate, Hz
params.tmin     = -2;
params.tmax     = 5;
params.dt       = 1/300;
params.smooth   = 15;
params.quality  = {'good'};   % good units only (findClusters trims blanks and ignores case)

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

%% CONSISTENCY CONFIG
% One place for the settings that were hard-coded in several spots, or that
% differ between the five scripts.

% ENGAGEMENT MODE. The Methods describe a 400 ms window either side of the
% transition; the code has always used 300 ms. Set this to 0.400 to follow the
% Methods, or change the Methods to 300 ms -- but the two have to agree.
cfg.engModeWin_s = 0.400;   % +/-400 ms around the transition

% PROJECTION. Methods: p(t) = sum_i w_i r_i(t), weights normalised so
% sum|w_i| = 1. 'mean' divides that by the unit count, which is what the code
cfg.projMode = 'sum';   % 'sum' (Methods) | 'mean' (previous behaviour)

%% SPECIFY DATA TO LOAD
datapth = '';   % raw data folder not used (was: datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';)

%% SESSIONS (one table, one source of truth)
% {loader, date, HMM folder for probe 1, HMM folder for probe 2}
% so slots 6, 7 and 8 combined one day's spikes with another day's transition
% times and a third day's probe assignment, silently.
% The list below is the one the tongue-length decoding script uses, matched by
% animal and date. TD4d is 02-21, 02-24, 02-25, 03-19.
% Where the two HMM folders on a row are the same, only one result exists for
% that session and both regions are aligned with it -- the consistency check
% below reports when that disagrees with the probe map.
hmmRoot = fullfile(v2Root, 'Disengagement Times');   % HMM-GLM disengagement times, one folder per session and probe
assert(exist(hmmRoot, 'dir') == 7, 'No disengagement-time folder: %s', hmmRoot);

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

% 25 numbered meta variables replaced by one loop.
all_meta = [];
for s = 1:nSessions
    all_meta = [all_meta; slimMeta(sessionTable{s,1}, sessionTable{s,2})];   %#ok<AGROW>
end

dataDirs = cell(nSessions, 2);
for s = 1:nSessions
    dataDirs{s,1} = fullfile(hmmRoot, sessionTable{s,3});
    dataDirs{s,2} = fullfile(hmmRoot, sessionTable{s,4});
end

%% PROBE MAP (defined ONCE)
% Which probe carries each region, per session. 0 = region not recorded that
% session. Position i refers to session i, so a length mismatch silently
% reassigns probes -- hence the assert.
% Probe maps cross-checked against tongue_r14.m (spec.groupMaps, matched by
% animal and date): all 25 sessions and every M1/ALM probe agree. One change:
% Session 25 (TD23d 2025-06-21) ALM is back on probe 2: the HMM Alignments
% folder now holds a TD23d_2025_06_21_P2 fit, so ALM is aligned to its own
% probe's transitions. HMM panels: M1 21 sessions / 6 animals, ALM 15 sessions /
% 4 animals.
m1  = [1 1 0 1  2 0 2 1  1 2 2  2 1 1 2  0 1 2 2 1  2 2 2 1 0];
alm = [2 2 2 2  0 2 0 2  0 0 0  0 0 0 0  2 2 1 1 0  1 1 1 2 2];   % session 25 ALM restored to probe 2 (a _P2 fit now exists)

assert(numel(m1) == nSessions && numel(alm) == nSessions, ...
    ['probe map has %d (m1) / %d (alm) entries but there are %d sessions. ' ...
     'Position i refers to session i, so a mismatch reassigns probes silently.'], ...
    numel(m1), numel(alm), nSessions);

regions = struct('name', {'M1','ALM'}, 'map', {m1, alm});

%% dataDir / probe consistency check
% The HMM folder names end in _P1 or _P2. For every (session, region) the
% figures actually use, check the suffix matches the probe the map asks for. A
% mismatch means that region is aligned with the OTHER probe's transition times.
fprintf('%s\n', repmat('-',1,72));
logf(' session list: %d sessions | dataDir / probe consistency\n', nSessions);
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
            logf('  MISMATCH sess %2d %-4s: map says probe %d, folder is %s\n', ...
                s, regions(r).name, pc, leaf);
        end
    end
end
if nMismatch == 0
    logf('  all used (session, region) pairs match their folder suffix\n');
end
for r = 1:numel(regions)
    fprintf('  %-4s sessions used: %d of %d\n', regions(r).name, ...
        sum(regions(r).map ~= 0), nSessions);
end
fprintf('%s\n\n', repmat('-',1,72));

% Define a common time axis relative to disengagement (in seconds)
tPre_dis  = -2.0;   % seconds before disengagement
tPost_dis =  2.0;   % seconds after disengagement
nBins_dis = round((tPost_dis - tPre_dis) / params.dt) + 1;
tAxis_dis = linspace(tPre_dis, tPost_dis, nBins_dis);

% Pre-allocate save structures
AllProj_R1     = cell(length(all_meta), 2);
AllProj_R4     = cell(length(all_meta), 2);
AllDtBins_R1   = cell(length(all_meta), 2);   % save dt offsets for alignment
AllDtBins_R4   = cell(length(all_meta), 2);
All_R1Modes    = cell(length(all_meta), 2);
All_R4Modes    = cell(length(all_meta), 2);
AllTime        = cell(length(all_meta), 2);

params.behav_only = 0;

%% MAIN LOOP
for sessnum = 1:length(all_meta)

    clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin

    meta = all_meta(sessnum,1);

    params.probe = {meta.probe};
    params.cluid = {};

    [obj, params, kin] = slimToLegacy(meta, params);

    for sessix = 1:numel(meta)
    end

    nSessions = numel(meta);
    for sessix = 1:numel(meta)
        message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), ...
            {' '},'out of',{' '},num2str(nSessions),'----');
        disp(message)
    end

    sessix = 1;

% Build allRegions from cluid
    [reg1, reg2, isSingleProbe] = regionSplit(obj(1), params(1));   % was: params.cluid{1,1} with no guard
    allRegions = {reg1, reg2};

    for aa = 1:2

% analysed as M1 and the real M1 probe skipped; on 10 and 15 probe 1 is tjS1.
        if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end

        brainRegion = allRegions{aa};
        dataDir     = dataDirs{sessnum, aa};

        fileBases = {'R1_Trial_Track','R4_Trial_Track','R4_dt','R1_dt'};
        HMM = readHMM(dataDir, fileBases);   % was: assignin into the base workspace

        r1Trials = table2array(HMM.R1_Trial_Track);
        r4Trials = table2array(HMM.R4_Trial_Track);

% dt files in 0.001s units -> neural bins
        r1DtsIdx = round(table2array(HMM.R1_dt) * 0.001 / params.dt);
        r4DtsIdx = round(table2array(HMM.R4_dt) * 0.001 / params.dt);
% HMM dt files: ms from the go cue (confirmed against Fig. 2K / S6B via the
% bigPlot session log). Printed so every run shows the same check.
        logf('  [hmm] sess %2d probe %d | median transition re: GC  R1 %.3f s  R4 %.3f s\n', ...
            sessnum, aa, median(r1DtsIdx,'omitnan')*params.dt, median(r4DtsIdx,'omitnan')*params.dt);

% SETUP
        neurons = brainRegion;
        allDat  = obj(sessix).trialdat(:, neurons, :);
        T       = size(allDat, 1);
        N       = size(allDat, 2);

        trials1 = r1Trials(:);
        dtBins1 = r1DtsIdx(:);
        trials4 = r4Trials(:);
        dtBins4 = r4DtsIdx(:);

        bad1 = isnan(dtBins1);
        if any(bad1)
            warning('%d R1-trials have NaN shift -> removing them', sum(bad1));
            trials1(bad1) = [];  dtBins1(bad1) = [];
        end

        bad4 = isnan(dtBins4);
        if any(bad4)
            warning('%d R4-trials have NaN shift -> removing them', sum(bad4));
            trials4(bad4) = [];  dtBins4(bad4) = [];
        end

% Align neural data to go cue to compute the mode weight vector
        aligned1 = doAlign(allDat, trials1, dtBins1);
        aligned4 = doAlign(allDat, trials4, dtBins4);

        nT1 = numel(trials1);
        nT4 = numel(trials4);

        psth1 = nan(T, N);
        psth4 = nan(T, N);
        for nn = 1:N
            psth1(:,nn) = nanmean(aligned1(:,nn,trials1), 3);
            psth4(:,nn) = nanmean(aligned4(:,nn,trials4), 3);
        end

% ENGAGEMENT MODE WEIGHT VECTOR (from R4, go-cue aligned)
        win1 = obj(sessix).time >= -cfg.engModeWin_s & obj(sessix).time < 0;
        win2 = obj(sessix).time >= 0    & obj(sessix).time <= cfg.engModeWin_s;
        m1_4 = mean(psth4(win1,:), 1);
        m2_4 = mean(psth4(win2,:), 1);

        w = m1_4 - m2_4;
        w = w / sum(abs(w));

% PER-TRIAL PROJECTIONS (go-cue aligned, raw trial data)
% Then re-index each trial relative to its disengagement bin
        nHalf   = floor(nBins_dis / 2);
        preBins = round(abs(tPre_dis)  / params.dt);
        postBins= round(abs(tPost_dis) / params.dt);

        proj1_dis = nan(nBins_dis, nT1);
        proj4_dis = nan(nBins_dis, nT4);

        for k = 1:nT1
            tr      = trials1(k);
            dtBin   = dtBins1(k);   % bin offset from go cue to disengagement
            X       = allDat(:,:,tr);   % T x N, go-cue-aligned raw trial
            projFull= projectMode(X, w, cfg.projMode);   % T x 1

% Index of disengagement within this trial's time axis
            disIdx  = dtBin;   % dtBin bins after go cue (which is at idx_zero)
            [~, goIdx] = min(abs(obj(sessix).time - 0));
            absDisIdx  = goIdx + disIdx - 1;

% Extract window around disengagement
            iStart = absDisIdx - preBins;
            iEnd   = absDisIdx + postBins;

            if iStart < 1 || iEnd > T, continue; end
            proj1_dis(:, k) = projFull(iStart:iEnd);
        end

        for k = 1:nT4
            tr      = trials4(k);
            dtBin   = dtBins4(k);
            X       = allDat(:,:,tr);
            projFull= projectMode(X, w, cfg.projMode);

            [~, goIdx] = min(abs(obj(sessix).time - 0));
            absDisIdx  = goIdx + dtBin - 1;

            iStart = absDisIdx - preBins;
            iEnd   = absDisIdx + postBins;

            if iStart < 1 || iEnd > T, continue; end
            proj4_dis(:, k) = projFull(iStart:iEnd);
        end

% Trial-averaged disengagement-aligned projections
        mod1_dis = nanmean(proj1_dis, 2);
        mod4_dis = nanmean(proj4_dis, 2);

% SAVE
        All_R1Modes{sessnum, aa}  = mod1_dis;
        All_R4Modes{sessnum, aa}  = mod4_dis;
        AllProj_R1{sessnum, aa}   = proj1_dis;
        AllProj_R4{sessnum, aa}   = proj4_dis;
        AllDtBins_R1{sessnum, aa} = dtBins1;
        AllDtBins_R4{sessnum, aa} = dtBins4;
        AllTime{sessnum, aa}      = obj(sessix).time;

    end   % aa loop
end   % sessnum loop

%% POST-LOOP PLOTTING — disengagement-aligned, R1 vs R4, combined M1 + ALM

allMeans_R1 = [];
allMeans_R4 = [];

for sessnum = 1:length(all_meta)
    for aa = 1:2
        if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end   % aa = probe
        if isempty(AllProj_R1{sessnum,aa}), continue; end

        p1 = AllProj_R1{sessnum, aa};
        p4 = AllProj_R4{sessnum, aa};

        m1sess = nanmean(p1, 2);
        m4sess = nanmean(p4, 2);

        if all(isnan(m1sess)) || all(isnan(m4sess)), continue; end

        combined = [m1sess; m4sess];
        cMin     = nanmin(combined);
        cMax     = nanmax(combined);
        if cMax == cMin, continue; end

        allMeans_R1(:, end+1) = (m1sess - cMin) / (cMax - cMin);
        allMeans_R4(:, end+1) = (m4sess - cMin) / (cMax - cMin);
    end
end

% Smooth each session trace
for k = 1:size(allMeans_R1, 2)
    allMeans_R1(:,k) = smoothdata(allMeans_R1(:,k), 'gaussian', 15);
    allMeans_R4(:,k) = smoothdata(allMeans_R4(:,k), 'gaussian', 15);
end

nSess        = size(allMeans_R1, 2);
grandMean_R1 = nanmean(allMeans_R1, 2);
grandMean_R4 = nanmean(allMeans_R4, 2);

se_R1 = nanstd(allMeans_R1, 0, 2) / sqrt(nSess);
se_R4 = nanstd(allMeans_R4, 0, 2) / sqrt(nSess);

tval  = tinv(0.975, nSess - 1);
ci_R1 = tval * se_R1;
ci_R4 = tval * se_R4;

% Colors — flipped: R4 cyan-blue, R1 dark blue
col_R1 = [0.2  0.2  0.7];   % dark blue  -> R1
col_R4 = [0.0  0.6  1.0];   % cyan-blue  -> R4

% Shared y limits
ylo = min([grandMean_R1 - ci_R1; grandMean_R4 - ci_R4]) - 0.03;
yhi = max([grandMean_R1 + ci_R1; grandMean_R4 + ci_R4]) + 0.03;
ylo = max(ylo, 0);
yhi = min(yhi, 1);

% Scale bar sizes
xbar_left  = 0.250;
xbar_right = 0.050;
ybar       = 0.1;

figure('Color','w','Units','normalized','Position',[0.1 0.2 0.7 0.55]);

% LEFT PANEL
ax1 = axes('Position',[0.06 0.15 0.52 0.75]);
hold(ax1,'on');

fill(ax1, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R1 - ci_R1; flipud(grandMean_R1 + ci_R1)], ...
     col_R1, 'FaceAlpha',0.25, 'EdgeColor','none');
fill(ax1, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R4 - ci_R4; flipud(grandMean_R4 + ci_R4)], ...
     col_R4, 'FaceAlpha',0.25, 'EdgeColor','none');

plot(ax1, tAxis_dis, grandMean_R1, '-', 'Color',col_R1, 'LineWidth',2.5);
plot(ax1, tAxis_dis, grandMean_R4, '-', 'Color',col_R4, 'LineWidth',2.5);

xline(ax1, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

xlim(ax1, [-1.3 0.45]);
ylim(ax1, [ylo yhi]);
axis(ax1,'off');

% X scale bar (250 ms)
xbar_x = -1.25;
xbar_y = ylo + 0.01;
plot(ax1, [xbar_x, xbar_x + xbar_left], [xbar_y xbar_y], 'k-', 'LineWidth',2);
text(ax1, xbar_x + xbar_left/2, xbar_y - 0.04, '250 ms', ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

% Y scale bar (0.1 a.u.)
ybar_x = xbar_x - 0.04;
ybar_y = xbar_y + 0.01;
plot(ax1, [ybar_x ybar_x], [ybar_y, ybar_y + ybar], 'k-', 'LineWidth',2);
text(ax1, ybar_x - 0.06, ybar_y + ybar/2, sprintf('%.1f\n(a.u.)', ybar), ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

% State transition label
text(ax1, 0, ylo - 0.08, 'State transition', ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

hold(ax1,'off');

% RIGHT PANEL — zoomed
ax2 = axes('Position',[0.62 0.15 0.32 0.75]);
hold(ax2,'on');

% Gray background patch
patch(ax2, [-0.15 0.2 0.2 -0.15], [ylo ylo yhi yhi], ...
    [0.93 0.93 0.93], 'EdgeColor','none', 'FaceAlpha',1.0);

fill(ax2, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R1 - ci_R1; flipud(grandMean_R1 + ci_R1)], ...
     col_R1, 'FaceAlpha',0.25, 'EdgeColor','none');
fill(ax2, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R4 - ci_R4; flipud(grandMean_R4 + ci_R4)], ...
     col_R4, 'FaceAlpha',0.25, 'EdgeColor','none');

plot(ax2, tAxis_dis, grandMean_R1, '-', 'Color',col_R1, 'LineWidth',2.5);
plot(ax2, tAxis_dis, grandMean_R4, '-', 'Color',col_R4, 'LineWidth',2.5);

xline(ax2, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

xlim(ax2, [-0.15 0.2]);
ylim(ax2, [ylo yhi]);
axis(ax2,'off');

% X scale bar (50 ms)
xbar_x2 = -0.14;
xbar_y2 = ylo + 0.01;
plot(ax2, [xbar_x2, xbar_x2 + xbar_right], [xbar_y2 xbar_y2], 'k-', 'LineWidth',2);
text(ax2, xbar_x2 + xbar_right/2, xbar_y2 - 0.04, '50 ms', ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

hold(ax2,'off');

ylim(ax1, [ylo yhi]);
ylim(ax2, [ylo yhi]);

%% 90-TO-10% DECAY OF THE ENGAGEMENT MODE
% How fast the mode falls once it starts falling. Each session contributes two
% traces (R1 and R4). Within a +/-80 ms window around the transition, the peak
% and the following trough are found, the 90% and 10% levels of that drop are
% crossed by linear interpolation, and the decay is the time between the two
% crossings. Traces whose peak comes after their trough, or that never cross,
% are left out.

decayWin_s  = 0.080;                                  % +/- window searched
inWin       = tAxis_dis >= -decayWin_s & tAxis_dis <= decayWin_s;
tWin        = tAxis_dis(inWin);

traces = [allMeans_R1(inWin,:), allMeans_R4(inWin,:)]';   % (2*sessions) x bins
decay  = nan(size(traces,1), 1);

for k = 1:size(traces,1)
    a = traces(k,:);
    if all(isnan(a)), continue; end

    [aPeak,   iPeak  ] = max(a);
    [aTrough, iTrough] = min(a);
    if iPeak >= iTrough, continue; end

    drop = aPeak - aTrough;
    if drop < 1e-6, continue; end

    tFall = tWin(iPeak:end);
    aFall = a(iPeak:end);

    t90 = crossTime(tFall, aFall, aPeak - 0.10*drop);
    t10 = crossTime(tFall, aFall, aPeak - 0.90*drop);

    if ~isnan(t90) && ~isnan(t10) && t10 > t90
        decay(k) = t10 - t90;
    end
end

ok = ~isnan(decay);
n  = sum(ok);
mu = mean(decay(ok));
ci = std(decay(ok)) / sqrt(n) * tinv(0.975, n-1);

fprintf('\n90-to-10%% decay of the engagement mode\n');
fprintf('  %.1f +/- %.1f ms  (mean +/- 95%% CI, n = %d traces from %d sessions)\n', ...
    1000*mu, 1000*ci, n, nSess);

%% ADDED LOCAL FUNCTIONS

function [reg1, reg2, isSingle] = regionSplit(obj, params)
% Cluster indices for the two probes, with BOTH always assigned.
% The original indexed params.cluid{1,1} unconditionally, so a session recorded
% boundary that does not exist -- and the second region then silently held
% whatever the previous session left behind.
    Ncells = size(obj.trialdat, 2);
    isSingle = ~iscell(params.cluid) || isempty(params.cluid) || numel(params.cluid) < 2;
    if iscell(params.cluid) && ~isempty(params.cluid)
        n1 = size(params.cluid{1,1}, 1);
    else
        n1 = Ncells;
    end
    n1 = min(n1, Ncells);
    reg1 = 1:n1;
    if isSingle || n1 >= Ncells
        reg2 = [];
        isSingle = true;
    else
        reg2 = (n1+1):Ncells;
    end
end

function HMM = readHMM(dataDir, fileBases)
% Read the HMM-GLM result tables into a STRUCT.
% bare name. That silently reuses the previous session's table if a read is ever
% skipped, and it does not work inside a function at all.
    HMM = struct();
    for k = 1:numel(fileBases)
        fn = fullfile(dataDir, [fileBases{k} '.csv']);
        if ~isfile(fn)
            error('Missing HMM result file: %s', fn);
        end
        HMM.(fileBases{k}) = readtable(fn);
    end
end

function p = projectMode(X, w, mode)
% Project single-trial activity onto the engagement mode.
%   Methods:  p(t) = sum_i w_i * r_i(t),   sum_i |w_i| = 1
% 'mean' is the previous behaviour: the same sum divided by the number of
% units. Because the weights are already normalised, that extra division just
% rescales each session by its own unit count, which matters as soon as
% sessions are averaged together.
    N = size(X,2);
    switch lower(mode)
        case 'sum',  p = sum(X .* reshape(w,1,N), 2);
        case 'mean', p = mean(X .* reshape(w,1,N), 2);
        otherwise,   error('cfg.projMode must be ''sum'' or ''mean''.');
    end
end

function t = crossTime(tv, av, level)
% First time the falling trace av crosses down through level, interpolated
% between the two samples that straddle it. NaN if it never gets there.
    i = find(av <= level, 1);
    if isempty(i)
        t = NaN;
    elseif i == 1
        t = tv(1);
    else
        t = tv(i-1) + (level - av(i-1)) * (tv(i) - tv(i-1)) / (av(i) - av(i-1));
    end
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

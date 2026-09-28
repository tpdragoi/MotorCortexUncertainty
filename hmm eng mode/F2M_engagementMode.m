%% F2M_engagementMode.m
%  Engaged and disengaged motor cortical states, Delayed Reward Task.
%  Per-trial disengagement times come from the HMM-GLM fits in the
%  Disengagement Times folder. Neural activity is aligned to those times to
%  define the engagement mode and measure how fast the state switches.
%  Produces Fig. 2M.
%
%  Engagement mode: each neuron is weighted by d' (pre-minus-post difference
%  over the pooled across-trial SD of the two windows); see engagementModeWeights.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; Data = dataRoot in setPaths.m
%    Disengagement Times\<session>\   HMM-GLM state fits
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/300
%    params.smooth      15
%    params.quality     {'good','excellent'}
%    params.lowFR       0.01
%    params.window      -2 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  Run the whole file.
clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


sz = 26;

params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

params.behav_only          = 0;

% ---- paths ----
% Sessions are read from the exported obj/kin files in the data folder set in
% setPaths.m, through slimMeta / slimToLegacy / loadSlimSession (in shared/).
% slimToLegacy rebuilds obj/params/kin for this script's params (dt, window,
% alignment, smoothing, unit quality, lowFR, conditions) from the exported spike
% times. Pipeline functions it needs are in shared\pipelineCopies.
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
assert(exist(fullfile(repoRoot, 'setPaths.m'), 'file') == 2, ...
    'Cannot find setPaths.m. Run this script from its file, or cd to the repository root first.');
addpath(repoRoot);
cfgPaths = setPaths();   % data locations are set once, in setPaths.m

%% PARAMETERS
params.timeWarp = 0;
params.nLicks   = 20;
params.lowFR    = 0.01;   % minimum mean firing rate, Hz
params.tmin     = -2;
params.tmax     = 5;
params.dt       = 1/300;
params.smooth   = 15;
params.quality  = {'good','excellent',' good','good '};   % good + excellent units; exactly the export's list, so slimToLegacy uses the exported units directly

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

%% CONFIG

% Half-width of the pre and post windows around the transition, in s.
cfg.engModeWin_s = 0.400;   % +/-400 ms around the transition

% A neuron whose rate barely varies across trials has a near-zero pooled SD and
% would otherwise dominate the weight vector. Its SD is floored here, in
% spikes/s, before the division.
cfg.sdFloor = 0.5;


% Projection: p(t) = sum_i w_i r_i(t) with sum|w_i| = 1. 'mean' divides that by the
% unit count, which rescales every session by its own N.
cfg.projMode = 'sum';   % 'sum' | 'mean'

%% SPECIFY DATA TO LOAD
datapth = '';   % raw data folder not used

%% SESSIONS
% {loader, date, HMM folder for probe 1, HMM folder for probe 2}
% Where the two HMM folders on a row are the same, only one result exists for
% that session and both regions are aligned with it -- the consistency check
% below reports when that disagrees with the probe map.
hmmRoot = cfgPaths.hmmRoot;   % HMM-GLM disengagement times, one folder per session and probe
assert(exist(hmmRoot, 'dir') == 7, 'No disengagement-time folder: %s (set hmmRoot in setPaths.m)', hmmRoot);

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
    dataDirs{s,1} = fullfile(hmmRoot, sessionTable{s,3});
    dataDirs{s,2} = fullfile(hmmRoot, sessionTable{s,4});
end

%% PROBE MAP
% Which probe carries each region, per session. 0 = region not recorded that
% session. Position i refers to session i, so a length mismatch silently
% reassigns probes -- hence the assert.
% M1: 21 sessions / 6 animals; ALM: 15 sessions / 4 animals.
m1  = [1 1 0 1  2 0 2 1  1 2 2  2 1 1 2  0 1 2 2 1  2 2 2 1 0];
alm = [2 2 2 2  0 2 0 2  0 0 0  0 0 0 0  2 2 1 1 0  1 1 1 2 2];

assert(numel(m1) == nSessions && numel(alm) == nSessions, ...
    ['probe map has %d (m1) / %d (alm) entries but there are %d sessions. ' ...
     'Position i refers to session i, so a mismatch reassigns probes silently.'], ...
    numel(m1), numel(alm), nSessions);

regions = struct('name', {'M1','ALM'}, 'map', {m1, alm});

%% dataDir / probe consistency check
% The HMM folder names end in _P1 or _P2. For every (session, region) the
% figures actually use, check the suffix matches the probe the map asks for. A
% mismatch means that region is aligned with the OTHER probe's transition times.
logf('%s\n', repmat('-',1,72));
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
logf('%s\n\n', repmat('-',1,72));

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

% trials whose disengagement shift is NaN are dropped; they are tallied here
% and reported once after the loop instead of one warning per session.
nanDropR1 = 0;  nanDropR4 = 0;
%% MAIN LOOP
for sessnum = 1:length(all_meta)

    clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin

    meta = all_meta(sessnum,1);

    params.probe = {meta.probe};
    params.cluid = {};

    [obj, params, kin] = slimToLegacy(meta, params);

    fprintf('Session %d\n', sessnum);

    sessix = 1;

% Build allRegions from cluid
    [reg1, reg2, isSingleProbe] = regionSplit(obj(1), params(1));
    allRegions = {reg1, reg2};

    for aa = 1:2

% Only probes the maps use; unmapped probes (e.g. tjS1 on probe 1 in
% sessions 10 and 15) are skipped.
        if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end

        brainRegion = allRegions{aa};
        dataDir     = dataDirs{sessnum, aa};

        fileBases = {'R1_Trial_Track','R4_Trial_Track','R4_dt','R1_dt'};
        HMM = readHMM(dataDir, fileBases);

        r1Trials = table2array(HMM.R1_Trial_Track);
        r4Trials = table2array(HMM.R4_Trial_Track);

% dt files in 0.001s units -> neural bins
        r1DtsIdx = round(table2array(HMM.R1_dt) * 0.001 / params.dt);
        r4DtsIdx = round(table2array(HMM.R4_dt) * 0.001 / params.dt);
% Median transition time per session, logged as a check on the dt units.
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
            nanDropR1 = nanDropR1 + sum(bad1);
            trials1(bad1) = [];  dtBins1(bad1) = [];
        end

        bad4 = isnan(dtBins4);
        if any(bad4)
            nanDropR4 = nanDropR4 + sum(bad4);
            trials4(bad4) = [];  dtBins4(bad4) = [];
        end

% Shift every trial by its own disengagement bin, so t = 0 on obj.time is now
% the state transition rather than the go cue.
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

% ENGAGEMENT MODE WEIGHT VECTOR (from R4 trials)
        win1 = obj(sessix).time >= -cfg.engModeWin_s & obj(sessix).time < 0;
        win2 = obj(sessix).time >= 0    & obj(sessix).time <= cfg.engModeWin_s;
        w = engagementModeWeights(aligned4, trials4, win1, win2, cfg.sdFloor);

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

if nanDropR1 > 0 || nanDropR4 > 0
    fprintf('Dropped %d R1 and %d R4 trials with a NaN disengagement shift.\n', ...
        nanDropR1, nanDropR4);
end

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

% Colors: R1 dark blue, R4 cyan-blue
col_R1 = [0.2  0.2  0.7];   % dark blue  -> R1
col_R4 = [0.0  0.6  1.0];   % cyan-blue  -> R4

% Shared y limits
ylo = min([grandMean_R1 - ci_R1; grandMean_R4 - ci_R4]) - 0.03;
yhi = max([grandMean_R1 + ci_R1; grandMean_R4 + ci_R4]) + 0.03;
ylo = max(ylo, 0);
yhi = min(yhi, 1);

% Inset window; it also sets the grey patch that marks it on the main panel.
zoomWin  = [-0.100 0.100];   % s, relative to the state transition
xlimMain = [-1.300 0.450];   % s, the main panel

% Scale bar lengths; the labels are generated from these values.
xbar_left  = 0.250;   % s, main panel
xbar_right = 0.050;   % s, inset
ybar       = 0.1;     % a.u.

assert(xbar_left  < diff(xlimMain), ...
    'The %g ms scale bar does not fit in the main panel.', 1000*xbar_left);
assert(xbar_right < diff(zoomWin), ...
    'The %g ms scale bar does not fit in the inset.', 1000*xbar_right);
assert(zoomWin(1) >= xlimMain(1) && zoomWin(2) <= xlimMain(2), ...
    'The inset window is not inside the main panel, so the patch would be clipped.');

fprintf('Fig. 2M inset: %+.0f to %+.0f ms around the transition\n', ...
    1000*zoomWin(1), 1000*zoomWin(2));

colPatch = [0.93 0.93 0.93];
patchX   = [zoomWin(1) zoomWin(2) zoomWin(2) zoomWin(1)];
patchY   = [ylo ylo yhi yhi];

figure('Color','w','Units','normalized','Position',[0.1 0.2 0.7 0.55]);

% LEFT PANEL
ax1 = axes('Position',[0.06 0.15 0.52 0.75]);
hold(ax1,'on');

% the stretch the inset enlarges, drawn first so it sits behind the traces
patch(ax1, patchX, patchY, colPatch, 'EdgeColor','none', 'FaceAlpha',1.0);

fill(ax1, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R1 - ci_R1; flipud(grandMean_R1 + ci_R1)], ...
     col_R1, 'FaceAlpha',0.25, 'EdgeColor','none');
fill(ax1, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R4 - ci_R4; flipud(grandMean_R4 + ci_R4)], ...
     col_R4, 'FaceAlpha',0.25, 'EdgeColor','none');

plot(ax1, tAxis_dis, grandMean_R1, '-', 'Color',col_R1, 'LineWidth',2.5);
plot(ax1, tAxis_dis, grandMean_R4, '-', 'Color',col_R4, 'LineWidth',2.5);

xline(ax1, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

xlim(ax1, xlimMain);
ylim(ax1, [ylo yhi]);
axis(ax1,'off');

% X scale bar
xbar_x = xlimMain(1) + 0.05;
xbar_y = ylo + 0.01;
plot(ax1, [xbar_x, xbar_x + xbar_left], [xbar_y xbar_y], 'k-', 'LineWidth',2);
text(ax1, xbar_x + xbar_left/2, xbar_y - 0.04, sprintf('%g ms', 1000*xbar_left), ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

% Y scale bar
ybar_x = xbar_x - 0.04;
ybar_y = xbar_y + 0.01;
plot(ax1, [ybar_x ybar_x], [ybar_y, ybar_y + ybar], 'k-', 'LineWidth',2);
text(ax1, ybar_x - 0.06, ybar_y + ybar/2, sprintf('%.1f\n(a.u.)', ybar), ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

% State transition label
text(ax1, 0, ylo - 0.08, 'State transition', ...
    'HorizontalAlignment','center','FontSize',10,'Color','k');

hold(ax1,'off');

% RIGHT PANEL -- the grey stretch above, enlarged
ax2 = axes('Position',[0.62 0.15 0.32 0.75]);
hold(ax2,'on');

patch(ax2, patchX, patchY, colPatch, 'EdgeColor','none', 'FaceAlpha',1.0);

fill(ax2, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R1 - ci_R1; flipud(grandMean_R1 + ci_R1)], ...
     col_R1, 'FaceAlpha',0.25, 'EdgeColor','none');
fill(ax2, [tAxis_dis(:); flipud(tAxis_dis(:))], ...
     [grandMean_R4 - ci_R4; flipud(grandMean_R4 + ci_R4)], ...
     col_R4, 'FaceAlpha',0.25, 'EdgeColor','none');

plot(ax2, tAxis_dis, grandMean_R1, '-', 'Color',col_R1, 'LineWidth',2.5);
plot(ax2, tAxis_dis, grandMean_R4, '-', 'Color',col_R4, 'LineWidth',2.5);

xline(ax2, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

xlim(ax2, zoomWin);
ylim(ax2, [ylo yhi]);
axis(ax2,'off');

% X scale bar
xbar_x2 = zoomWin(1) + 0.01;
xbar_y2 = ylo + 0.01;
plot(ax2, [xbar_x2, xbar_x2 + xbar_right], [xbar_y2 xbar_y2], 'k-', 'LineWidth',2);
text(ax2, xbar_x2 + xbar_right/2, xbar_y2 - 0.04, sprintf('%g ms', 1000*xbar_right), ...
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

%% LOCAL FUNCTIONS

function [reg1, reg2, isSingle] = regionSplit(obj, params)
% Cluster indices for the two probes. A single-probe session returns an empty
% second region and isSingle = true.
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
% Read the HMM-GLM result tables (one CSV per name) into a struct.
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
%   p(t) = sum_i w_i * r_i(t),  with sum_i |w_i| = 1
% 'mean' divides by the unit count, which rescales each session by its own N.
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


function w = engagementModeWeights(aligned4, trials4, win1, win2, sdFloor)
% Weights defining the engagement mode, from the C4-reward trials.
%
% For each neuron: the difference between its mean firing rate in the window
% before the state transition and in the window after it, divided by the pooled
% across-trial standard deviation of those two windows. That is d', so a neuron
% counts for more only when it separates the two states reliably, rather than
% merely by a large number of spikes. Weights are normalized so sum|w| = 1.
%
% Each trial contributes one value per window (the bins inside the window are
% averaged first), so the standard deviations are taken across trials.

    N   = size(aligned4, 2);
    nTr = numel(trials4);

    preTr  = reshape(nanmean(aligned4(win1,:,trials4), 1), N, nTr);
    postTr = reshape(nanmean(aligned4(win2,:,trials4), 1), N, nTr);

% pooled SD: square each window's across-trial SD, average the two, square root
    sdPooled = sqrt((nanstd(preTr,0,2).^2 + nanstd(postTr,0,2).^2) / 2)';

    w = (nanmean(preTr,2)' - nanmean(postTr,2)') ./ max(sdPooled, sdFloor);

    w(~isfinite(w)) = 0;
    w = w / sum(abs(w));
end

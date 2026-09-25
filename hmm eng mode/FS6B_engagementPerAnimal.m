%% FS6B_engagementPerAnimal.m
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
%    params.smooth      0
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -2 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


sz = 26;

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

%% PARAMETERS
params.alignEvent          = 'goCue';
params.behav_only          = 0;
params.timeWarp            = 0;
params.nLicks              = 20;
params.lowFR               = 0.01;   % minimum mean firing rate, Hz

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1'};

params.tmin   = -2;
params.tmax   = 5;
params.dt     = 1/100;
params.smooth = 0;
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

y1_all = []; y2_all = []; y3_all = []; y4_all = [];
allP1 = {}; allP4 = {};
allmoveP1 = []; allmoveP4 = [];
sess_y = {}; sess_y4 = {}; sess_yhat = {}; sess_yhat4 = {};
r1time = {};
r4time = {};

params.behav_only = 0;

for sessnum = 1:length(all_meta)

clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin

meta = all_meta(sessnum,1);

params.probe = {meta.probe};
params.cluid = {};

[obj, params, kin] = slimToLegacy(meta, params);

% define region indices from the two probes
[reg1, reg2, isSingleProbe] = regionSplit(obj(1), params(1));   % was: params.cluid{1,1} with no guard
allRegions = {reg1, reg2};

for sessix = 1:numel(meta)
end

nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
end

conds2use = [1];
kinfeat   = 'tongue_length';
sessix    = 1;

kinix  = find(strcmp(kin(sessix).featLeg, kinfeat));
Length = kin.dat(:, :, kinix);

for aa = 1:2

% Only probes the maps use (same rule as tongue_r14's probesOf). aa is the
% PROBE index; m1/alm give the probe each region is on.
    if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end

    brainRegion = allRegions{aa};
    dataDir     = dataDirs{sessnum, aa};

    fileBases = {'R1_Trial_Track','R4_Trial_Track','R4_dt','R1_dt'};
    HMM = readHMM(dataDir, fileBases);   % was: assignin into the base workspace

    r1Trials = table2array(HMM.R1_Trial_Track);
    r4Trials = table2array(HMM.R4_Trial_Track);

% R1_dt and R4_dt are indices at 0.001s resolution (1000Hz)
% convert to seconds: * 0.001
% convert to bins at 100Hz: * 100
% net: * 0.001 * 100 = * 0.1
    r1DtsSec = table2array(HMM.R1_dt) * 0.001;   % time in seconds from go cue to event
    r4DtsSec = table2array(HMM.R4_dt) * 0.001;
% bins at THIS script's sampling rate, not a hard-coded 100 Hz.
    r1DtsIdx = round(r1DtsSec / params.dt);
    r4DtsIdx = round(r4DtsSec / params.dt);
% HMM dt files: ms from the go cue (confirmed against Fig. 2K / S6B via the
% bigPlot session log). Printed so every run shows the same check.
    logf('  [hmm] sess %2d probe %d | median transition re: GC  R1 %.3f s  R4 %.3f s\n', ...
        sessnum, aa, median(r1DtsIdx,'omitnan')*params.dt, median(r4DtsIdx,'omitnan')*params.dt);

    neurons = brainRegion;
    allDat  = obj.trialdat(:,neurons,:);
    T       = size(allDat,1);
    N       = size(allDat,2);

    trials1 = r1Trials(:);
    dtBins1 = r1DtsIdx(:);
    dtSecs1 = r1DtsSec(:);   % keep seconds version for r1time
    trials4 = r4Trials(:);
    dtBins4 = r4DtsIdx(:);
    dtSecs4 = r4DtsSec(:);   % keep seconds version for r4time

    bad1 = isnan(dtBins1);
    if any(bad1)
        warning('%d R1-trials have NaN shift -> removing them', sum(bad1));
        trials1(bad1) = [];
        dtBins1(bad1) = [];
        dtSecs1(bad1) = [];
    end

    bad4 = isnan(dtBins4);
    if any(bad4)
        warning('%d R4-trials have NaN shift -> removing them', sum(bad4));
        trials4(bad4) = [];
        dtBins4(bad4) = [];
        dtSecs4(bad4) = [];
    end

% --- subtract first lick latency to get time from first contact ---
% dtSecs is currently time from go cue to disengagement event
% subtract first post-GC lick latency → time from first contact to event

    dtFromFirstLick1 = NaN(size(dtSecs1));
    for ii = 1:numel(trials1)
        tr        = trials1(ii);
        gocueTime = obj.bp.ev.goCue(tr);
        lickTimes = obj.bp.ev.lickL{tr};
        idx = find(lickTimes > gocueTime, 1, 'first');
        if ~isempty(idx)
            lickLatency          = lickTimes(idx) - gocueTime;
            dtFromFirstLick1(ii) = dtSecs1(ii) - lickLatency;
        end
    end

    dtFromFirstLick4 = NaN(size(dtSecs4));
    for ii = 1:numel(trials4)
        tr        = trials4(ii);
        gocueTime = obj.bp.ev.goCue(tr);
        lickTimes = obj.bp.ev.lickL{tr};
        idx = find(lickTimes > gocueTime, 1, 'first');
        if ~isempty(idx)
            lickLatency          = lickTimes(idx) - gocueTime;
            dtFromFirstLick4(ii) = dtSecs4(ii) - lickLatency;
        end
    end

    aligned1 = doAlign(allDat, trials1, dtBins1);
    aligned4 = doAlign(allDat, trials4, dtBins4);

    nT1 = numel(trials1);
    nT4 = numel(trials4);

    psth1 = nan(T,N);
    psth4 = nan(T,N);
    for nn = 1:N
        psth1(:,nn) = nanmean(aligned1(:,nn,trials1),3);
        psth4(:,nn) = nanmean(aligned4(:,nn,trials4),3);
    end

    win1 = obj.time >= -cfg.engModeWin_s & obj.time < 0;
    win2 = obj.time >= 0    & obj.time <= cfg.engModeWin_s;
    m1_4 = mean(psth4(win1,:), 1);
    m2_4 = mean(psth4(win2,:), 1);
    w = m1_4 - m2_4;
    w = w / sum(abs(w));

    proj1_goCue = nan(T, nT1);
    proj4_goCue = nan(T, nT4);
    for k = 1:nT1
        tr = trials1(k);
        X  = allDat(:,:,tr);
        proj1_goCue(:,k) = projectMode(X, w, cfg.projMode);
    end
    for k = 1:nT4
        tr = trials4(k);
        X  = allDat(:,:,tr);
        proj4_goCue(:,k) = projectMode(X, w, cfg.projMode);
    end

    mod1 = mean(proj1_goCue, 2);
    mod4 = mean(proj4_goCue, 2);

    All_R1Modes{sessnum, aa} = mod1;
    All_R4Modes{sessnum, aa} = mod4;
    AllProj_R1{sessnum, aa}  = proj1_goCue;
    AllProj_R4{sessnum, aa}  = proj4_goCue;
    AllTime{sessnum, aa}     = obj.time;

% store time from first contact in seconds
    r1time{sessnum, aa} = dtFromFirstLick1;
    r4time{sessnum, aa} = dtFromFirstLick4;

end   % aa loop
end   % sessnum loop


rowLengths = [4, 4, 3, 4, 5, 5];
% The session list is grouped by animal in this order, so rowLengths has to add
% up to the number of sessions. Nothing checked that, and a session added to
% all_meta without touching this line would have silently reassigned every
% animal boundary after it.
assert(sum(rowLengths) == numel(all_meta), ...
    ['rowLengths sums to %d but all_meta has %d sessions. These vectors index ' ...
     'sessions by position, so a mismatch reassigns animals silently.'], ...
    sum(rowLengths), numel(all_meta));
cumSess    = [0 cumsum(rowLengths)];
nAnimals   = numel(rowLengths);
alpha      = 0.05;

meanR1 = nan(nAnimals, 1);
ciR1   = nan(nAnimals, 1);
meanR4 = nan(nAnimals, 1);
ciR4   = nan(nAnimals, 1);
nSessPerAnimal = nan(nAnimals, 1);

for a = 1:nAnimals
    sessStart = cumSess(a) + 1;
    sessEnd   = cumSess(a+1);
    sessR1 = [];   % ONE VALUE PER SESSION PER PROBE, not a pool of every trial
    sessR4 = [];

    for j = sessStart:sessEnd
        code1 = m1(j);
        code2 = alm(j);

        if code1 > 0
            v1 = r1time{j, code1}(:);
            v4 = r4time{j, code1}(:);
            v1(isnan(v1)) = [];
            v4(isnan(v4)) = [];
            sessR1(end+1) = mean(v1);   %#ok<AGROW>
            sessR4(end+1) = mean(v4);   %#ok<AGROW>
        end
        if code2 > 0
            v1 = r1time{j, code2}(:);
            v4 = r4time{j, code2}(:);
            v1(isnan(v1)) = [];
            v4(isnan(v4)) = [];
            sessR1(end+1) = mean(v1);   %#ok<AGROW>
            sessR4(end+1) = mean(v4);   %#ok<AGROW>
        end
    end

% The animal's value is the mean over ITS SESSIONS, and the CI is over
% trials as independent replicates, making the interval far too narrow.
% The unit of analysis for this panel is the animal.
    meanR1(a) = mean(sessR1);
    meanR4(a) = mean(sessR4);
    n1 = numel(sessR1);   n4 = numel(sessR4);
    nSessPerAnimal(a) = n1;
    ciR1(a) = tinv(1-alpha/2, max(n1-1,1)) * std(sessR1) / sqrt(max(n1,1));
    ciR4(a) = tinv(1-alpha/2, max(n4-1,1)) * std(sessR4) / sqrt(max(n4,1));
end

% --- print summary ---
for a = 1:nAnimals
    fprintf('Animal %d: meanR1=%.4fs  meanR4=%.4fs\n', a, meanR1(a), meanR4(a));
end
fprintf('Grand mean R1=%.4fs  R4=%.4fs\n', mean(meanR1), mean(meanR4));

%% Fig. S6B : THE TEST THE PANEL REPORTS
% Legend: "Disengagement times were significantly later for C4-reward trials
% compared to C1-reward trials (asterisk; p<0.05, one-tailed Wilcoxon
% signed-rank test)."
% Paired, because the same animal contributes both values. One-tailed for
% C4 > C1, the direction the Results state in advance. One test, so there is
% nothing to correct for.
% THIS DID NOT EXIST BEFORE. The panel drew a bracket and a literal '*' with no
% test anywhere in the file.
okA = isfinite(meanR1) & isfinite(meanR4);
nA  = sum(okA);

fprintf('\n%s\n', repmat('=',1,74));
fprintf('Fig. S6B | disengagement time re: first contact, C1 vs C4, paired by animal\n');
fprintf('%s\n', repmat('-',1,74));
fprintf('%-8s %-10s %-12s %-12s %s\n', 'animal', 'sessions', 'C1 (s)', 'C4 (s)', 'C4 - C1 (s)');
for a = 1:nAnimals
    fprintf('%-8d %-10d %-12.4f %-12.4f %+0.4f\n', a, nSessPerAnimal(a), ...
        meanR1(a), meanR4(a), meanR4(a)-meanR1(a));
end
fprintf('%-8s %-10d %-12.4f %-12.4f %+0.4f\n', 'MEAN', sum(nSessPerAnimal(okA)), ...
    mean(meanR1(okA)), mean(meanR4(okA)), mean(meanR4(okA)-meanR1(okA)));
fprintf('%s\n', repmat('-',1,74));

pSR = NaN;
if nA >= 2
    pSR  = signrank(meanR4(okA), meanR1(okA), 'tail', 'right');   % C4 later than C1
    pSR2 = signrank(meanR4(okA), meanR1(okA));
    fprintf('Wilcoxon signed-rank, one-tailed (C4 > C1) : p = %.4f%s\n', ...
        pSR, repmat('   *', 1, pSR < alpha));
    fprintf('two-tailed reference                      : p = %.4f\n', pSR2);
    fprintf('exact floor at n = %d pairs                : %.4f one-tailed\n', nA, 2^(-nA));
    if 2^(-nA) >= alpha
        fprintf(['!! At n = %d the smallest one-tailed p an exact signed-rank can return is\n' ...
                 '   %.4f, which is not below alpha. No effect size can produce a star at\n' ...
                 '   this n with this test.\n'], nA, 2^(-nA));
    end
else
    fprintf('test not run: only %d animals with finite values.\n', nA);
end
fprintf('%s\n', repmat('=',1,74));

figure('Color','w', 'Units','inches', 'Position',[1 1 2.8 4.5]);
hold on

% --- Individual animal lines (gray) ---
for a = 1:nAnimals
    plot([1 2], [meanR1(a) meanR4(a)], '-', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5, ...
        'HandleVisibility','off');
end

% --- Grand mean line (thick black) ---
grandMeanR1 = mean(meanR1);
grandMeanR4 = mean(meanR4);
plot([1 2], [grandMeanR1, grandMeanR4], '-k', ...
    'LineWidth', 4.5, ...
    'HandleVisibility','off');

% --- Significance bracket ---
yBracket   = 0.78;
bracketPad = 0.02;
line([1 2], [yBracket yBracket],           'Color','k', 'LineWidth',1.2);
line([1 1], [yBracket-bracketPad yBracket],'Color','k', 'LineWidth',1.2);
line([2 2], [yBracket-bracketPad yBracket],'Color','k', 'LineWidth',1.2);
% unconditionally, with no test anywhere in the file.
if isfinite(pSR) && pSR < alpha
    text(1.5, yBracket+0.015, '*', 'HorizontalAlignment','center', ...
        'FontSize', 16, 'FontWeight','bold', 'Color','k');
else
    text(1.5, yBracket+0.015, 'n.s.', 'HorizontalAlignment','center', ...
        'FontSize', 10, 'Color',[0.35 0.35 0.35]);
end

% --- Axes formatting ---
xlim([0.5 2.5]);
ylim([0 0.8]);
yticks(0:0.2:0.8);
xticks([1 2]);
xticklabels({});
ylabel('Time from First Contact (s)', 'FontSize',11, 'FontWeight','normal');
set(gca, 'FontSize',10, 'LineWidth',1.2, 'Box','off', ...
    'TickDir','out', 'XColor','k', 'YColor','k');

% --- Colored x-axis labels (stacked, manual) ---
text(1, -0.07,  'C1',     'Units','data', ...
    'HorizontalAlignment','center', 'FontSize',11, ...
    'FontWeight','bold', 'Color',[0 0.2 0.9]);
text(1, -0.115, 'Trials', 'Units','data', ...
    'HorizontalAlignment','center', 'FontSize',11, ...
    'FontWeight','bold', 'Color',[0 0.2 0.9]);

text(2, -0.07,  'C4',     'Units','data', ...
    'HorizontalAlignment','center', 'FontSize',11, ...
    'FontWeight','bold', 'Color',[0 0.8 1]);
text(2, -0.115, 'Trials', 'Units','data', ...
    'HorizontalAlignment','center', 'FontSize',11, ...
    'FontWeight','bold', 'Color',[0 0.8 1]);

hold off

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

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

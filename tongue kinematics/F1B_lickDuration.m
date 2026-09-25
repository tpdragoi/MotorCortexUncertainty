%% F1B_lickDuration.m
%  Duration of successive tongue protrusions, Simple Reward Task.
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

clear, clc

% ---- SELF-CONTAINED: data and functions come only from this folder ----
% Data\R1 holds the exported <ANM>_<DATE>_obj.mat / _kin.mat; shared\ holds the
% pipeline functions it needs (shared\pipelineCopies). uninstructedMovements_v2-main
% and the raw data tree are not on the path.
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'loadBehavSession.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));
dataDir = fullfile(v2Root, 'Data', 'R1');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s', dataDir);
% loader name -> animal name as stored in the exported files (obj.pth.anm)
anmOfLoader = containers.Map({'loadTD10s_many', 'loadTD9s_many', 'loadTD27_many', 'loadTD26_many', 'loadTD3l_many', 'loadTD2l_many', 'loadTD4l_many', 'loadTD5l_many', 'loadTD1_many', 'loadTD4_many', 'loadTD13_many', 'loadTD15_many', 'loadTD8_many', 'loadTD22_many', 'loadTD23_many', 'loadYH2_many', 'loadYH1_many', 'loadTD4f_many', 'loadTD7f_many', 'loadTD24_many', 'loadTD25_many', 'loadTDv1_many', 'loadTDv4_many', 'loadTDv5_many', 'loadTDv6_many'}, ...
                             {'TD10si', 'TD9si', 'TD27d', 'TD26d', 'TDl3', 'TDl2', 'TDl4', 'TD5l', 'TD1d', 'TD4d', 'TD13d', 'TD15d', 'TD8d', 'TD22d', 'TD23d', 'YH2', 'YH1', 'TD4f', 'TD7f', 'TD24d', 'TD25d', 'TDv1', 'TDv4', 'TDv5', 'TDv6'});

%  LICK BOUT KINEMATICS -- R1 (reward at contact 1)
%  Per-contact lick duration, maximum tongue length, inter-contact interval
%  and within-contact shape variance, from the first lick bout after the go
%  cue. Behaviour only (params.behav_only = 1); no neural data is loaded.
%  WHAT CHANGED FROM THE PREVIOUS VERSION
%  1. SESSION LIST IS NOW A TABLE, NOT 36 meta VARIABLES.
%     The old m1/alm maps had 35 entries indexed by LOADER-CALL ORDER, but
%     rows. Indexing the map by sessnum would therefore have read the wrong
%     session's probe from session 7 onward. sessionTable below carries the
%     loader, the date and both probe numbers on one row, so the map cannot
%     drift out of step with the data no matter which rows are switched off.
%  2. REGION FILTER. Sessions with no M1 and no ALM probe are skipped when
%     CFG.regionFilter is true. (On the current active list this removes
%     nothing -- all 22 have at least one region -- but it matters as soon
%     as more rows are switched on.)
%  3. CONTACT 1 IS THE FIRST DETECTED CONTACT. The old code indexed columns
%     by position in the post-go-cue lick-event list, so if the first lick
%     event produced no detectable tongue bout, column 1 stayed NaN and
%     "contact 1" silently became a different lick on different trials.
%  4. FIRST-BOUT TRIAL FILTER (new). A trial is used only if its first bout
%     contains at least CFG.boutMinLicks contacts spaced under
%     CFG.boutMaxILI_s. The bout ends at the last contact not followed by
%     another within CFG.boutEndGap_s.
%  5. t-tests -> two-tailed Wilcoxon signed-rank across sessions, BH-corrected.
%  7. Only contact duration and its trial-to-trial SD are plotted, one
%     the workspace as allVals_len and allVals_ili.
%  6. Bugs fixed: the empty catch that hid every per-trial failure now
%     counts them; the duration filter used a hardcoded 0.005 instead of
%     params.dt; the ILI row conflated "next contact" with "next lick
%     number"; three of the four plots were labelled "Mean Lick Duration"
%  loadTD27_many 2025-07-27 is listed as m1 = 1 AND alm = 1. One probe
%  cannot be in both regions. The row is flagged at load time and its
%  region assignment should be treated as unknown until corrected.

%% OPTIONS

CFG.regionFilter  = true;   % skip sessions with neither M1 nor ALM
% CFG.contactOrigin removed: lick columns are the n-th post-go-cue lick,

% ---- first-bout trial filter ----
% The bout starts at the first contact after the go cue. It ends at the last
% contact that is NOT followed by another contact within CFG.boutEndGap_s.
% The trial is kept only if the leading run of contacts spaced at most
% CFG.boutMaxILI_s apart is at least CFG.boutMinLicks long.
CFG.applyBoutFilter = true;
CFG.boutMinLicks    = 8;
CFG.boutMaxILI_s    = 0.20;
CFG.boutEndGap_s    = 0.30;

CFG.maxCols   = 35;   % columns held per trial before trimming
CFG.nLicks    = 15;   % contacts carried into the per-session summary
CFG.showLicks = 1:8;   % contacts plotted, and the BH family

% ---- bout detection: SAME values and SAME code as VTA/R14/R16/Learning ----
cfg.nanBridge_sec   = 0.020;   % interpolate tracking gaps up to this long
cfg.minBout_sec     = 0.030;   % a contact run shorter than this is not a lick
cfg.maxBout_sec     = 0.400;   % artefact guard after the contact-driven split
cfg.contactTol_sec  = 0.020;   % a run must fall within this of a lick event
cfg.maxLickIdx      = CFG.maxCols;

% ---- statistics ----
STATS.test  = 'signrank';   % 'signrank' (paired, across sessions) | 'ttest'
STATS.tail  = 'left';   % contact N vs contact 1, ONE-tailed: shorter than contact 1 (as R16)
STATS.alpha = 0.05;
STATS.bh    = true;

PLOT.fontSize = 12;

%% PATHS

%% PARAMETERS

params.alignEvent = 'goCue';
params.behav_only = 1;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz

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

params.tmin   = -2.5;
params.tmax   = 4;
params.dt     = 1/200;
params.smooth = 1;

params.quality = {'good'};   % good units only (findClusters trims blanks and ignores case)

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};
params.feat_varToExplain = 80;
params.N_varToExplain    = 80;
params.advance_movement  = 0;

params.fcut   = 10;
params.cond   = 5;
params.method = 'xcorr';
params.fa     = false;
params.bctype = 'reflect';

%% SESSIONS AND PROBE MAP
% Columns: loader | date | M1 probe | ALM probe | use
% This table replaces meta1..meta36: the probe map travels with the row, so
% switching a session off cannot shift anyone else's region assignment.
% Values decoded from the m1/alm arrays in the previous version, which were
% carries 0/0 here, which is what confirms the alignment.

sessionTable = {
    'loadTD10s_many', '2024-07-09', 0, 2, true
    'loadTD10s_many', '2024-07-10', 1, 2, true
    'loadTD10s_many', '2024-07-11', 0, 2, true
    'loadTD10s_many', '2024-07-13', 2, 0, true
    'loadTD10s_many', '2024-07-14', 2, 1, true

    'loadTD9s_many',  '2024-07-05', 2, 1, true
    'loadTD9s_many',  '2024-07-06', 0, 0, false
    'loadTD9s_many',  '2024-07-07', 2, 0, true
    'loadTD9s_many',  '2024-07-08', 0, 1, true
    'loadTD9s_many',  '2024-07-09', 1, 0, true
    'loadTD9s_many',  '2024-07-10', 0, 1, true

    'loadTD3l_many',  '2025-02-18', 0, 0, false
    'loadTD3l_many',  '2025-02-19', 0, 0, false
    'loadTD3l_many',  '2025-02-20', 0, 0, false
    'loadTD3l_many',  '2025-02-22', 0, 0, false

    'loadTD2l_many',  '2025-02-18', 0, 0, false
    'loadTD2l_many',  '2025-02-19', 0, 0, false
    'loadTD2l_many',  '2025-02-20', 0, 0, false
    'loadTD2l_many',  '2025-02-21', 0, 0, false
    'loadTD2l_many',  '2025-02-24', 0, 0, false

    'loadTD27_many',  '2025-07-25', 2, 0, true
    'loadTD27_many',  '2025-07-26', 1, 0, true
    'loadTD27_many',  '2025-07-27', 1, 1, true   % <-- CONFLICT, see header
    'loadTD27_many',  '2025-07-28', 0, 1, true
    'loadTD27_many',  '2025-07-29', 2, 1, true
    'loadTD27_many',  '2025-07-30', 0, 0, false
    'loadTD27_many',  '2025-07-31', 1, 0, true

    'loadTD26_many',  '2025-07-30', 1, 2, true
    'loadTD26_many',  '2025-07-31', 0, 0, false
    'loadTD26_many',  '2025-08-01', 2, 0, true
    'loadTD26_many',  '2025-08-02', 0, 0, false
    'loadTD26_many',  '2025-08-03', 1, 0, true
    'loadTD26_many',  '2025-08-05', 2, 1, true
    'loadTD26_many',  '2025-08-06', 2, 1, true
    'loadTD26_many',  '2025-08-07', 1, 0, true
};

% ---- integrity check: one probe cannot serve two regions ----
for r = 1:size(sessionTable,1)
    m1r  = sessionTable{r,3};
    almr = sessionTable{r,4};
    if m1r > 0 && almr > 0 && m1r == almr
        warning('SessionTable:probeConflict', ...
            'Row %d (%s %s): M1 and ALM are both probe %d. Region assignment is unreliable for this session.', ...
            r, sessionTable{r,1}, sessionTable{r,2}, m1r);
    end
end

% ---- build all_meta, carrying the probe map alongside ----
all_meta  = [];
sessProbe = zeros(0,2);
sessName  = {};
nSkipReg  = 0;

for r = 1:size(sessionTable,1)
    if ~sessionTable{r,5}, continue; end
    if CFG.regionFilter && sessionTable{r,3} == 0 && sessionTable{r,4} == 0
        nSkipReg = nSkipReg + 1;
        continue
    end
    m = struct('anm', anmOfLoader(sessionTable{r,1}), 'date', sessionTable{r,2});   % was feval(loader)
    all_meta  = [all_meta; m];
    sessProbe = [sessProbe; sessionTable{r,3}, sessionTable{r,4}];
    sessName{end+1} = sprintf('%s %s', sessionTable{r,1}, sessionTable{r,2});   %#ok<SAGROW>
end

nSessTot = size(all_meta,1);
assert(nSessTot == size(sessProbe,1), 'probe map and all_meta are out of step');


%% PER-SESSION EXTRACTION

% trials after these indices are dropped in the named sessions
exclusions = {
    'TD13d', '2024-11-11', 278
    'TD8d',  '2024-09-07', 313
    'TD8d',  '2024-09-09', 298
    'TD26d', '2025-08-07', 187
};

allVals_dur = cell(1, nSessTot);
allVals_len = cell(1, nSessTot);
allVals_ili = cell(1, nSessTot);
allVals_sd  = cell(1, nSessTot);

nTrialKept = zeros(1, nSessTot);
nTrialSeen = zeros(1, nSessTot);
nTrialBout = zeros(1, nSessTot);
nTrialFail = zeros(1, nSessTot);

for sessnum = 1:nSessTot

    clear obj kin me Length

    meta         = all_meta(sessnum, 1);
    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);
    trialSet      = (1:obj.bp.Ntrials)';
    fprintf('Session %d\n', sessnum);

% ---- trial exclusions ----
    condtrix = trialSet;
    cutAt    = Inf;
    for ex = 1:size(exclusions,1)
        if strcmp(obj.pth.anm, exclusions{ex,1}) && strcmp(obj.pth.dt, exclusions{ex,2})
            cutAt = exclusions{ex,3};
            break
        end
    end
    condtrix(condtrix > cutAt) = [];

% Length's columns are condtrix in order. Because the exclusion only ever
% truncates the tail, column j still corresponds to trial j; do not change
% the exclusion to remove interior trials without reindexing here.
    kinix  = find(strcmp(kin.featLeg, 'tongue_length'));
    Length = kin.dat(:, condtrix, kinix);

% ---- R1 trials ----
    allTr = 1:obj.bp.Ntrials;
    hitTr = allTr(obj.bp.hit == 1);
    P8    = intersect(hitTr, allTr(obj.bp.rewardedLick == 1))';
    P8(P8 > cutAt) = [];

    dur_all = [];  len_all = [];  ili_all = [];

    for i = 1:numel(P8)

        trial = P8(i);
        nTrialSeen(sessnum) = nTrialSeen(sessnum) + 1;

        try
            val       = Length(:, trial);
            licks_abs = obj.bp.ev.lickL{trial, 1};
            gc        = obj.bp.ev.goCue(trial);

            if numel(licks_abs) <= 2 || all(isnan(val)), continue; end

            licks = licks_abs - gc;   % contact times relative to the go cue
            licks = licks(licks > 0);
            licks = sort(licks(:))';
            if isempty(licks), continue; end

% ---------- FIRST-BOUT FILTER ----------
% bout ends at the last contact not followed by one within
% CFG.boutEndGap_s; the run used for the criterion is the leading
% stretch of contacts spaced at most CFG.boutMaxILI_s apart.
            iliAll = diff(licks);

            kEnd = find(iliAll > CFG.boutEndGap_s, 1, 'first');
            if isempty(kEnd), boutN = numel(licks); else, boutN = kEnd; end

            kLoose = find(iliAll > CFG.boutMaxILI_s, 1, 'first');
            if isempty(kLoose), tightN = numel(licks); else, tightN = kLoose; end

            if CFG.applyBoutFilter && tightN < CFG.boutMinLicks, continue; end
            nTrialBout(sessnum) = nTrialBout(sessnum) + 1;

% ---------- shared extraction (lickFeaturesForTrial) ----------
            [dur_row, len_row, ili_row, why] = lickFeaturesForTrial( ...
                val, obj.time, licks_abs, gc, params.dt, cfg);
            if why > 0 && why < 5, continue; end

            dur_all = [dur_all; dur_row];
            len_all = [len_all; len_row];
            ili_all = [ili_all; ili_row];
            nTrialKept(sessnum) = nTrialKept(sessnum) + 1;

        catch ME
            nTrialFail(sessnum) = nTrialFail(sessnum) + 1;
            if nTrialFail(sessnum) == 1
                fprintf('     first trial error in this session: %s\n', ME.message);
            end
        end
    end

% ---- per-session summary over the first CFG.nLicks contacts ----
    nL = CFG.nLicks;
    if isempty(dur_all)
        allVals_dur{sessnum} = nan(nL,1);
        allVals_len{sessnum} = nan(nL,1);
        allVals_ili{sessnum} = nan(nL,1);
        allVals_sd{sessnum}  = nan(nL,1);
    else
        allVals_dur{sessnum} = nanmean(dur_all(:,1:nL), 1)';
        allVals_len{sessnum} = nanmean(len_all(:,1:nL), 1)';
        allVals_ili{sessnum} = nanmean(ili_all(:,1:nL), 1)';
        allVals_sd{sessnum}  = nanstd(dur_all(:,1:nL), [], 1)';
    end
end


%% FIGURES
% One figure per measure, in the original style: a shaded 95% CI box at each
% contact with a filled marker at the mean, in a tall narrow window.
% allVals_len and allVals_ili are still computed above and are available in
% the workspace; they are simply not plotted.

licks = CFG.showLicks;
nLk   = numel(licks);

panels = {allVals_dur, allVals_sd};
pName  = {'Contact duration', 'Trial-to-trial SD of contact duration'};
pYlab  = {'Mean Lick Duration (s)', 'SD of Lick Duration (s)'};
pCol   = {'r', 'r'};

barWidth = 0.2;

for pIdx = 1:numel(panels)

    M = cell2mat(panels{pIdx});   % nLicks x nSessions
    M = M(licks, :);
    n = sum(~isnan(M), 2)';

    mu  = nanmean(M, 2)';
    sem = nanstd(M, [], 2)' ./ sqrt(max(n,1));
    tc  = 1.96*ones(size(n));
    tc(n > 1) = tinv(1 - STATS.alpha/2, n(n > 1) - 1);
    ci  = sem .* tc;

% ---- contact N vs contact 1, paired across sessions ----
    pRaw = nan(1, nLk);
    for ii = 2:nLk
        aV = M(ii,:);  bV = M(1,:);
        ok = ~isnan(aV) & ~isnan(bV);
        if sum(ok) < 2, continue; end
        if strcmp(STATS.test, 'signrank')
            if any(aV(ok) - bV(ok) ~= 0)
                pRaw(ii) = signrank(aV(ok), bV(ok), 'tail', STATS.tail);
            end
        else
            [~, pRaw(ii)] = ttest(aV(ok), bV(ok), 'Tail', STATS.tail);
        end
    end

% ---- Benjamini-Hochberg over the contacts tested ----
    pAdj = nan(1, nLk);
    okp  = ~isnan(pRaw);
    if any(okp) && STATS.bh
        pv = pRaw(okp);
        [ps, ord] = sort(pv(:));
        mm  = numel(ps);
        ad  = ps .* (mm ./ (1:mm)');
        ad  = flipud(cummin(flipud(ad)));
        tmp = nan(mm,1); tmp(ord) = min(ad,1);
        pAdj(okp) = tmp';
    elseif ~STATS.bh
        pAdj = pRaw;
    end

% ---- draw ----
    figure; hold on;
    ax = gca;

    for ii = 1:nLk
        xShade = [licks(ii)-barWidth, licks(ii)+barWidth, ...
                  licks(ii)+barWidth, licks(ii)-barWidth];
        yShade = [mu(ii)-ci(ii), mu(ii)-ci(ii), mu(ii)+ci(ii), mu(ii)+ci(ii)];
        fill(xShade, yShade, pCol{pIdx}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(licks(ii), mu(ii), 'o', 'Color', pCol{pIdx}, ...
             'MarkerFaceColor', pCol{pIdx}, 'MarkerSize', 8);
    end

    yl = ylim(ax);  ylim(ax, [yl(1), yl(2) + 0.12*diff(yl)]);  yl = ylim(ax);
    for ii = 1:nLk
        if pAdj(ii) < STATS.alpha
            text(licks(ii), yl(2) - 0.045*diff(yl), '*', 'Color', 'k', ...
                 'FontSize', PLOT.fontSize+4, 'HorizontalAlignment', 'center');
        end
    end

    xlabel('Lick Number', 'FontSize', PLOT.fontSize);
    ylabel(pYlab{pIdx},   'FontSize', PLOT.fontSize);
    xlim([min(licks)-0.5, max(licks)+0.5]);
    set(gca, 'FontSize', PLOT.fontSize);
    set(gcf, 'Position', [75 + 330*(pIdx-1), 75, 300, 700]);
    box off
    hold off;

% ---- report ----
    fprintf('\n%s -- %s vs contact 1, %s-tailed', pName{pIdx}, STATS.test, STATS.tail);
    if STATS.bh
        fprintf(', BH over %d contacts\n', sum(okp));
    else
        fprintf(', uncorrected\n');
    end
    disp(table(licks', mu', ci', pRaw', pAdj', n', ...
        'VariableNames', {'Contact','mean','ci95','p','pAdj','n_sess'}));
    sigC = licks(pAdj < STATS.alpha);
    if isempty(sigC)
        fprintf('  starred contacts: none\n');
    else
        fprintf('  starred contacts: %s\n', num2str(sigC));
    end
end

%% HELPERS

% findConsecutiveSets removed; shared helpers copied verbatim from the R16 script.

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

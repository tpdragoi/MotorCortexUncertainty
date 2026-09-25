%% F3B_lickDuration.m
%  Duration of successive tongue protrusions, Double Reward Task.
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

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


% ---- SELF-CONTAINED: data and functions come only from this folder ----
% Data\R16 holds the exported <ANM>_<DATE>_obj.mat / _kin.mat; shared\ holds the
% pipeline functions it needs (shared\pipelineCopies). uninstructedMovements_v2-main
% and the raw data tree are not on the path.
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'loadBehavSession.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));
dataDir = fullfile(v2Root, 'Data', 'R16');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s', dataDir);

%% PARAMETERS
params.alignEvent = 'goCue';
params.behav_only = 1;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 6'};
params.condition(end+1) = {'hit==1'};

params.tmin   = -2.5;
params.tmax   = 4;
params.dt     = 1/200;
params.smooth = 1;
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

%% CONFIG
% BOUT DETECTION. Every threshold is a TIME; sample counts are derived from
% params.dt, so changing the sampling rate cannot silently change what counts
% as a lick.
cfg.nanBridge_sec   = 0.020;   % interpolate tracking gaps up to this long
cfg.minBout_sec     = 0.030;   % a contact run shorter than this is not a lick
cfg.maxBout_sec     = 0.400;   % ... and one still longer than this AFTER the
% contact-driven split is an artefact, not a lick.
% not the protrusion-duration filter the old
% 150/180 ms limits were acting as (see note 1).
cfg.contactTol_sec  = 0.020;   % a run must fall within this of a lick event
cfg.maxLickIdx      = 35;   % lick-indexed columns kept per trial
cfg.nLicksSummary   = 15;   % licks carried into the summary matrices

% ANALYSIS / FIGURE
cfg.nLicksPlot  = 8;   % licks shown and tested
cfg.alpha       = 0.05;
cfg.plotChoices = {'mean','sd'};   % both measures are plotted, one figure each
cfg.refLick     = 2;   % [3B] reference protrusion: every other protrusion is compared with lick 2
cfg.tail        = 'right';   % [3B] one-tailed: lick i LONGER than lick 2 -- lick 1 > lick 2 (drop after
% the first reward) and licks 3-8 > lick 2 (rise again before the second reward)
cfg.durScale    = 1000;   % seconds -> ms for display
cfg.durUnit     = 'ms';

% TRIAL CAPS: {animal, date, last usable trial}
% TRIAL SELECTION: hit trials of BOTH reward types (reward on lick 1 and on
% lick 6), and only trials with at least cfg.minContacts port contacts after the
% go cue, so every kept trial reaches the second reward contact.
cfg.minContacts = 6;

cfg.exclusions = {
    'TD13d', '2024-11-11', 278;
    'TD8d',  '2024-09-07', 313;
    'TD8d',  '2024-09-09', 298;
    'TD4f',  '2023-12-16', 203;
    'TD25d', '2025-07-17', 199;
};

%% SESSIONS TO LOAD
%% Expected for this panel: n = 26 sessions, 6 animals

% one empty placeholder per session slot; the ones a task uses are
% filled in below and the rest drop out of the all_meta concatenation
[meta, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11, ...
    meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19, meta20, meta21, ...
    meta22, meta23, meta24, meta25, meta26, meta27, meta28, meta29, meta30, meta31, ...
    meta32, meta33, meta34, meta35, meta36, meta37, meta38, meta39, meta40, ...
    meta41] = deal([]);
% TD1d R1
date = '2023-06-03';
meta1 = struct('anm','YH2','date',date);   % was loadYH2_many
date = '2023-06-04';
meta2 = struct('anm','YH2','date',date);   % was loadYH2_many
date = '2023-06-05';
meta3 = struct('anm','YH2','date',date);   % was loadYH2_many
date = '2023-06-06';
meta4 = struct('anm','YH2','date',date);   % was loadYH2_many

date = '2023-05-10';
meta11 = struct('anm','YH1','date',date);   % was loadYH1_many
date = '2023-05-11';
meta12 = struct('anm','YH1','date',date);   % was loadYH1_many
date = '2023-05-12';
meta13 = struct('anm','YH1','date',date);   % was loadYH1_many
date = '2023-06-03';
meta15 = struct('anm','YH1','date',date);   % was loadYH1_many
date = '2023-06-04';
meta16 = struct('anm','YH1','date',date);   % was loadYH1_many
date = '2023-12-06';
meta17 = struct('anm','TD4f','date',date);   % was loadTD4f_many
date = '2023-12-07';
meta18 = struct('anm','TD4f','date',date);   % was loadTD4f_many
date = '2023-12-09';
meta19 = struct('anm','TD4f','date',date);   % was loadTD4f_many
date = '2023-12-13';
meta20 = struct('anm','TD4f','date',date);   % was loadTD4f_many
date = '2023-12-16';
meta22 = struct('anm','TD4f','date',date);   % was loadTD4f_many

date = '2024-01-15';
meta23 = struct('anm','TD7f','date',date);   % was loadTD7f_many
date = '2024-01-17';
meta25 = struct('anm','TD7f','date',date);   % was loadTD7f_many
date = '2024-01-18';
meta26 = struct('anm','TD7f','date',date);   % was loadTD7f_many
date = '2024-01-22';
meta27 = struct('anm','TD7f','date',date);   % was loadTD7f_many
date = '2024-01-24';
meta29 = struct('anm','TD7f','date',date);   % was loadTD7f_many

% TD18d R1
date = '2025-06-25';
meta30 = struct('anm','TD24d','date',date);   % was loadTD24_many
date = '2025-06-27';
meta32 = struct('anm','TD24d','date',date);   % was loadTD24_many
date = '2025-06-28';
meta33 = struct('anm','TD24d','date',date);   % was loadTD24_many

date = '2025-07-14';
meta36 = struct('anm','TD25d','date',date);   % was loadTD25_many
date = '2025-07-15';
meta37 = struct('anm','TD25d','date',date);   % was loadTD25_many
date = '2025-07-16';
meta38 = struct('anm','TD25d','date',date);   % was loadTD25_many
date = '2025-07-17';
meta39 = struct('anm','TD25d','date',date);   % was loadTD25_many

y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12;meta13 ...
    ;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24;meta25;meta26 ...
    ;meta27;meta28;meta29;meta30;meta31;meta32;meta33;meta34;meta35;meta36;meta37;meta38;meta39;meta40];

%% MAIN LOOP
nSess    = numel(all_meta);
allDur   = cell(1, nSess);   allDurSD = cell(1, nSess);
allLen   = cell(1, nSess);   allILI   = cell(1, nSess);
allAnm   = cell(1, nSess);   allDate  = cell(1, nSess);
allDur6  = cell(1, nSess);   allDurSD6 = cell(1, nSess);   % reward-on-lick-6 trials
allLen6  = cell(1, nSess);   allILI6   = cell(1, nSess);


for sessnum = 1:nSess

    clear obj kin me Length

    meta         = all_meta(sessnum, 1);
    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);
    sessix = 1;

% ---- trial cap (per-session exclusions, applied once, in one place) ----
    capN     = trialCapFor(obj.pth.anm, obj.pth.dt, cfg.exclusions);
    condtrix = (1:obj.bp.Ntrials)';
    if ~isnan(capN)
        condtrix(condtrix > capN) = [];
        logf('  [cap] %s %s: keeping trials 1-%d of %d\n', ...
            obj.pth.anm, obj.pth.dt, capN, obj.bp.Ntrials);
    end

    kinix = find(strcmp(kin(sessix).featLeg, 'tongue_length'));
    assert(~isempty(kinix), 'tongue_length not found in kin.featLeg');
    Length = kin(sessix).dat(:, condtrix, kinix);

% ---- trial sets ----
% NOTE: Length has been indexed by condtrix, so its columns are 1..numel(condtrix).
% Trial lists must be expressed in the SAME indexing, which is what the
% ismember/find step below does. The old code compared raw trial numbers
% against a capped matrix, which silently shifted trials in capped sessions.
    allTr = (1:obj.bp.Ntrials)';
    hitTr = allTr(obj.bp.hit == 1);
    P1 = intersect(hitTr, allTr(obj.bp.rewardedLick == 1));   % single-reward trials
    P4 = intersect(hitTr, allTr(obj.bp.rewardedLick == 6));   % double-reward trials
% keep only trials with >= cfg.minContacts contacts after the go cue
    nPost = zeros(obj.bp.Ntrials, 1);
    for t = 1:obj.bp.Ntrials
        lks = obj.bp.ev.lickL{t};
        if isempty(lks), continue; end
        nPost(t) = sum(lks(:) > obj.bp.ev.goCue(t));
    end
    okC = allTr(nPost >= cfg.minContacts);
    P1 = intersect(P1, okC);
    P4 = intersect(P4, okC);
    P  = {P1, P4};   % both trial types

    for k = 1:numel(P)
        P{k} = find(ismember(condtrix, P{k}));   % -> columns of Length
    end

% ---- per-condition extraction ----
    for k = 1:numel(P)
        trialList = P{k}(:)';
        dur_all = [];  len_all = [];  ili_all = [];

        for tr = trialList
            trAbs = condtrix(tr);   % original trial number, for obj.bp.ev
            [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial( ...
                Length(:,tr), obj.time, obj.bp.ev.lickL{trAbs,1}, obj.bp.ev.goCue(trAbs), ...
                params.dt, cfg);
            if why > 0
                if why < 5, continue; end   % why == 5 is a collision warning only
            end
            dur_all = [dur_all; durRow];   %#ok<AGROW>
            len_all = [len_all; lenRow];   %#ok<AGROW>
            ili_all = [ili_all; iliRow];   %#ok<AGROW>
        end

        nL = cfg.nLicksSummary;
        if isempty(dur_all)
            warning('session %d (%s %s) condition %d produced no usable trials.', ...
                sessnum, obj.pth.anm, obj.pth.dt, k);
            dur_all = nan(1,nL); len_all = nan(1,nL); ili_all = nan(1,nL);
        end
        dur_all = padTo(dur_all, nL);
        len_all = padTo(len_all, nL);
        ili_all = padTo(ili_all, nL);

% ONE VALUE PER SESSION PER LICK: mean (and SD) ACROSS TRIALS. The tests
% below are over sessions, so the session is the unit of analysis and
% individual trials must not enter the test as independent samples.
        valDur   = mean(dur_all(:,1:nL), 1, 'omitnan')';
        valDurSD = std( dur_all(:,1:nL), 0, 1, 'omitnan')';
        valLen   = mean(len_all(:,1:nL), 1, 'omitnan')';
        valILI   = mean(ili_all(:,1:nL), 1, 'omitnan')';

        if k == 1
            allDur{sessnum}   = valDur;
            allDurSD{sessnum} = valDurSD;
            allLen{sessnum}   = valLen;
            allILI{sessnum}   = valILI;
        else
            allDur6{sessnum}   = valDur;   % reward-on-lick-6 trials
            allDurSD6{sessnum} = valDurSD;
            allLen6{sessnum}   = valLen;
            allILI6{sessnum}   = valILI;
        end
    end

    allAnm{sessnum}  = obj.pth.anm;
    allDate{sessnum} = obj.pth.dt;
    fprintf('Session %d\n', sessnum);
end

anmList = unique(allAnm(~cellfun(@isempty, allAnm)));

%% FIGURE + STATISTICS
% Fig. 3B: each protrusion compared to the FIRST, Wilcoxon signed-rank,
% paired within session (one value per session per lick), Benjamini-Hochberg
% FDR across the m = 7 comparisons.
% TAIL. cfg.tail is set in the config block at the top:
%   'both'  -- two-tailed
%   'left'  -- one-tailed, later protrusion SHORTER than the first
%   'right' -- one-tailed, later protrusion LONGER than the first
% READ THIS FOR Fig. 3B. The Results text makes TWO directional claims about
% this panel: duration decreases after the first reward, then increases again
% before the second reward at C6. A single one-tailed test cannot support
% both. cfg.tail = 'left' tests the decrease, which is the direction
% established earlier in the paper (Fig. 1B) and so is the one a pre-stated
% hypothesis covers. If the pre-C6 increase is to be claimed with an
% asterisk, it needs its own pre-stated hypothesis and its own family; the
% two-tailed column printed below is the safe fallback for the whole panel.
% The two-tailed p is printed alongside as a reference column regardless, so
% the cost of the one-sided choice is visible. The tail must be fixed by the
% hypothesis BEFORE looking at the data -- picking it per lick from the sign of
% the difference is not a one-tailed test, it is a two-tailed test with the
% p-value halved.
% BOTH trial types are plotted and tested. Within each type, every
% protrusion is compared with that type's SECOND (cfg.refLick = 2; paired
% signed-rank across sessions, one-tailed lick i > lick 2, BH over m = 7). Reward on lick 1 = red,
% reward on lick 6 = blue, as in the R14 figure. A paired two-tailed
% lick-1-reward vs lick-6-reward comparison at each lick is printed as well.
for pc = 1:numel(cfg.plotChoices)
    plotChoice = cfg.plotChoices{pc};

    switch lower(plotChoice)
        case 'mean', M = {cell2mat(allDur)   * cfg.durScale, cell2mat(allDur6)   * cfg.durScale};
                     yLab = sprintf('Mean protrusion duration (%s)', cfg.durUnit);
        case 'sd',   M = {cell2mat(allDurSD) * cfg.durScale, cell2mat(allDurSD6) * cfg.durScale};
                     yLab = sprintf('SD of protrusion duration (%s)', cfg.durUnit);
        otherwise,   error('plotChoice must be ''mean'' or ''sd''.');
    end
    tName = {'reward on lick 1', 'reward on lick 6'};
    tCol  = {[0.90 0.45 0.70], [0.45 0.15 0.65]};   % R1 pink, R16 purple

    switch lower(cfg.tail)
        case 'both',  tailA = {};                 tailWord = 'two-tailed';
        case 'left',  tailA = {'tail','left'};    tailWord = sprintf('one-tailed (lick i < lick %d)', cfg.refLick);
        case 'right', tailA = {'tail','right'};   tailWord = sprintf('one-tailed (lick i > lick %d)', cfg.refLick);
        otherwise,    error('cfg.tail must be ''both'', ''left'' or ''right''.');
    end

    licksTested = setdiff(1:cfg.nLicksPlot, cfg.refLick);   % [3B] every protrusion except the reference (m = 7)
    mu = cell(1,2);  ci = cell(1,2);  sig = cell(1,2);  nSessUsed = 0;
    for k = 1:2
        mat = M{k}(1:cfg.nLicksPlot, :);   % licks x sessions
        M{k} = mat;
        nSessUsed = max(nSessUsed, sum(any(isfinite(mat),1)));
        nOK   = sum(isfinite(mat), 2);
        mu{k} = mean(mat, 2, 'omitnan');
        ci{k} = std(mat, 0, 2, 'omitnan') ./ sqrt(max(nOK,1)) .* tinv(1 - cfg.alpha/2, max(nOK - 1, 1));

        p = nan(1, numel(licksTested));  p2 = p;  nPair = zeros(1, numel(licksTested));
        for ii = 1:numel(licksTested)
            a1 = mat(licksTested(ii), :);  b1 = mat(cfg.refLick, :);   % [3B] vs lick 2
            ok = isfinite(a1) & isfinite(b1);
            nPair(ii) = sum(ok);
            if nPair(ii) < 2, continue; end
            p(ii)  = signrank(a1(ok), b1(ok), tailA{:});
            p2(ii) = signrank(a1(ok), b1(ok));
        end
        [p, pAdj, sig{k}] = bhCorrect(p, cfg.alpha);

        fprintf('\n%s\n', repmat('=',1,78));
        fprintf('Fig. 3B | %s | Wilcoxon signed-rank, each protrusion vs lick %d (%s)\n', tName{k}, cfg.refLick, plotChoice);
        fprintf('%s | m = %d comparisons | trials with >= %d contacts\n', tailWord, numel(licksTested), cfg.minContacts);
        fprintf('%s\n', repmat('-',1,78));
        fprintf('%-6s %-8s %-10s %-11s %-11s %-6s %s\n', 'Lick', 'n sess', 'mean', 'p (raw)', 'p (BH)', 'sig', 'p two-tailed');
        for ii = 1:numel(licksTested)
            fprintf('%-6d %-8d %-10.2f %-11.4f %-11.4f %-6s %.4f\n', licksTested(ii), nPair(ii), ...
                mu{k}(licksTested(ii)), p(ii), pAdj(ii), string(sig{k}(ii)), p2(ii));
        end
    end

% reference: lick-1-reward vs lick-6-reward trials at each lick, paired by session
    pB = nan(1, cfg.nLicksPlot);
    for lick = 1:cfg.nLicksPlot
        a1 = M{1}(lick,:);  b1 = M{2}(lick,:);  ok = isfinite(a1) & isfinite(b1);
        if sum(ok) >= 2, pB(lick) = signrank(a1(ok), b1(ok)); end
    end
    [pB, pBadj] = bhCorrect(pB, cfg.alpha);
    fprintf('\n%s\nreference: reward on lick 1 vs reward on lick 6, paired signed-rank, two-tailed, BH m = %d\n', ...
        repmat('-',1,78), cfg.nLicksPlot);
    fprintf('  %s\n', strjoin(arrayfun(@(L,q) sprintf('L%d q=%.4f', L, q), 1:cfg.nLicksPlot, pBadj, 'UniformOutput', false), ' | '));
    fprintf('%s\n', repmat('=',1,78));

    figure('Color','w','Position',[75 + 380*(pc-1), 75, 360, 620]); hold on
    h = gobjects(1,2);
    for k = 1:2
        for i = 1:cfg.nLicksPlot
            fill([i-0.22 i+0.22 i+0.22 i-0.22], ...
                [mu{k}(i)-ci{k}(i) mu{k}(i)-ci{k}(i) mu{k}(i)+ci{k}(i) mu{k}(i)+ci{k}(i)], ...
                tCol{k}, 'FaceAlpha',0.30, 'EdgeColor','none');
            h(k) = plot(i, mu{k}(i), 'o', 'Color',tCol{k}, 'MarkerFaceColor',tCol{k}, 'MarkerSize',8);
        end
    end
    yl = ylim;
    for k = 1:2
        for ii = 1:numel(licksTested)
            if sig{k}(ii)
                text(licksTested(ii) + (k-1.5)*0.3, yl(2) - 0.03*range(yl), '*', 'Color', tCol{k}, ...
                    'HorizontalAlignment','center', 'FontSize',16, 'FontWeight','bold');
            end
        end
    end
    xlabel('Lick number'); ylabel(yLab);
    xlim([0.5 cfg.nLicksPlot+0.5]); xticks(1:cfg.nLicksPlot);
    legend(h, tName, 'Location','best', 'Box','off');
    title(sprintf('Fig. 3B  (n = %d sessions, %d animals)', nSessUsed, numel(anmList)), ...
        'FontSize',10, 'FontWeight','normal');
    set(gca,'TickDir','out'); box off

end

%% LOCAL FUNCTIONS

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

function M = padTo(M, n)
    if size(M,2) < n, M(:, end+1:n) = NaN; end
end

function capN = trialCapFor(anm, dte, caps)
    capN = NaN;
    for r = 1:size(caps,1)
        if strcmpi(caps{r,1}, anm) && strcmpi(caps{r,2}, dte)
            capN = caps{r,3};  return
        end
    end
end

function [p, pAdj, sig] = bhCorrect(p, alpha)
% Benjamini-Hochberg adjusted p-values. m is numel(p) -- the number of licks in
% this family -- which is the number the Methods should quote.
    m = numel(p);
    [ps, ix] = sort(p);
    pAdj = nan(1, m);
    for i = 1:m
        pAdj(ix(i)) = min(ps(i:end) .* m ./ (i:m));
    end
    pAdj = min(pAdj, 1);
    sig  = pAdj < alpha;
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

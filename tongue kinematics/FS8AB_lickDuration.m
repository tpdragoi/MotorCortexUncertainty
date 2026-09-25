%% FS8AB_lickDuration.m
%  Duration of successive tongue protrusions, Simple Reward Task, first five days of training.
%  Protrusions are detected from tongue length, summarized per session, then
%  averaged across sessions.
%  READS
%    Data_1\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data_1\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/200
%    params.smooth      1
%    params.quality     {'good'}
%    params.lowFR       0.001
%    params.window      -2.5 to 4 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


% ---- SELF-CONTAINED: data and functions come only from this folder ----
% Data\Learning holds the exported <ANM>_<DATE>_obj.mat / _kin.mat; shared\ holds the
% pipeline functions it needs (shared\pipelineCopies). uninstructedMovements_v2-main
% and the raw data tree are not on the path.
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'loadBehavSession.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));
dataDir = fullfile(v2Root, 'Data_1', 'Learning');   % updated obj/kin files
if ~exist(dataDir, 'dir'), dataDir = fullfile(v2Root, 'Data', 'Learning'); end

%% PARAMETERS
params.alignEvent = 'goCue';
params.behav_only = 1;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.001;   % minimum mean firing rate, Hz

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1'};

params.tmin   = -2.5;
params.tmax   = 4;
params.dt     = 1/200;   % set to 1/100 to match the original Learning script (the other four use 1/200)
params.smooth = 1;   % was 10 (neural-only setting; no effect on kinematics)
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
cfg.nanBridge_sec   = 0.020;
cfg.minBout_sec     = 0.030;
cfg.maxBout_sec     = 0.400;
cfg.contactTol_sec  = 0.020;

% WHICH EXTRACTION. 'original' reproduces the lick detection of the
% original Learning script exactly: contacts re-referenced to the FIRST CONTACT
% and looked up on obj.time (which is go-cue aligned), findConsecutiveSets with
% fixed-size chopping, 30-180 ms runs, a run kept if it comes within
% cfg.origTolSamps samples of a lick, columns = order of the surviving runs.
% 'fixed' = the shared contactRuns extraction used by VTA/R1/R14/R16.
cfg.extraction   = 'original';   % 'original' | 'fixed'
cfg.origMin_ms   = 30;   % original: minn = 30/params.dt/1000
cfg.origMax_ms   = 180;   % original: maxx = 180/params.dt/1000
cfg.origTolSamps = 2;   % original: abs(tval - closest_indices) <= 2 (samples)
cfg.origStaleLicks = false;   % true = as the ORIGINAL: closest_indices is never cleared, so
% lick positions left over from earlier trials (and sessions)
% with more licks still count as matches. false = rebuilt per trial.
cfg.maxLickIdx      = 35;
cfg.nLicksSummary   = 15;

% TRIAL SELECTION: both outcomes, as requested.
cfg.hitOnly = false;   % false = hit == 1 | hit == 0

% ANALYSIS / FIGURE
cfg.nLicksPlot  = 8;
cfg.alpha       = 0.05;
cfg.test        = 'ranksum';   % 'ranksum' (unpaired) | 'signrank' (paired)

% One row per measure, each with its own tail. Both are tested ACROSS ANIMALS:
%   row 1  mean duration   two-tailed   (day 5 differs from day 1)
%   row 2  SD of duration  one-tailed   (day 5 less variable than day 1)
cfg.panels      = { 'mean', 'both' ; 'sd', 'right' };
cfg.figPos      = [0.06 0.18 0.88 0.62];
cfg.durScale    = 1000;
cfg.durUnit     = 'ms';
cfg.yLimMean    = [40 90];   % [] = auto. The original pinned these.
cfg.yLimSD      = [0  25];   % [] = auto
cfg.nSegGrad    = 100;   % segments in the black-to-red group-mean gradient

cfg.exclusions = {
    'TD4l', '2025-06-05', 123;
};
% NOTE: the old script wrote this cap for animal 'TDl4', but the loader is
% loadTD4l_many, so the string never matched and the cap never applied. It is
% spelled 'TD4l' here. Confirm against obj.pth.anm -- if that session prints
% no [cap] line when switched on, the spelling is still wrong.

%% SESSIONS TO LOAD
% loader | date | animal label | day (1 or 5) | use
% implicit in which meta slot a date landed in.
sessionTable = {
    'loadTD3l_many', '2025-02-03', 'TD3l', 1, true
    'loadTD3l_many', '2025-02-04', 'TD3l', 2, false
    'loadTD3l_many', '2025-02-05', 'TD3l', 3, false
    'loadTD3l_many', '2025-02-06', 'TD3l', 4, false
    'loadTD3l_many', '2025-02-07', 'TD3l', 5, true

    'loadTD2l_many', '2025-02-03', 'TD2l', 1, true
    'loadTD2l_many', '2025-02-04', 'TD2l', 2, false
    'loadTD2l_many', '2025-02-05', 'TD2l', 3, false
    'loadTD2l_many', '2025-02-06', 'TD2l', 4, false
    'loadTD2l_many', '2025-02-07', 'TD2l', 5, true

    'loadTD4l_many', '2025-06-03', 'TD4l', 1, true
    'loadTD4l_many', '2025-06-04', 'TD4l', 2, false
    'loadTD4l_many', '2025-06-05', 'TD4l', 3, false
    'loadTD4l_many', '2025-06-06', 'TD4l', 4, false
    'loadTD4l_many', '2025-06-07', 'TD4l', 5, true

    'loadTD5l_many', '2025-08-21', 'TD5l', 1, true
    'loadTD5l_many', '2025-08-22', 'TD5l', 2, false
    'loadTD5l_many', '2025-08-23', 'TD5l', 3, false
    'loadTD5l_many', '2025-08-24', 'TD5l', 4, false
    'loadTD5l_many', '2025-08-25', 'TD5l', 5, true
};

% ---- build all_meta from the table; the metadata travels with the row ----
all_meta = [];
sessAnm  = {};
sessDate = {};
sessTag  = [];
for r = 1:size(sessionTable,1)
    if ~sessionTable{r,5}, continue; end
    m        = struct('anm', sessionTable{r,3}, 'date', sessionTable{r,2});   % was feval(loader)
    all_meta = [all_meta; m];   %#ok<AGROW>
    sessAnm{end+1}  = sessionTable{r,3};   %#ok<SAGROW>
    sessDate{end+1} = sessionTable{r,2};   %#ok<SAGROW>
    sessTag(end+1)  = sessionTable{r,4};   %#ok<SAGROW>
end
nSess = size(all_meta,1);
assert(nSess == numel(sessAnm), 'session labels and all_meta are out of step');

%% MAIN LOOP
allDur   = cell(1, nSess);   allDurSD = cell(1, nSess);
allLen   = cell(1, nSess);   allILI   = cell(1, nSess);
allAnm   = cell(1, nSess);   allDate  = cell(1, nSess);
allDurCI = cell(1, nSess);   allSDCI  = cell(1, nSess);
allNtr   = cell(1, nSess);


origClosest = [];   % carried across trials AND sessions, as closest_indices was
for sessnum = 1:nSess

    clear obj kin me Length

    meta         = all_meta(sessnum, 1);
    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);
    sessix = 1;

% ---- trial cap ----
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

% ---- TRIAL SELECTION: both outcomes (note B) ----
    allTr = (1:obj.bp.Ntrials)';
    if cfg.hitOnly
        P1 = allTr(obj.bp.hit == 1);
    else
        P1 = allTr(obj.bp.hit == 1 | obj.bp.hit == 0);
    end
    P = {P1};

% Length is indexed by condtrix, so its columns are 1..numel(condtrix).
% Trial lists must be expressed in the SAME indexing.
    for k = 1:numel(P)
        P{k} = find(ismember(condtrix, P{k}));
    end

    for k = 1:numel(P)
        trialList = P{k}(:)';
        dur_all = [];  len_all = [];  ili_all = [];

        for tr = trialList
            trAbs = condtrix(tr);
            if strcmpi(cfg.extraction, 'original')
                [durRow, lenRow, iliRow, why, nDrop, origClosest] = lickFeaturesOriginal( ...
                    Length(:,tr), obj.time, obj.bp.ev.lickL{trAbs,1}, obj.bp.ev.goCue(trAbs), ...
                    params.dt, cfg, origClosest);
            else
                [durRow, lenRow, iliRow, why, nDrop] = lickFeaturesForTrial( ...
                    Length(:,tr), obj.time, obj.bp.ev.lickL{trAbs,1}, obj.bp.ev.goCue(trAbs), ...
                    params.dt, cfg);
            end
            if why > 0 && why < 5, continue; end
            dur_all = [dur_all; durRow];   %#ok<AGROW>
            len_all = [len_all; lenRow];   %#ok<AGROW>
            ili_all = [ili_all; iliRow];   %#ok<AGROW>
        end

        nL = cfg.nLicksSummary;
        if isempty(dur_all)
            warning('session %d (%s %s) produced no usable trials.', ...
                sessnum, obj.pth.anm, obj.pth.dt);
            dur_all = nan(1,nL); len_all = nan(1,nL); ili_all = nan(1,nL);
        end
        dur_all = padTo(dur_all, nL);
        len_all = padTo(len_all, nL);
        ili_all = padTo(ili_all, nL);

% ONE VALUE PER SESSION PER LICK. The tests below are over sessions, so
% the session is the unit of analysis and individual trials must not
% enter the test as independent samples.
        allDur{sessnum}   = mean(dur_all(:,1:nL), 1, 'omitnan')';
        allDurSD{sessnum} = std( dur_all(:,1:nL), 0, 1, 'omitnan')';
        allLen{sessnum}   = mean(len_all(:,1:nL), 1, 'omitnan')';
        allILI{sessnum}   = mean(ili_all(:,1:nL), 1, 'omitnan')';

% Within-session spread, for the per-animal error bars on the figure.
% Same definitions the original used: a t-based 95% CI for the mean and
% the SD/sqrt(2(n-1)) approximation for the CI of the SD.
        nTr = sum(~isnan(dur_all(:,1:nL)), 1);
        sdv = std(dur_all(:,1:nL), 0, 1, 'omitnan');
        allNtr{sessnum}   = nTr';
        allDurCI{sessnum} = (tinv(0.975, max(nTr-1,1)) .* sdv ./ sqrt(max(nTr,1)))';
        allSDCI{sessnum}  = (sdv ./ sqrt(2 .* max(nTr-1,1)))';
    end

    allAnm{sessnum}  = obj.pth.anm;
    allDate{sessnum} = obj.pth.dt;
    fprintf('Session %d\n', sessnum);
end

anmList = unique(allAnm(~cellfun(@isempty, allAnm)));

%% FIGURE + STATISTICS
% Plotted exactly as the original: one tile per lick, day 1 on the left and
% day 5 on the right, one grey line per animal carrying that session's
% within-session 95% CI, a black-to-red gradient connecting the group means,
% and a Benjamini-Hochberg corrected asterisk.
% The TEST is the UNPAIRED rank-sum (note F). The paired signed-rank is
% computed and printed alongside, with both floors, so the choice stays
% visible -- but the FIGURE is unchanged from the original.

% ---- pair the sessions by animal, from the table rather than by position ----
% The original hard-coded pairs = [1 2; 3 4; 5 6; 7 8], which is only correct
anmU  = unique(sessAnm, 'stable');
pairs = nan(numel(anmU), 2);
for a = 1:numel(anmU)
    i1 = find(strcmp(sessAnm, anmU{a}) & sessTag == 1, 1);
    i5 = find(strcmp(sessAnm, anmU{a}) & sessTag == 5, 1);
    if ~isempty(i1) && ~isempty(i5), pairs(a,:) = [i1 i5]; end
end
dropped = anmU(any(isnan(pairs),2));
if ~isempty(dropped)
    warning('animals without both days, excluded: %s', strjoin(dropped, ', '));
end
pairs(any(isnan(pairs),2), :) = [];
nAnimals = size(pairs,1);
assert(nAnimals >= 2, 'need at least two animals with both a day 1 and a day 5 session');

% ---- the two rows of the figure -------------------------------------------
% One row per measure, each with its own test. BOTH tests are ACROSS ANIMALS:
% each session contributes one value per lick, so every lick compares 4 day-1
% values with 4 day-5 values (n = 8 sessions, 4 animals). The asterisk is that
% test, Benjamini-Hochberg corrected across the plotted licks. Nothing else
% sets the asterisk.
c2v = @(c) horzcat(c{:})';   % -> sessions x licks

nLicksToPlot = min(cfg.nLicksPlot, size(c2v(allDur), 2));

% ---- floors, so the n = 4 ceiling is on the record ----
fprintf('\nEXACT FLOORS with n = %d animals per day\n', nAnimals);
fprintf('  paired signed-rank : one-tailed %.4f | two-tailed %.4f\n', ...
    2^(-nAnimals), 2^(1-nAnimals));
fprintf('  unpaired rank-sum  : one-tailed %.4f | two-tailed %.4f\n', ...
    1/nchoosek(2*nAnimals,nAnimals), 2/nchoosek(2*nAnimals,nAnimals));

% ---- colours, as in the original ----
col_grey = [0.65 0.65 0.65];
col_blk  = [0    0    0   ];
col_red  = [0.85 0.15 0.15];
nSeg     = cfg.nSegGrad;

nPanels = size(cfg.panels, 1);

figure('Color', 'w', ...
       'Name',        sprintf('Learning | mean and SD of lick duration | licks 1-%d', nLicksToPlot), ...
       'NumberTitle', 'off', ...
       'Units','normalized', 'Position', cfg.figPos);
tiledlayout(nPanels, nLicksToPlot, 'Padding','compact', 'TileSpacing','compact');

for pan = 1:nPanels
    plotChoice = cfg.panels{pan,1};
    tailChoice = cfg.panels{pan,2};

% ---- what this row plots ----
    switch lower(plotChoice)
        case 'mean'
            data_mu   = c2v(allDur)   * cfg.durScale;
            data_err  = c2v(allDurCI) * cfg.durScale;
            mainLabel = 'Mean duration (ms)';
            yLimRange = cfg.yLimMean;
        case 'sd'
            data_mu   = c2v(allDurSD) * cfg.durScale;
            data_err  = c2v(allSDCI)  * cfg.durScale;
            mainLabel = 'SD (ms)';
            yLimRange = cfg.yLimSD;
        otherwise
            error('cfg.panels column 1 must be ''mean'' or ''sd''.');
    end

    if isempty(yLimRange)
        lo  = min(data_mu(:) - data_err(:), [], 'omitnan');
        hi  = max(data_mu(:) + data_err(:), [], 'omitnan');
        pad = 0.10 * (hi - lo);
        yLimRange = [lo - pad, hi + pad];
    end

% ---- this row's tail ----
    switch lower(tailChoice)
        case 'both',  tailA = {};               tailWord = 'two-tailed';
        case 'left',  tailA = {'tail','left'};  tailWord = 'one-tailed (day1 < day5)';
        case 'right', tailA = {'tail','right'}; tailWord = 'one-tailed (day1 > day5)';
        otherwise,    error('cfg.panels column 2 must be ''both'', ''left'' or ''right''.');
    end

% ---- per-lick test across animals: 4 day-1 values vs 4 day-5 values ----
    pUnp  = nan(1, nLicksToPlot);
    pPair = nan(1, nLicksToPlot);
    for lick = 1:nLicksToPlot
        d1 = data_mu(pairs(:,1), lick);
        d5 = data_mu(pairs(:,2), lick);
        ok = isfinite(d1) & isfinite(d5);
        if sum(ok) < 2, continue; end
        pUnp(lick) = ranksum(d1(ok), d5(ok), tailA{:});
        if any(d1(ok) - d5(ok) ~= 0)
            pPair(lick) = signrank(d1(ok), d5(ok), tailA{:});
        end
    end
    [pUnp,  pUnpAdj ] = bhCorrect(pUnp,  cfg.alpha);
    [pPair, pPairAdj] = bhCorrect(pPair, cfg.alpha);

    if strcmpi(cfg.test, 'signrank')
        pShowAdj = pPairAdj;  testWord = 'paired signed-rank';
    else
        pShowAdj = pUnpAdj;   testWord = 'unpaired rank-sum';
    end
    h_bh = pShowAdj < cfg.alpha;

% ---- printed summary for this row ----
    fprintf('\n%s\n', repmat('=',1,84));
    fprintf('Learning | day 1 vs day 5 at each lick (%s)\n', plotChoice);
    fprintf('reported: %s, %s | across animals | m = %d | n = %d sessions per day\n', ...
        testWord, tailWord, nLicksToPlot, nAnimals);
    fprintf('%s\n', repmat('-',1,84));
    fprintf('%-6s %-10s %-10s %-11s %-11s %-11s %-11s %s\n', ...
        'Lick','day1','day5','p unpaired','pAdj unp','p paired','pAdj pair','sig');
    for lick = 1:nLicksToPlot
        fprintf('%-6d %-10.2f %-10.2f %-11.4f %-11.4f %-11.4f %-11.4f %s\n', lick, ...
            mean(data_mu(pairs(:,1),lick),'omitnan'), mean(data_mu(pairs(:,2),lick),'omitnan'), ...
            pUnp(lick), pUnpAdj(lick), pPair(lick), pPairAdj(lick), string(h_bh(lick)));
    end
    fprintf('%s\n', repmat('=',1,84));

% ---- draw this row ----
    for lick = 1:nLicksToPlot
        nexttile; hold on;

% per-animal grey lines + within-session 95% CI
        for a = 1:nAnimals
            val1 = data_mu (pairs(a,1), lick);
            val2 = data_mu (pairs(a,2), lick);
            ci1  = data_err(pairs(a,1), lick);
            ci2  = data_err(pairs(a,2), lick);

            errorbar(1, val1, ci1, '.', 'Color',col_grey, 'CapSize',6, 'LineWidth',1.2);
            errorbar(2, val2, ci2, '.', 'Color',col_grey, 'CapSize',6, 'LineWidth',1.2);
            plot([1 2],[val1 val2], '-', 'Color',col_grey, 'LineWidth',1.2);
        end

% group means, black -> red gradient
        grp_val1 = mean(data_mu(pairs(:,1), lick), 'omitnan');
        grp_val2 = mean(data_mu(pairs(:,2), lick), 'omitnan');

        xgrad = linspace(1, 2, nSeg+1);
        ygrad = linspace(grp_val1, grp_val2, nSeg+1);
        for sg = 1:nSeg
            frac   = (sg-1) / (nSeg-1);
            segCol = (1-frac)*col_blk + frac*col_red;
            plot(xgrad(sg:sg+1), ygrad(sg:sg+1), '-', 'Color', segCol, 'LineWidth', 3);
        end

        plot(1, grp_val1, 'o', 'Color',col_blk, 'MarkerFaceColor',col_blk, ...
             'MarkerSize',8, 'LineWidth',1.5);
        plot(2, grp_val2, 'o', 'Color',col_red, 'MarkerFaceColor',col_red, ...
             'MarkerSize',8, 'LineWidth',1.5);

        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels({'',''});

        text(1, yLimRange(1) - 0.12*(yLimRange(2)-yLimRange(1)), '1', ...
             'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','bold', ...
             'Color',col_blk, 'Clipping','off');
        text(2, yLimRange(1) - 0.12*(yLimRange(2)-yLimRange(1)), '5', ...
             'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','bold', ...
             'Color',col_red, 'Clipping','off');

        if pan == 1
            title(sprintf('Lick %d', lick));
        end
        if lick == 1
            ylabel(sprintf('%s  (95%% CI)', mainLabel));
        end

        ylim(yLimRange);

        if h_bh(lick)
            ystar = yLimRange(2) - 0.05*(yLimRange(2)-yLimRange(1));
            text(1.5, ystar, '*', 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top', 'FontSize',20, 'FontWeight','bold', 'Color','k');
        end

        box off; hold off;
    end
end

%% LOCAL FUNCTIONS

function [durRow, lenRow, iliRow, why, nDrop, closestOut] = lickFeaturesOriginal(lenTrace, tAxis, licksAbs, gcTime, dt, cfg, closestPrev)
% The per-trial lick extraction of the ORIGINAL Learning script, line
%   keep a set if any of its samples is within cfg.origTolSamps of a closest index
% closest_indices: with cfg.origStaleLicks = true it behaves as in the original --
% it is never cleared, entries 1..nLicks are overwritten and any entries beyond
% that from an earlier trial (or session) remain and still count as matches.
% With false it is rebuilt for every trial. Not copied: a trial with more than 35
% kept runs is trimmed to 35 instead of being dropped by the empty try.
    durRow = nan(1, cfg.maxLickIdx);
    lenRow = nan(1, cfg.maxLickIdx);
    iliRow = nan(1, cfg.maxLickIdx);
    why = 0;  nDrop = [0 0];
    closestOut = closestPrev;   % unchanged unless this trial reaches the lookup, as in the original

    postCue = licksAbs(licksAbs > gcTime);
    if isempty(postCue),       why = 3; return; end
    licks = postCue - postCue(1);
    if ~(numel(licks) > 2),    why = 1; return; end
    val = lenTrace(:);
    if all(isnan(val)),        why = 2; return; end

    minn = cfg.origMin_ms / dt / 1000;
    maxx = cfg.origMax_ms / dt / 1000;
    realIndices = find(~isnan(val));
    sets = findConsecutiveSets(realIndices, minn, maxx);

    if cfg.origStaleLicks
        closest_indices = closestPrev(:)';   % stale entries kept, as the original
    else
        closest_indices = zeros(1, numel(licks));   % rebuilt per trial
    end
    for j = 1:numel(licks)
        [~, closest_indices(j)] = min(abs(tAxis - licks(j)));
    end
    closestOut = closest_indices;
    keep = false(1, numel(sets));
    for jj = 1:numel(sets)
        keep(jj) = any(any(abs(sets{jj}(:)' - closest_indices(:)) <= cfg.origTolSamps));
    end
    sets = sets(keep);
    if isempty(sets),          why = 4; return; end

    nK = min(numel(sets), cfg.maxLickIdx);
    first_vals = nan(1, nK);
    for j = 1:nK
        durRow(j)     = numel(sets{j}) * dt;
        lenRow(j)     = max(val(sets{j}));
        first_vals(j) = sets{j}(1);
    end
    if nK > 1
        iliRow(1:nK-1) = diff(first_vals) * dt;
    end
end

function sets = findConsecutiveSets(indices, minLength, maxLength)
% a run longer than maxLength is CHOPPED into pieces of maxLength and every
% piece of at least minLength is kept.
    sets = {};
    currentSet = [];
    i = 1;
    while i <= length(indices)-1
        if indices(i+1) - indices(i) == 1
            currentSet = [currentSet, indices(i)];   %#ok<AGROW>
        else
            currentSet = [currentSet, indices(i)];   %#ok<AGROW>
            while length(currentSet) >= minLength
                truncatedSet = currentSet(1:min(length(currentSet), maxLength));
                sets{end+1} = truncatedSet;   %#ok<AGROW>
                currentSet = currentSet(min(length(currentSet), maxLength)+1:end);
            end
            currentSet = [];
        end
        i = i + 1;
    end
    currentSet = [currentSet, indices(end)];
    while length(currentSet) >= minLength
        truncatedSet = currentSet(1:min(length(currentSet), maxLength));
        sets{end+1} = truncatedSet;   %#ok<AGROW>
        currentSet = currentSet(min(length(currentSet), maxLength)+1:end);
    end
end

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

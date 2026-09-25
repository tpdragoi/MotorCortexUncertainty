%% F3JKL_opto.m
%  Tongue kinematics on control and photoinhibition trials.
%  Bilateral photoinhibition of tjM1 and ALM triggered at the fourth port contact.
%  Behaviour only; no spike data is loaded.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat and _laser.mat
%    through shared\loadBehavSession; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'firstLick'
%    params.dt          1/200
%    params.smooth      50
%    params.quality     {'ok','good','mua','great'}
%    params.lowFR       0.01
%    params.window      -1.5 to 3 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear, clc

%% self-contained root (_clean): reads only <v2>\Data and <v2>\shared
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root,'shared'),'dir')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root,'shared'));
addpath(fullfile(v2Root,'shared','pipelineCopies'));
dataDir = fullfile(v2Root, 'Data', 'C4Stim_R16');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s', dataDir);
% (needs <stem>_laser.mat next to each session: run exportLaserTrig.m once)

%  C4 STIM -- Reward 1/6 sessions
%  Stimulation is delivered at contact 4 (stimTrialGroups(obj, side, 4, task), local function at the end of this file).
%    FIG 1  lick angle (4 groups) + lick duration (Control vs Stim)
%    FIG 2  P(lick) Control vs Stim, zoomed contact range
%    FIG 3  P(lick) Control vs Stim, full contact range, with significance
%  Each figure shows BOTH error-bar conventions: bars over pooled TRIALS (the
%  original) and bars over SESSIONS (the unit the tests use). Same stats and
%  same stars in both; only the bars and the plotted central value change.
%  The per-session imagesc panels (one figure per session) are gone.
%  STATISTICS -- every t-test and the proportion z-test are replaced with
%  Wilcoxon tests computed at the SESSION level:
%        independent, so the session is the unit of analysis)
%  Benjamini-Hochberg correction runs across the contacts each figure shows.
%  BUGS FIXED FROM THE ORIGINAL (all of these changed results):
%   1. y_*a(end:30) = NaN overwrote the LAST detected lick with NaN on every
%      trial, because the range started at end rather than end+1. Every trial
%      silently lost its final contact. Set EXTRACT.reproduceLastLickNaN =
%      pattern happens to be identical so the P(lick) panels were unaffected,
%      but y4_length held the wrong numbers.
%   3. The block labelled "Two-tailed t-tests for ANGLE" read all_data{2,*},
%      which is DURATION, not angle. Its table header also said y1_vs_y3 and
%      y2_vs_y4 while the code compared y1 vs y2 and y3 vs y4. The comparison
%      the data row and the labels were wrong.
%      trial with no visible tongue reused the previous trial's lick indices.
%      row kept its natural length and the vertical concatenation would have
%      thrown a dimension error. Padding is now unconditional.
%   6. S21c = S21c(randperm(length(S21c), ctrlTrials)) errors whenever a
%      session has fewer than ctrlTrials control trials, and its result fed
%      only dead variables. Removed.
%  animal-level analysis would be n = 4 and would be.

%% OPTIONS

PLOT.fig1 = true;   % angle + duration vs contact number
PLOT.fig2 = false;   % P(lick), zoomed 2:8 -- OFF; set true to get it back
%            (fig3 below is the one kept: full range, with stars)
PLOT.fig3 = true;   % P(lick), full range, with stars

PLOT.angleLicks     = 1:12;   % contacts shown in fig 1 (and its BH family)
PLOT.probLicksZoom  = 2:8;   % contacts shown in fig 2 (and its BH family)
PLOT.probLicks      = 2:9;   % contacts shown in fig 3 (and its BH family)
% was 2:12; the BH family is now these 8 contacts

STIM.contact = 4;   % the contact stimulation is delivered on

% Control / Stim colours match the go-cue figures (black / light blue).
% The original used red / blue for the P(lick) panels; swap back here.
PLOT.colCtrl  = [0   0   0  ];
PLOT.colStim  = [0.5 0.7 1  ];
PLOT.colLctrl = [0   0   1  ];   % L_ctrl  dark blue
PLOT.colLstim = [0.5 0.7 1  ];   % L_stim  light blue
PLOT.colRctrl = [1   0   0  ];   % R_ctrl  dark red
PLOT.colRstim = [1   0.6 0.6];   % R_stim  light red

PLOT.fontSize = 12;

% ---- statistics -------------------------------------------------------
STATS.alpha   = 0.05;
STATS.bh      = true;   % BH across the contacts each figure shows
STATS.figTest = 'signrank';   % which p drives the stars: 'signrank' | 'ranksum'
STATS.tail    = 'both';   % applies to every comparison below

% ---- extraction -------------------------------------------------------
EXTRACT.startBin   = 290;   % first sample used, index into the trial trace
EXTRACT.minSeg     = 5;   % findConsecutiveSets minimum run length
EXTRACT.maxSeg     = 20;   % findConsecutiveSets maximum run length
EXTRACT.gapThresh  = 125;   % samples; a longer gap ends the bout
EXTRACT.maxLicks   = 30;   % columns kept per trial
EXTRACT.reproduceLastLickNaN = false;   % true = restore bug 1 above

%% PATHS

% d = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main';

% addpath(genpath(fullfile(d,'utils')))
% addpath(genpath(fullfile(d,'DataLoadingScripts')))
% addpath(genpath(fullfile(d,'funcs')))
% rmpath(genpath(fullfile(d,'fig1')));

% addpath 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\base code\functions_td'
% addpath 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\ObjVis\warp'
% addpath 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\base code\other_codes\functions_td'

%% PARAMETERS

params.alignEvent = 'firstLick';

params.timeWarp = 0;
params.nLicks   = 20;
params.lowFR    = 0.01;

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 1 & rewardedLick == 6'};
params.condition(end+1) = {'hit==1 & trialTypes == 2 & rewardedLick == 6'};
params.condition(end+1) = {'hit==1 & trialTypes == 3 & rewardedLick == 6'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 6'};
params.condition(end+1) = {'hit==1'};

params.tmin   = -1.5;
params.tmax   = 3;
params.dt     = 1/200;
params.smooth = 50;

params.quality    = {'ok','good','mua','great'};
params.behav_only = 1;

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

%% SESSIONS TO LOAD

% datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';

meta1 = []; meta2 = []; meta3 = []; meta4 = []; meta5 = [];
meta6 = []; meta7 = []; meta8 = []; meta9 = []; meta10 = [];
meta11 = []; meta12 = []; meta13 = []; meta14 = []; meta15 = [];
meta16 = []; meta17 = []; meta18 = []; meta19 = []; meta20 = [];
meta21 = []; meta22 = [];

date = '2023-11-18';
meta1 = struct('anm','TD4f','date',date);
date = '2023-11-19';
meta2 = struct('anm','TD4f','date',date);
date = '2023-11-20';
meta3 = struct('anm','TD4f','date',date);
date = '2023-11-25';
meta5 = struct('anm','TD4f','date',date);
date = '2023-11-26';
meta6 = struct('anm','TD4f','date',date);

date = '2023-11-18';
meta7 = struct('anm','TD5f','date',date);
date = '2023-11-19';
meta8 = struct('anm','TD5f','date',date);
date = '2023-11-20';
meta9 = struct('anm','TD5f','date',date);
date = '2023-11-21';
meta10 = struct('anm','TD5f','date',date);
date = '2023-11-25';
meta11 = struct('anm','TD5f','date',date);

date = '2023-12-08';
meta12 = struct('anm','TD6f','date',date);
date = '2023-12-09';
meta13 = struct('anm','TD6f','date',date);

date = '2023-12-08';
meta14 = struct('anm','TD7f','date',date);
date = '2023-12-09';
meta15 = struct('anm','TD7f','date',date);
date = '2023-12-14';
meta16 = struct('anm','TD7f','date',date);
date = '2023-12-15';
meta17 = struct('anm','TD7f','date',date);
date = '2023-12-16';
meta18 = struct('anm','TD7f','date',date);
date = '2023-12-17';
meta19 = struct('anm','TD7f','date',date);

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11; ...
            meta12;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22];

%% EXTRACT KINEMATICS
% Group order everywhere below:  1 = L_ctrl, 2 = L_stim, 3 = R_ctrl, 4 = R_stim

maxLicks = EXTRACT.maxLicks;
nSessTot = size(all_meta,1);
task     = 16;

SESS = struct('angle',{},'length',{},'duration',{},'lastLick',{});

y1_angle = []; y2_angle = []; y3_angle = []; y4_angle = [];
y1_length = []; y2_length = []; y3_length = []; y4_length = [];
y1_duration = []; y2_duration = []; y3_duration = []; y4_duration = [];

for sessnum = 1:nSessTot

    clear L_ctrl L_stim R_ctrl R_stim Length angle obj idxHit kin me

    meta         = all_meta(sessnum,1);

    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);   % [_clean]
    assert(isfield(obj,'sglx') && isfield(obj.sglx,'laserTrigIX'), ...
        ['%s: no laser trigger file (%s). The exported obj/kin files do not hold ' ...
         'obj.sglx.laserTrigIX -- run exportLaserTrig.m (v2 root) once to write it.'], ...
        params.slimStem, fullfile(dataDir, [params.slimStem '_laser.mat']));

    trialSet = (1:obj.bp.Ntrials)';

    sessix   = 1;
    condtrix = trialSet;

    kinix  = find(strcmp(kin(sessix).featLeg,'tongue_length'));
    Length = kin.dat(:,condtrix,kinix);

    kinix = find(strcmp(kin(sessix).featLeg,'tongue_angle'));
    angle = kin.dat(:,condtrix,kinix);

% ---- trial groups ----
    f = 0;

    [S21, S21c] = stimTrialGroups(obj, 1, STIM.contact, task);   % was find_StimTrials(obj, 1, f, STIM.contact, 6, task)
    idxHit = find(obj.bp.hit == 1);
    S21c   = S21c(ismember(S21c, idxHit));
    L_stim = S21;
    L_ctrl = S21c;

    [S21, S21c] = stimTrialGroups(obj, 3, STIM.contact, task);   % was find_StimTrials(obj, 3, f, STIM.contact, 6, task)
    idxHit = find(obj.bp.hit == 1);
    S21c   = S21c(ismember(S21c, idxHit));
    R_stim = S21;
    R_ctrl = S21c;

    groups_trials = {L_ctrl, L_stim, R_ctrl, R_stim};

    for g = 1:4

        trials = groups_trials{g};
        rows_a = [];
        rows_l = [];
        rows_d = [];
        ll     = nan(numel(trials),1);

        for i = 1:numel(trials)

            val_a = angle(EXTRACT.startBin:end, trials(i));
            val_l = Length(EXTRACT.startBin:end, trials(i));

            if all(isnan(val_a)), continue; end

            realIdx = find(~isnan(val_a));
            if isempty(realIdx), continue; end

            sets = findConsecutiveSets(realIdx, EXTRACT.minSeg, EXTRACT.maxSeg);
            nL   = numel(sets);
            if nL == 0, continue; end

            startBins = cellfun(@(s) s(1), sets);

            row_a = nan(1, maxLicks);
            row_l = nan(1, maxLicks);
            row_d = nan(1, maxLicks);

            nKeep = min(nL, maxLicks);
            for k = 1:nKeep
                seg = sets{k};
                sa  = val_a(seg);
                sl  = val_l(seg);
                [~, pk] = max(abs(sl));
                row_a(k) = sa(pk);
                row_l(k) = max(sl);
                row_d(k) = numel(seg);
            end

            if EXTRACT.reproduceLastLickNaN
                row_a(nKeep) = NaN;
                row_l(nKeep) = NaN;
                row_d(nKeep) = NaN;
            end

            rows_a = [rows_a; row_a];
            rows_l = [rows_l; row_l];
            rows_d = [rows_d; row_d];

            li = find(diff(startBins) > EXTRACT.gapThresh, 1);
            if isempty(li), li = nL; else, li = li + 1; end
            ll(i) = li;
        end

        SESS(sessnum).angle{g}    = rows_a;
        SESS(sessnum).length{g}   = rows_l;
        SESS(sessnum).duration{g} = rows_d;
        SESS(sessnum).lastLick{g} = ll;
    end

    y1_angle = [y1_angle; SESS(sessnum).angle{1}];
    y2_angle = [y2_angle; SESS(sessnum).angle{2}];
    y3_angle = [y3_angle; SESS(sessnum).angle{3}];
    y4_angle = [y4_angle; SESS(sessnum).angle{4}];

    y1_length = [y1_length; SESS(sessnum).length{1}];
    y2_length = [y2_length; SESS(sessnum).length{2}];
    y3_length = [y3_length; SESS(sessnum).length{3}];
    y4_length = [y4_length; SESS(sessnum).length{4}];

    y1_duration = [y1_duration; SESS(sessnum).duration{1}];
    y2_duration = [y2_duration; SESS(sessnum).duration{2}];
    y3_duration = [y3_duration; SESS(sessnum).duration{3}];
    y4_duration = [y4_duration; SESS(sessnum).duration{4}];

end

nS = numel(SESS);

fprintf('\nEXACT FLOORS with n = %d sessions\n', nS);
fprintf('  paired signed-rank : one-tailed %.3g | two-tailed %.3g\n', 2^(-nS), 2^(1-nS));
fprintf('  unpaired rank-sum  : one-tailed %.3g | two-tailed %.3g\n', ...
    1/nchoosek(2*nS,nS), 2/nchoosek(2*nS,nS));
if 2^(1-nS) > STATS.alpha
    fprintf('  *** the paired two-tailed test cannot reach p < %.3f at this n ***\n', STATS.alpha);
end

%% FIG 1 -- LICK ANGLE AND LICK DURATION
% Error bars are over TRIALS. The Wilcoxon tests below use SESSIONS as the
% unit, so the bars describe trial-to-trial spread while the stars describe
% consistency across sessions.

if PLOT.fig1

    licks   = PLOT.angleLicks;
    nLicksP = numel(licks);

% ---- per-session means ----
    A_sess = nan(nS, nLicksP, 4);
    D_sess = nan(nS, nLicksP, 2);
    for s = 1:nS
        for g = 1:4
            a = SESS(s).angle{g};
            if ~isempty(a), A_sess(s,:,g) = nanmean(a(:,licks), 1); end
        end
        dC = [SESS(s).duration{1}; SESS(s).duration{3}];
        dS = [SESS(s).duration{2}; SESS(s).duration{4}];
        dC(dC < 0) = NaN;
        dS(dS < 0) = NaN;
        if ~isempty(dC), D_sess(s,:,1) = nanmean(dC(:,licks), 1) * params.dt; end
        if ~isempty(dS), D_sess(s,:,2) = nanmean(dS(:,licks), 1) * params.dt; end
    end

% ---- pooled trials ----
    angle_pooled = {y1_angle, y2_angle, y3_angle, y4_angle};
    durC = [y1_duration; y3_duration];  durC(durC < 0) = NaN;
    durS = [y2_duration; y4_duration];  durS(durS < 0) = NaN;
    dur_pooled = {durC * params.dt, durS * params.dt};

    angCol  = {PLOT.colLctrl, PLOT.colLstim, PLOT.colRctrl, PLOT.colRstim};
    angName = {'L_{ctrl}','L_{stim}','R_{ctrl}','R_{stim}'};
    durCol  = {PLOT.colCtrl, PLOT.colStim};
    durName = {'Control','Stim'};

% ---- statistics, at the session level ----
    [p_angL, q_angL, p_angR, q_angR, p_dur, q_dur] = deal(nan(1,nLicksP));
    [n_angL, n_angR, n_dur] = deal(zeros(1,nLicksP));

    for i = 1:nLicksP
        xL = A_sess(:,i,1); yL = A_sess(:,i,2);
        ok = ~isnan(xL) & ~isnan(yL);  n_angL(i) = sum(ok);
        if sum(ok) >= 2
            if any(xL(ok) - yL(ok) ~= 0)
                p_angL(i) = signrank(xL(ok), yL(ok), 'tail', STATS.tail);
            end
            q_angL(i) = ranksum(xL(ok), yL(ok), 'tail', STATS.tail);
        end

        xR = A_sess(:,i,3); yR = A_sess(:,i,4);
        ok = ~isnan(xR) & ~isnan(yR);  n_angR(i) = sum(ok);
        if sum(ok) >= 2
            if any(xR(ok) - yR(ok) ~= 0)
                p_angR(i) = signrank(xR(ok), yR(ok), 'tail', STATS.tail);
            end
            q_angR(i) = ranksum(xR(ok), yR(ok), 'tail', STATS.tail);
        end

        x = D_sess(:,i,1); y = D_sess(:,i,2);
        ok = ~isnan(x) & ~isnan(y);  n_dur(i) = sum(ok);
        if sum(ok) >= 2
            if any(x(ok) - y(ok) ~= 0)
                p_dur(i) = signrank(x(ok), y(ok), 'tail', STATS.tail);
            end
            q_dur(i) = ranksum(x(ok), y(ok), 'tail', STATS.tail);
        end
    end

% ---- Benjamini-Hochberg across the displayed contacts ----
    pAll = {p_angL, q_angL, p_angR, q_angR, p_dur, q_dur};
    aAll = cell(size(pAll));
    for z = 1:numel(pAll)
        pv = pAll{z}(:); av = nan(size(pv)); okz = ~isnan(pv);
        if any(okz) && STATS.bh
            [ps, ordz] = sort(pv(okz));
            mm = numel(ps);
            ad = ps .* (mm ./ (1:mm)');
            ad = flipud(cummin(flipud(ad)));
            tmp = nan(mm,1); tmp(ordz) = min(ad,1);
            av(okz) = tmp;
        elseif ~STATS.bh
            av = pv;
        end
        aAll{z} = av';
    end
    pAdjL = aAll{1}; qAdjL = aAll{2};
    pAdjR = aAll{3}; qAdjR = aAll{4};
    pAdjD = aAll{5}; qAdjD = aAll{6};

    if strcmpi(STATS.figTest,'ranksum')
        sL = qAdjL; sR = qAdjR; sD = qAdjD;
    else
        sL = pAdjL; sR = pAdjR; sD = pAdjD;
    end

% ---- draw ----
    figure('Color','w','Units','normalized','Position',[0.06 0.08 0.80 0.82]);
    euName = {'trial','session'};

    for eu = 1   % error bars over TRIALS only

% ------------------------- ANGLE --------------------------------
        axA = subplot(1,2,1); hold on;

        for g = 1:4
            if eu == 2
                M  = A_sess(:,:,g);
                mu = nanmean(M,1);  N = sum(~isnan(M),1);
                se = nanstd(M,0,1) ./ sqrt(max(N,1));
            else
                y = angle_pooled{g};
                if isempty(y), continue; end
                y  = y(:, licks);
                mu = nanmean(y,1);  N = sum(~isnan(y),1);
                se = nanstd(y,0,1) ./ sqrt(max(N,1));
            end
            tc = 1.96*ones(size(N));
            tc(N > 1) = tinv(1 - STATS.alpha/2, N(N > 1) - 1);
            ci = se .* tc;
            v  = ~isnan(mu);
            xv = licks(v); mv = mu(v); cv = ci(v);
            plot(xv, mv, '-', 'Color', angCol{g}, 'LineWidth', 2);
            for i = 1:numel(xv)
                plot([xv(i) xv(i)], [mv(i)-cv(i) mv(i)+cv(i)], '-', ...
                     'Color', angCol{g}, 'LineWidth', 5);
                plot(xv(i), mv(i), 'o', 'MarkerSize', 6, ...
                     'MarkerFaceColor', angCol{g}, 'MarkerEdgeColor','none');
            end
        end

        yl = ylim(axA); ylim(axA, [yl(1), yl(2) + 0.12*diff(yl)]); yl = ylim(axA);
        for i = 1:nLicksP
            if sL(i) < STATS.alpha
                text(licks(i)-0.15, yl(2)-0.04*diff(yl), '*', 'Color', PLOT.colLctrl, ...
                     'FontSize', PLOT.fontSize+3, 'HorizontalAlignment','center');
            end
            if sR(i) < STATS.alpha
                text(licks(i)+0.15, yl(2)-0.04*diff(yl), '*', 'Color', PLOT.colRctrl, ...
                     'FontSize', PLOT.fontSize+3, 'HorizontalAlignment','center');
            end
        end
        if eu == 1
            for g = 1:4
                text(0.70, 0.06 + 0.07*(4-g), angName{g}, 'Units','normalized', ...
                     'Color', angCol{g}, 'FontSize', PLOT.fontSize-1, 'FontWeight','bold');
            end
        end

        line([STIM.contact STIM.contact], ylim(axA), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        xlim([0.5 max(licks)+0.5]); xticks(licks); xtickangle(0);
        xlabel('Contact number'); ylabel('Lick angle (deg)');
        title(sprintf('Angle  --  error bars over %ss', euName{eu}), 'FontWeight','normal');
        set(axA, 'TickDir','out', 'FontSize', PLOT.fontSize); box off

% ------------------------ DURATION ------------------------------
        axD = subplot(1,2,2); hold on;
        x_off = 0.15;

        for g = 1:2
            if eu == 2
                M  = D_sess(:,:,g);
                mu = nanmean(M,1);  N = sum(~isnan(M),1);
                se = nanstd(M,0,1) ./ sqrt(max(N,1));
            else
                y = dur_pooled{g};
                if isempty(y), continue; end
                y  = y(:, licks);
                mu = nanmean(y,1);  N = sum(~isnan(y),1);
                se = nanstd(y,0,1) ./ sqrt(max(N,1));
            end
            tc = 1.96*ones(size(N));
            tc(N > 1) = tinv(1 - STATS.alpha/2, N(N > 1) - 1);
            ci = se .* tc;
            xv = licks + ((-1)^g)*x_off;
            v  = ~isnan(mu);
            plot(xv(v), mu(v), '-', 'Color', durCol{g}, 'LineWidth', 2);
            for i = find(v)
                plot([xv(i) xv(i)], [max(mu(i)-ci(i),0) mu(i)+ci(i)], '-', ...
                     'Color', durCol{g}, 'LineWidth', 5);
                plot(xv(i), mu(i), 'o', 'MarkerSize', 6, ...
                     'MarkerFaceColor', durCol{g}, 'MarkerEdgeColor','none');
            end
        end

        yl = ylim(axD); ylim(axD, [yl(1), yl(2) + 0.12*diff(yl)]); yl = ylim(axD);
        for i = 1:nLicksP
            if sD(i) < STATS.alpha
                text(licks(i), yl(2)-0.04*diff(yl), '*', 'Color', 'k', ...
                     'FontSize', PLOT.fontSize+3, 'HorizontalAlignment','center');
            end
        end
        if eu == 1
            for g = 1:2
                text(0.72, 0.13 - 0.07*(g-1), durName{g}, 'Units','normalized', ...
                     'Color', durCol{g}, 'FontSize', PLOT.fontSize-1, 'FontWeight','bold');
            end
        end

        line([STIM.contact STIM.contact], ylim(axD), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        xlim([0.5 max(licks)+0.5]); xticks(licks); xtickangle(0);
        xlabel('Contact number'); ylabel('Lick duration (s)');
        title(sprintf('Duration  --  error bars over %ss', euName{eu}), 'FontWeight','normal');
        set(axD, 'TickDir','out', 'FontSize', PLOT.fontSize); box off
    end

    fprintf('\nFIG 1 -- p = paired signed-rank across sessions | q = unpaired rank-sum\n');
    fprintf('  tail = %s | stars use %s', STATS.tail, STATS.figTest);
    if STATS.bh
        fprintf(', BH over m = %d contacts\n', nLicksP);
    else
        fprintf(', UNCORRECTED\n');
    end

    disp(table(licks', p_angL', pAdjL', q_angL', qAdjL', n_angL', ...
        'VariableNames', {'Contact','p_L','pAdj_L','q_L','qAdj_L','n_L'}));
    disp(table(licks', p_angR', pAdjR', q_angR', qAdjR', n_angR', ...
        'VariableNames', {'Contact','p_R','pAdj_R','q_R','qAdj_R','n_R'}));
    disp(table(licks', p_dur', pAdjD', q_dur', qAdjD', n_dur', ...
        'VariableNames', {'Contact','p_dur','pAdj_dur','q_dur','qAdj_dur','n_dur'}));
    fprintf('  starred contacts -- L angle: %s | R angle: %s | duration: %s\n', ...
        num2str(licks(sL < STATS.alpha)), num2str(licks(sR < STATS.alpha)), ...
        num2str(licks(sD < STATS.alpha)));

end

%% FIG 2 AND FIG 3 -- P(LICK), CONTROL vs STIM
% compared the two proportions with a z-test; it is now a per-session
% proportion compared with a Wilcoxon test across sessions.
% Each figure is drawn twice side by side: left panel error bars over TRIALS
% (binomial on the pooled trials, the original), right panel over SESSIONS.

for figIdx = 1:2

    if figIdx == 1
        if ~PLOT.fig2, continue; end
        lickSet = PLOT.probLicksZoom;
        figTag  = 'FIG 2';
    else
        if ~PLOT.fig3, continue; end
        lickSet = PLOT.probLicks;
        figTag  = 'FIG 3';
    end

    nLk = numel(lickSet);

% ---- per-session P(lick) ----
    P_sess = nan(nS, nLk, 2);
    for s = 1:nS
        rC = [SESS(s).length{1}; SESS(s).length{3}];
        rS = [SESS(s).length{2}; SESS(s).length{4}];
        if ~isempty(rC), P_sess(s,:,1) = mean(~isnan(rC(:,lickSet)), 1); end
        if ~isempty(rS), P_sess(s,:,2) = mean(~isnan(rS(:,lickSet)), 1); end
    end

% ---- pooled-trial proportions ----
    poolC = [y1_length; y3_length];
    poolS = [y2_length; y4_length];
    pool  = {poolC, poolS};
    pCol  = {PLOT.colCtrl, PLOT.colStim};
    pName = {'Control','Stim'};

    pPool = nan(2, nLk);
    nTot  = zeros(1,2);
    for g = 1:2
        dat     = pool{g};
        nTot(g) = size(dat,1);
        for li = 1:nLk
            pPool(g,li) = mean(~isnan(dat(:, lickSet(li))));
        end
    end

% ---- statistics ----
    p_prob = nan(1,nLk); q_prob = nan(1,nLk); n_prob = zeros(1,nLk);
    for li = 1:nLk
        x = P_sess(:,li,1); y = P_sess(:,li,2);
        ok = ~isnan(x) & ~isnan(y);
        n_prob(li) = sum(ok);
        if sum(ok) >= 2
            if any(x(ok) - y(ok) ~= 0)
                p_prob(li) = signrank(x(ok), y(ok), 'tail', STATS.tail);
            end
            q_prob(li) = ranksum(x(ok), y(ok), 'tail', STATS.tail);
        end
    end

    pAll = {p_prob, q_prob};
    aAll = cell(size(pAll));
    for z = 1:numel(pAll)
        pv = pAll{z}(:); av = nan(size(pv)); okz = ~isnan(pv);
        if any(okz) && STATS.bh
            [ps, ordz] = sort(pv(okz));
            mm = numel(ps);
            ad = ps .* (mm ./ (1:mm)');
            ad = flipud(cummin(flipud(ad)));
            tmp = nan(mm,1); tmp(ordz) = min(ad,1);
            av(okz) = tmp;
        elseif ~STATS.bh
            av = pv;
        end
        aAll{z} = av';
    end
    pAdjP = aAll{1};
    qAdjP = aAll{2};

    if strcmpi(STATS.figTest,'ranksum'), sP = qAdjP; else, sP = pAdjP; end
    starContacts = lickSet(sP < STATS.alpha);

% ---- draw ----
    figure('Color','w','Units','normalized','Position',[0.10 0.25 0.72 0.50]);
    euName  = {'trial','session'};
    offsets = [-0.1, 0.1];

    for eu = 1   % error bars over TRIALS only

        axP = subplot(1,1,1); hold on;

        pHat = nan(2, nLk);
        pCI  = nan(2, nLk);
        for g = 1:2
            if eu == 2
                M          = P_sess(:,:,g);
                pHat(g,:)  = nanmean(M, 1);
                nn         = sum(~isnan(M), 1);
                tc         = 1.96*ones(size(nn));
                tc(nn > 1) = tinv(1 - STATS.alpha/2, nn(nn > 1) - 1);
                pCI(g,:)   = tc .* nanstd(M, 0, 1) ./ sqrt(max(nn,1));
            else
                tc        = tinv(1 - STATS.alpha/2, max(nTot(g)-1,1));
                pHat(g,:) = pPool(g,:);
                pCI(g,:)  = tc .* sqrt(pPool(g,:).*(1-pPool(g,:))./nTot(g));
            end
        end

        for g = 1:2
            xv = lickSet + offsets(g);
            lo = max(pHat(g,:) - pCI(g,:), 0);
            hi = min(pHat(g,:) + pCI(g,:), 1);
            plot(xv, pHat(g,:), '-', 'Color', pCol{g}, 'LineWidth', 2);
            for li = 1:nLk
                plot([xv(li) xv(li)], [lo(li) hi(li)], '-', 'Color', pCol{g}, 'LineWidth', 2.5);
                plot(xv(li), pHat(g,li), 'o', 'MarkerSize', 8, ...
                     'MarkerFaceColor', pCol{g}, 'MarkerEdgeColor','none');
            end
        end

        if figIdx == 2
            for li = 1:nLk
                if sP(li) < STATS.alpha
                    text(lickSet(li), 1.04, '*', 'FontSize', PLOT.fontSize+4, ...
                         'HorizontalAlignment','center', 'Color','k');
                end
            end
            ylim([0 1.10]);
        else
            ylim([0 1]);
        end

        if eu == 1
            for g = 1:2
                text(0.06, 0.14 - 0.07*(g-1), pName{g}, 'Units','normalized', ...
                     'Color', pCol{g}, 'FontSize', PLOT.fontSize, 'FontWeight','bold');
            end
        end

        if STIM.contact >= min(lickSet) && STIM.contact <= max(lickSet)
            line([STIM.contact STIM.contact], ylim(axP), 'Color', [0.6 0.6 0.6], ...
                 'LineStyle','--', 'LineWidth', 1.5);
        end

        xlim([min(lickSet)-0.5, max(lickSet)+0.5]);
        xticks(lickSet); xtickangle(0);
        xlabel('Contact number'); ylabel('P(lick)');
        title(sprintf('error bars over %ss', euName{eu}), 'FontWeight','normal');
        set(axP, 'TickDir','out', 'FontSize', PLOT.fontSize); box off
    end

    fprintf('\n%s -- P(lick) Control vs Stim, per-session proportions\n', figTag);
    fprintf('  trials pooled: Control %d, Stim %d\n', nTot(1), nTot(2));
    fprintf('  tail = %s | stars use %s', STATS.tail, STATS.figTest);
    if STATS.bh
        fprintf(', BH over m = %d contacts\n', nLk);
    else
        fprintf(', UNCORRECTED\n');
    end
    disp(table(lickSet', pPool(1,:)', pPool(2,:)', p_prob', pAdjP', q_prob', qAdjP', n_prob', ...
        'VariableNames', {'Contact','P_ctrl','P_stim','p','pAdj','q','qAdj','n_sess'}));
    if isempty(starContacts)
        fprintf('  starred contacts: none\n');
    else
        fprintf('  starred contacts: %s\n', num2str(starContacts));
    end

end

%% HELPERS

function sets = findConsecutiveSets(indices, minLength, maxLength)
    sets = {};
    currentSet = [];

    i = 1;
    while i <= length(indices)-1
        if indices(i+1) - indices(i) == 1
            currentSet = [currentSet, indices(i)];
        else
            currentSet = [currentSet, indices(i)];
            while length(currentSet) >= minLength
                truncatedSet = currentSet(1:min(length(currentSet), maxLength));
                sets{end+1} = truncatedSet;
                currentSet = currentSet(min(length(currentSet), maxLength)+1:end);
            end
            currentSet = [];
        end
        i = i + 1;
    end

    currentSet = [currentSet, indices(end)];
    while length(currentSet) >= minLength
        truncatedSet = currentSet(1:min(length(currentSet), maxLength));
        sets{end+1} = truncatedSet;
        currentSet = currentSet(min(length(currentSet), maxLength)+1:end);
    end
end

function [stimTr, ctrlTr] = stimTrialGroups(obj, p, stimContact, task)
% Replaces find_StimTrials (base code\other_codes\functions_td). Returns
% only the two outputs this script used -- S21 (stim) and S21c (control) -- for
% the task == 16 branch, with every expression kept exactly as in find_StimTrials:
%                  otherwise        : index of the lick contact closest to the laser
%                                     of contacts before the go cue; 0 = no laser
    if task ~= 16
        error('stimTrialGroups: only task 16 is implemented (the value this script uses).');
    end

    StimTrials = zeros(1, obj.bp.Ntrials);
    AbsIndices = cell(1, obj.bp.Ntrials);
    lickContact = cell(1, obj.bp.Ntrials);
    GC = zeros(1, obj.bp.Ntrials);

    for i = 1:obj.bp.Ntrials
        temp = obj.sglx.laserTrigIX{i,1};
        if ~isempty(temp)
            temp = temp ./ 25000;
            AbsIndices{i} = temp(1) - 0.5;
        else
            AbsIndices{i} = [];
        end
        lickContact{i} = obj.bp.ev.lickL{i};
        GC(i) = obj.bp.ev.goCue(i);

        if stimContact == 1
            if ~isempty(AbsIndices{i})
                StimTrials(i) = 1;
            else
                StimTrials(i) = 0;
            end
        else
            if ~isempty(AbsIndices{i})
                differences = abs(lickContact{i} - AbsIndices{i});
                [~, closest_index] = min(differences);
                if ~isempty(closest_index)
                    StimTrials(i) = closest_index - sum(lickContact{i} < GC(i));
                else
                    StimTrials(i) = 0;
                end
            end
        end
    end

    Position    = obj.bp.trialTypes;
    LickedOrNot = double(~cellfun(@isempty, obj.bp.ev.lickL))';

    if stimContact == 1
        stimTr = find(StimTrials == stimContact & Position == p);
        ctrlTr = find(StimTrials == 0 & Position == p);
    else
        LorN = 1;
        stimTr = find(StimTrials == stimContact & Position == p & LickedOrNot == LorN );
        ctrlTr = find(StimTrials == 0 & Position == p & LickedOrNot == LorN );
    end
end

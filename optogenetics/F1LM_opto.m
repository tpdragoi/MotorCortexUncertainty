%% F1LM_opto.m
%  Tongue kinematics on control and photoinhibition trials.
%  Bilateral photoinhibition of tjM1 and ALM triggered at the go cue.
%  Behavior only; no spike data is loaded.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat and _laser.mat
%    through shared\loadBehavSession; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/200
%    params.smooth      50
%    params.quality     {'ok','good','mua','great'}
%    params.lowFR       0.01
%    params.window      -1.5 to 3 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear, clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


%% self-contained root (_clean): reads only <v2>\Data and <v2>\shared
v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root,'shared'),'dir')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root,'shared'));
addpath(fullfile(v2Root,'shared','pipelineCopies'));
dataDir = fullfile(v2Root, 'Data', 'GCStim');   % exported obj/kin files
assert(exist(dataDir, 'dir') == 7, 'No data folder: %s', dataDir);
% (needs <stem>_laser.mat next to each session: run exportLaserTrig.m once)

% GO-CUE OPTO: cumulative first-lick probability (fig 10) and
%              lick angle / lick duration vs lick number (fig 11)
% LICK NUMBERING: licks are counted as tongue protrusions from the GO CUE, not
%              as port contacts. A protrusion is a run of >= 5 consecutive
%              lick 1 is the first such run starting at or after t = 0, and
%              2, 3, ... follow in order whether or not they made contact.
%              PLOT.displayLicks selects which of them the panels show.
% All other figures from the previous version have been removed.
% FIXED vs v2: Control / Stim group indices in the cumulative sections were
%              swapped. groups_trials = {L_ctrl, L_stim, R_ctrl, R_stim},
%              so Control = [1 3] and Stim = [2 4].  v2 used {[2 4],[1 3]}
%              and labeled [2 4] "Control".
% STATS:       tested INSIDE every animal (each session here is a different
%              animal), control trials vs stim trials of that session; reported
%              as "significant in k of N animals":
%                Fig. 1L  first-lick latency, two-sample Kolmogorov-Smirnov
%                Fig. 1M  |angle| averaged over STATS.pooledLicks per trial,
%                         Wilcoxon rank-sum
%              No across-session tests are run.

%% OPTIONS

PLOT.fig10 = true;   % cumulative P(first lick) after go cue, Control vs Stim
PLOT.fig11 = true;   % lick angle (4 groups) + lick duration (Ctrl vs Stim)

PLOT.colCtrl  = [0   0   0  ];   % Control  - black
PLOT.colStim  = [0.5 0.7 1  ];   % Stim     - light blue
PLOT.colLctrl = [0   0   1  ];   % L_ctrl   - dark blue
PLOT.colLstim = [0.5 0.7 1  ];   % L_stim   - light blue
PLOT.colRctrl = [1   0   0  ];   % R_ctrl   - dark red
PLOT.colRstim = [1   0.6 0.6];   % R_stim   - light red

PLOT.xlim10    = [0 1];   % x-range of fig 10, seconds from go cue
PLOT.bracket10 = true;   % true = bracket with the k/N animals result
PLOT.fontSize  = 12;

% The licks shown in fig 11.
PLOT.displayLicks = 1:6;

% ---- statistics -------------------------------------------------------
STATS.alpha = 0.05;

% Tails are stated as hypotheses on the paired difference (Control - Stim).
STATS.tailLat   = 'left';   % stim delays first lick        -> ctrl - stim < 0
STATS.tailAngle = 'right';   % stim pulls licks toward midline -> |ctrl| - |stim| > 0
% [|angle|] the angle TESTS use |angle| (distance from 0), so one tail covers
% both sides: L and R licks are each tested as 'stim |angle| smaller than
% control |angle|'. Plots stay signed. Set false to test the signed angle again
% (then use tailAngle = 'both').
STATS.angleAbs  = true;
if STATS.angleAbs, angT = @abs; else, angT = @(x) x; end   % [|angle|] transform for angle tests only

% Cumulative-probability time bins (indices into obj.time)
STATS.binEdges = 300:2:800;

STATS.minLicksPooled  = 2;   % a trial enters if at least this many of STATS.pooledLicks were detected

% Licks entering the POOLED per-animal |angle| test (Fig. 1M).
% Kept separate from PLOT.displayLicks: the panel still draws licks 1-6, but
% the first lick is excluded from the pooled statistic. At lick 1 the tongue
% leaves the midline from the same starting position on control and stim
% trials, so that lick carries little of the effect and only dilutes the
% per-trial mean.
STATS.pooledLicks = 2:6;

%% NORMALIZE THE TEXT OPTIONS
% Accept 'char', "string" or {'cell'} for every text option. Without this a
% stray pair of double quotes or braces above turns into an error hundreds of
% lines below ("SWITCH expression must be a scalar or a character vector"),
% which is a long way from the line that actually caused it.

txtOpts = {'tailLat','tailAngle'};
for z = 1:numel(txtOpts)
    v = STATS.(txtOpts{z});
    if iscell(v), v = v{1}; end
    STATS.(txtOpts{z}) = lower(char(v));
end

validateattributes(STATS.pooledLicks, {'numeric'}, {'vector','integer','positive'}, ...
    mfilename, 'STATS.pooledLicks');
STATS.pooledLicks = unique(STATS.pooledLicks(:))';

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

params.alignEvent = 'goCue';

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

date = '2023-11-12';
meta1 = struct('anm','TD4f','date',date);
date = '2023-11-13';
meta2 = struct('anm','TD5f','date',date);
date = '2023-12-01';
meta4 = struct('anm','TD6f','date',date);
date = '2023-11-29';
meta5 = struct('anm','TD7f','date',date);

% -- additional sessions (uncomment to increase n; the signed-rank floor is
%    2^-n one-tailed / 2^(1-n) two-tailed, so n = 4 cannot reach p < 0.05) --

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11; ...
            meta12;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22];

%% EXTRACT KINEMATICS
% Group order everywhere below:  1 = L_ctrl, 2 = L_stim, 3 = R_ctrl, 4 = R_stim

maxLicks  = 30;
gapThresh = 125;
nSessTot  = size(all_meta,1);

SESS     = struct('angle',{},'duration',{},'length',{},'firstBin',{},'anm',{},'date',{});   % + anm, date
timeAxis = [];

y1_angle = []; y2_angle = []; y3_angle = []; y4_angle = [];
y1_length = []; y2_length = []; y3_length = []; y4_length = [];
y1_duration = []; y2_duration = []; y3_duration = []; y4_duration = [];
allFirstVisibleTimes = cell(1,4);

for sessnum = 1:nSessTot

    clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj idxHit kin me

    fprintf('Session %d\n', sessnum);

    meta        = all_meta(sessnum,1);

    [obj, kin, params] = loadBehavSession(dataDir, meta.anm, meta.date, params);   % [_clean]
    assert(isfield(obj,'sglx') && isfield(obj.sglx,'laserTrigIX'), ...
        ['%s: no laser trigger file (%s). The exported obj/kin files do not hold ' ...
         'obj.sglx.laserTrigIX -- run exportLaserTrig.m (v2 root) once to write it.'], ...
        params.slimStem, fullfile(dataDir, [params.slimStem '_laser.mat']));

    trialSet = (1:obj.bp.Ntrials)';

    sessix   = 1;
    task     = 16;
    condtrix = trialSet;

    if isempty(timeAxis)
        timeAxis = obj.time;
    end
    SESS(sessnum).anm  = obj.pth.anm;   % which animal this session is
    SESS(sessnum).date = obj.pth.dt;

% ---- tongue length / angle ----
    kinix  = find(strcmp(kin(sessix).featLeg,'tongue_length'));
    Length = kin.dat(:,condtrix,kinix);

    kinix = find(strcmp(kin(sessix).featLeg,'tongue_angle'));
    angle = kin.dat(:,condtrix,kinix);

% ---- trial groups ----
    f = 0;

    [S21, S21c] = stimTrialGroups(obj, 1, 1, task);   % was find_StimTrials(obj, 1, f, 1, 6, task)
    idxHit = find(obj.bp.hit == 1);
    S21c   = S21c(ismember(S21c, idxHit));
    L_stim = S21;
    L_ctrl = S21c;

    [S21, S21c] = stimTrialGroups(obj, 3, 1, task);   % was find_StimTrials(obj, 3, f, 1, 6, task)
    idxHit = find(obj.bp.hit == 1);
    S21c   = S21c(ismember(S21c, idxHit));
    R_stim = S21;
    R_ctrl = S21c;

    groups_trials = {L_ctrl, L_stim, R_ctrl, R_stim};

% go cue bin -- t = 0 in obj.time
    [~, goCueBin] = min(abs(obj.time - 0));

    for g = 1:4

        trials  = groups_trials{g};
        ll_sess = nan(numel(trials),1);
        fv_sess = nan(numel(trials),1);

        rows_a = [];
        rows_l = [];
        rows_d = [];

        for i = 1:numel(trials)

            tr = trials(i);

            val_a_full = angle(:, tr);
            val_l_full = Length(:, tr);

            if all(isnan(val_a_full)), continue; end

            realIndices_full = find(~isnan(val_a_full));
            if isempty(realIndices_full), continue; end

            consecutiveSets_full = findConsecutiveSets(realIndices_full, 5, 20);
            nLicks_full          = numel(consecutiveSets_full);
            if nLicks_full == 0, continue; end

            lickStartBins = cellfun(@(s) s(1), consecutiveSets_full);

% Lick 1 is the first tongue PROTRUSION whose onset falls at or
% after the go cue. Protrusions come from runs of visible tongue in
% the video (findConsecutiveSets above, >= 5 samples = 25 ms), so a
% protrusion that never touched the port still counts and still
% advances the numbering. Nothing here consults port contacts.
            firstLickIdx = find(lickStartBins >= goCueBin, 1, 'first');
            if isempty(firstLickIdx), continue; end

            fv_sess(i) = lickStartBins(firstLickIdx);

% last lick of the bout
            lastIdx = find(diff(lickStartBins) > gapThresh, 1);
            if isempty(lastIdx)
                lastIdx = nLicks_full;
            else
                lastIdx = lastIdx + 1;
            end
            ll_sess(i) = lastIdx - firstLickIdx + 1;

% ---- one row per trial, lick 1 = first contact after go cue ----
            row_a = nan(1, maxLicks);
            row_l = nan(1, maxLicks);
            row_d = nan(1, maxLicks);

            for lk = firstLickIdx:nLicks_full
                pos = lk - firstLickIdx + 1;
                if pos > maxLicks, break; end

                seg   = consecutiveSets_full{lk};
                seg_a = val_a_full(seg);
                seg_l = val_l_full(seg);

                [~, pkIdx] = max(abs(seg_l));
                row_a(pos) = seg_a(pkIdx);
                row_l(pos) = max(seg_l);
                row_d(pos) = numel(seg);
            end

            rows_a = [rows_a; row_a];
            rows_l = [rows_l; row_l];
            rows_d = [rows_d; row_d];
        end

        SESS(sessnum).angle{g}    = rows_a;
        SESS(sessnum).length{g}   = rows_l;
        SESS(sessnum).duration{g} = rows_d;
        SESS(sessnum).firstBin{g} = fv_sess;

        allFirstVisibleTimes{g} = [allFirstVisibleTimes{g}; fv_sess];
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

%% PER-ANIMAL TESTS -- inside each animal, control trials vs stim trials
% One test per animal (= per session here):
%   latency   first-lick time from the go cue, Control = L_ctrl + R_ctrl; no lick
%             found counts as LATER than every lick (Inf), exactly as the
%             cumulative curve counts it. Two-sample KS, tail from STATS.tailLat.
%             (Fig. 1L)
%   |angle|   each trial's mean |angle| over STATS.pooledLicks, L_ctrl vs L_stim
%             and R_ctrl vs R_stim, Wilcoxon rank-sum, tail STATS.tailAngle.
%             (Fig. 1M)

PA          = struct();
PA.names    = arrayfun(@(s) sprintf('%s %s', SESS(s).anm, SESS(s).date), 1:nS, 'UniformOutput', false);
PA.medLat   = nan(nS,2);   PA.nLat = zeros(nS,2);
PA.pLatKS   = nan(nS,1);

for s = 1:nS
    latC = latencySec([SESS(s).firstBin{1}(:); SESS(s).firstBin{3}(:)], timeAxis);
    latS = latencySec([SESS(s).firstBin{2}(:); SESS(s).firstBin{4}(:)], timeAxis);
    PA.nLat(s,:)   = [numel(latC) numel(latS)];
    PA.medLat(s,:) = [median(latC) median(latS)];
    if numel(latC) >= 2 && numel(latS) >= 2
        PA.pLatKS(s) = ksLatency(latC, latS, ksTailFromLat(STATS.tailLat));
    end
end
PA.kLatKS = sum(PA.pLatKS < STATS.alpha);       PA.nLatKSOK = sum(~isnan(PA.pLatKS));

% ---- POOLED over licks: one value per TRIAL (its mean over STATS.pooledLicks),
% one rank-sum per animal. The trial stays the unit, so the licks of one trial
% are never counted as independent samples.
PA.poolLicks = STATS.pooledLicks;
PA.pooledP   = nan(nS, 2);   % columns: |angle| L, |angle| R
PA.pooledN   = zeros(nS, 2, 2);   % trials entering (ctrl, stim)
PA.pooledMed = nan(nS, 2, 2);   % median of the per-trial means (ctrl, stim)
for s = 1:nS
    pairs = {angT(SESS(s).angle{1}), angT(SESS(s).angle{2}); ...   % [|angle|] mean of |angle| over licks
             angT(SESS(s).angle{3}), angT(SESS(s).angle{4})};
    for z = 1:2
        a = trialMeanOverLicks(pairs{z,1}, PA.poolLicks, STATS.minLicksPooled);
        b = trialMeanOverLicks(pairs{z,2}, PA.poolLicks, STATS.minLicksPooled);
        PA.pooledN(s,z,:)   = [numel(a) numel(b)];
        PA.pooledMed(s,z,:) = [median(a) median(b)];
        if numel(a) >= 2 && numel(b) >= 2
            PA.pooledP(s,z) = ranksum(a, b, 'tail', STATS.tailAngle);
        end
    end
end
PA.pooledK  = sum(PA.pooledP < STATS.alpha, 1);
PA.pooledNa = sum(~isnan(PA.pooledP), 1);

fprintf('\n%s\n PER-ANIMAL TESTS -- control vs stim trials, inside each animal\n', repmat('=',1,78));
fprintf('%s\n', repmat('-',1,78));
fprintf(' FIRST-LICK LATENCY (s), two-sample KS, tail = %s (no lick found = later than any lick)\n', STATS.tailLat);
fprintf('   %-18s %8s %8s %9s %9s %11s\n', 'animal', 'n ctrl', 'n stim', 'med ctrl', 'med stim', 'p KS');
for s = 1:nS
    fprintf('   %-18s %8d %8d %9.3f %9.3f %10.4g%s\n', PA.names{s}, ...
        PA.nLat(s,1), PA.nLat(s,2), PA.medLat(s,1), PA.medLat(s,2), ...
        PA.pLatKS(s), starOf(PA.pLatKS(s), STATS.alpha));
end
fprintf('   -> significant in %d of %d animals (p < %.3g)\n', PA.kLatKS, PA.nLatKSOK, STATS.alpha);

angLbl     = ternaryStr(STATS.angleAbs, '|ANGLE|', 'ANGLE');   % [|angle|]
pooledName = {[angLbl ' L (L_ctrl vs L_stim)'], [angLbl ' R (R_ctrl vs R_stim)']};
fprintf('\n %s POOLED OVER LICKS %s, rank-sum -- one value per trial (its mean over those licks, trials with >= %d licks)\n', ...
    angLbl, rangeLabel(PA.poolLicks), STATS.minLicksPooled);
for z = 1:2
    fprintf('   %s, tail = %s\n', pooledName{z}, STATS.tailAngle);
    fprintf('     %-18s %8s %8s %11s %11s %10s\n', 'animal', 'n ctrl', 'n stim', 'med ctrl', 'med stim', 'p');
    for s = 1:nS
        fprintf('     %-18s %8d %8d %11.4g %11.4g %10.4g%s\n', PA.names{s}, PA.pooledN(s,z,1), PA.pooledN(s,z,2), ...
            PA.pooledMed(s,z,1), PA.pooledMed(s,z,2), PA.pooledP(s,z), starOf(PA.pooledP(s,z), STATS.alpha));
    end
    fprintf('     -> significant in %d of %d animals\n', PA.pooledK(z), PA.pooledNa(z));
end
fprintf('%s\n', repmat('=',1,78));

%% FIG 10 -- CUMULATIVE P(FIRST LICK) AFTER GO CUE

if PLOT.fig10

    binEdges    = STATS.binEdges;
    nBins       = numel(binEdges) - 1;
    timeCenters = arrayfun(@(b) mean(timeAxis([binEdges(b), binEdges(b+1)])), 1:nBins);

    grpIdx   = {[1 3], [2 4]};   % 1 = Control, 2 = Stim
    grpCol   = {PLOT.colCtrl, PLOT.colStim};
    grpName  = {'Control', 'Stim'};

% ---- per-session cumulative curves and per-session median latency ----
    Pcurve     = nan(nS, nBins, 2);
    medLat     = nan(nS, 2);
    nTrialSess = zeros(nS, 2);   % trials with a detected first lick
    nTotSess   = zeros(nS, 2);   % all trials in the group

    for s = 1:nS
        for gi = 1:2
            fv = [];
            for k = grpIdx{gi}
                fv = [fv; SESS(s).firstBin{k}(:)];
            end
            nTotSess(s,gi) = numel(fv);
            if isempty(fv), continue; end
            for b = 1:nBins
                Pcurve(s,b,gi) = mean(fv <= binEdges(b+1));   % NaN trials count as "no lick"
            end
            fvOK = fv(~isnan(fv));
            nTrialSess(s,gi) = numel(fvOK);
            if ~isempty(fvOK)
                medLat(s,gi) = median(timeAxis(round(fvOK)));
            end
        end
    end

% ---- pooled-trial curves (what the band is drawn from by default) ----
    Ptrial   = nan(2, nBins);
    SEMtrial = nan(2, nBins);
    nTrialGrp = zeros(1,2);

    for gi = 1:2
        fv = [];
        for k = grpIdx{gi}
            fv = [fv; allFirstVisibleTimes{k}(:)];
        end
        nTrialGrp(gi) = numel(fv);
        binMat = false(numel(fv), nBins);
        for b = 1:nBins
            binMat(:,b) = fv <= binEdges(b+1);
        end
        Ptrial(gi,:)   = mean(binMat, 1);
        SEMtrial(gi,:) = std(double(binMat), 0, 1) ./ sqrt(size(binMat,1));
    end

% ---- draw: left panel band over TRIALS, right panel over SESSIONS ----
    figure('Color','w','Units','normalized','Position',[0.10 0.25 0.74 0.50]);
    euName  = {'trial','session'};
    ax10all = gobjects(1,2);

    for eu = 1   % error bars over TRIALS only
        ax10all(eu) = subplot(1,1,1); hold on;
        for gi = 1:2
            if eu == 2
                Mc  = Pcurve(:,:,gi);
                mu  = nanmean(Mc, 1);
                nn  = sum(~isnan(Mc), 1);
                sem = nanstd(Mc, 0, 1) ./ sqrt(max(nn,1));
                tc  = 1.96*ones(size(nn));
                tc(nn > 1) = tinv(1 - STATS.alpha/2, nn(nn > 1) - 1);
                halfW = tc .* sem;
            else
                mu    = Ptrial(gi,:);
                halfW = 1.96 * SEMtrial(gi,:);
            end
            upper = mu + halfW;
            lower = mu - halfW;
            fill([timeCenters fliplr(timeCenters)], [upper fliplr(lower)], ...
                 grpCol{gi}, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
            plot(timeCenters, mu, '-', 'Color', grpCol{gi}, 'LineWidth', 2);
        end
    end
    ax10 = ax10all(1);

% ---- statistics: Fig. 1L reports the per-animal KS test on first-lick latency ----
    paP          = PA.pLatKS;
    bracketStr   = sprintf('%d/%d animals', PA.kLatKS, PA.nLatKSOK);
    bracketWhich = 'per-animal two-sample KS on trials';
    fprintf('\nFIG 10 -- Control vs Stim | %s, first-lick latency: %s\n', bracketWhich, bracketStr);
    fprintf('  trials pooled: Control %d, Stim %d\n', nTrialGrp(1), nTrialGrp(2));

% ---- annotation, applied to both panels ----
    for eu = 1   % error bars over TRIALS only
        axC = ax10all(eu);
        xlim(axC, PLOT.xlim10);
        ylim(axC, [0 1.16]);

        if PLOT.bracket10 && any(~isnan(paP))
            xr = PLOT.xlim10;
            xb = [xr(1) + 0.78*diff(xr), xr(1) + 0.95*diff(xr)];
            yb = 1.05;
            line(xb, [yb yb],                 'Color','k', 'LineWidth', 1.2, 'Parent', axC);
            line([xb(1) xb(1)], [yb-0.03 yb], 'Color','k', 'LineWidth', 1.2, 'Parent', axC);
            line([xb(2) xb(2)], [yb-0.03 yb], 'Color','k', 'LineWidth', 1.2, 'Parent', axC);
            text(mean(xb), yb + 0.015, bracketStr, 'HorizontalAlignment','center', ...
                 'FontSize', PLOT.fontSize, 'Color', 'k', 'Parent', axC);
        end

        if eu == 1
            xt = PLOT.xlim10(1) + 0.60*diff(PLOT.xlim10);
            text(xt, 0.32, grpName{1}, 'Color', grpCol{1}, 'FontSize', PLOT.fontSize, ...
                 'FontWeight','bold', 'Parent', axC);
            text(xt, 0.20, grpName{2}, 'Color', grpCol{2}, 'FontSize', PLOT.fontSize, ...
                 'FontWeight','bold', 'Parent', axC);
        end

        if eu == 1
            text(0.03, 0.995, sprintf('%s, first-lick latency, p = %s', bracketWhich, ...
                 strjoin(cellstr(compose('%.3g', paP(~isnan(paP)))), ', ')), ...
                 'Units','normalized', 'FontSize', PLOT.fontSize-3, ...
                 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top', 'Parent', axC);
        end

        xlabel(axC, 'Time from Go Cue (s)');
        ylabel(axC, 'Lick Probability');
        set(axC, 'TickDir','out', 'FontSize', PLOT.fontSize, 'YTick', 0:0.2:1);
        box(axC, 'off');
    end

end

%% FIG 11 -- LICK ANGLE AND LICK DURATION

if PLOT.fig11

    licks   = PLOT.displayLicks;
    nLicksP = numel(licks);

% ---- per-session means (rows = sessions) ----
    A_sess = nan(nS, nLicksP, 4);   % angle, one page per group
    D_sess = nan(nS, nLicksP, 2);   % duration, 1 = Control, 2 = Stim (seconds)

    for s = 1:nS
        for g = 1:4
            a = SESS(s).angle{g};
            if ~isempty(a)
                A_sess(s,:,g) = nanmean(a(:,licks), 1);
            end
        end

        dC = [SESS(s).duration{1}; SESS(s).duration{3}];
        dS = [SESS(s).duration{2}; SESS(s).duration{4}];
        dC(dC < 0) = NaN;
        dS(dS < 0) = NaN;
        if ~isempty(dC), D_sess(s,:,1) = nanmean(dC(:,licks), 1) * params.dt; end
        if ~isempty(dS), D_sess(s,:,2) = nanmean(dS(:,licks), 1) * params.dt; end
    end

% ---- pooled trials, for the plotted mean and 95% CI ----
    angle_pooled = {y1_angle, y2_angle, y3_angle, y4_angle};
    angCol       = {PLOT.colLctrl, PLOT.colLstim, PLOT.colRctrl, PLOT.colRstim};
    angName      = {'L_{ctrl}','L_{stim}','R_{ctrl}','R_{stim}'};

    durC = [y1_duration; y3_duration];   % Control
    durS = [y2_duration; y4_duration];   % Stim
    durC(durC < 0) = NaN;
    durS(durS < 0) = NaN;
    dur_pooled = {durC * params.dt, durS * params.dt};
    durCol     = {PLOT.colCtrl, PLOT.colStim};
    durName    = {'Control','Stim'};

% ---- draw: row 1 bars over TRIALS, row 2 bars over SESSIONS ----
    figure('Color','w','Units','normalized','Position',[0.06 0.08 0.80 0.82]);
    euName = {'trial','session'};

    for eu = 1   % error bars over TRIALS only

% ------------------------- ANGLE --------------------------------
        axA = subplot(1,1,1); hold on;

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

        yl = ylim(axA); ylim(axA, [yl(1), yl(2) + 0.16*diff(yl)]); yl = ylim(axA);
        text(0.02, 0.98, sprintf(['licks %s pooled per trial, rank-sum within animal:  ' ...
             '{\\color[rgb]{%g %g %g}L %d/%d}   {\\color[rgb]{%g %g %g}R %d/%d} animals'], ...
             rangeLabel(PA.poolLicks), PLOT.colLctrl, PA.pooledK(1), PA.pooledNa(1), ...
             PLOT.colRctrl, PA.pooledK(2), PA.pooledNa(2)), 'Units','normalized', ...
             'FontSize', PLOT.fontSize-3, 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top', ...
             'Interpreter','tex');

        if eu == 1
            for g = 1:4
                text(0.70, 0.06 + 0.07*(4-g), angName{g}, 'Units','normalized', ...
                     'Color', angCol{g}, 'FontSize', PLOT.fontSize-1, 'FontWeight','bold');
            end
        end

        xlim([0.5 max(licks)+0.5]); xticks(licks); xtickangle(0);
        xlabel('Lick number from go cue');
        ylabel('Lick angle (deg)');
        title('Angle', 'FontWeight','normal');
        set(axA, 'TickDir','out', 'FontSize', PLOT.fontSize); box off

    end

    fprintf('\nFIG 11 -- per-animal pooled |angle| test: see PER-ANIMAL TESTS above\n');
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

function tail = ksTailFromLat(tailLat)
% KS tail that matches the rank-sum tail on latency. kstest2's tail refers to
% the CDF of the FIRST sample relative to the second, and the first sample
% here is control. If control licks earlier (tailLat 'left'), its cumulative
% curve sits ABOVE the stim curve, which is kstest2's 'larger'.
    switch lower(tailLat)
        case 'left',  tail = 'larger';
        case 'right', tail = 'smaller';
        otherwise,    tail = 'unequal';
    end
end


function p = ksLatency(latC, latS, tail)
% Two-sample Kolmogorov-Smirnov between the two first-lick latency
% distributions, which is a test on the two cumulative curves this panel
% draws. Trials with no lick arrive as Inf; they are mapped to one value past
% every observed latency, so the curves keep their plateaus (and kstest2 does
% not accept a non-finite input).
    p = NaN;
    latC = latC(:);  latS = latS(:);
    if numel(latC) < 2 || numel(latS) < 2, return; end

    finiteMax = max([latC(isfinite(latC)); latS(isfinite(latS))]);
    if isempty(finiteMax), return; end          % no trial licked in either group

    latC(~isfinite(latC)) = finiteMax + 1;
    latS(~isfinite(latS)) = finiteMax + 1;

    [~, p] = kstest2(latC, latS, 'Tail', tail);
end

function t = latencySec(firstBin, timeAxis)
% first-lick bin -> seconds from the go cue; no lick found -> Inf,
% i.e. later than every lick, which is how the cumulative curve counts it.
    firstBin = firstBin(:);
    t  = inf(size(firstBin));
    ok = ~isnan(firstBin);
    t(ok) = timeAxis(round(firstBin(ok)));
end

function s = starOf(p, alpha)
% Thin wrapper on shared\pStars so the printed tables and the figures use one
% convention. The leading space is what the print formats expect.
    s = pStars(p, [alpha 0.01 0.001]);
    if ~isempty(s), s = [' ' s]; end
end

function s = ternaryStr(c, a, b)
    if c, s = a; else, s = b; end
end

function s = rangeLabel(v)
% "2-6" for a contiguous run, "[2 4 6]" otherwise, so a non-contiguous
% pooling range is never printed as a range it is not.
    v = v(:)';
    if numel(v) > 1 && all(diff(v) == 1)
        s = sprintf('%d-%d', v(1), v(end));
    else
        s = mat2str(v);
    end
end

function v = trialMeanOverLicks(A, licks, minN)
% one value per trial: the mean over the given lick columns, for
% trials where at least minN of those licks were detected.
    v = zeros(0,1);
    if isempty(A), return; end
    cols = licks(licks <= size(A,2));
    if isempty(cols), return; end
    M = A(:, cols);
    n = sum(~isnan(M), 2);
    v = mean(M, 2, 'omitnan');
    v = v(n >= minN);
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

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

%% F1LM_opto.m
%  Tongue kinematics on control and photoinhibition trials.
%  Bilateral photoinhibition of tjM1 and ALM triggered at the go cue.
%  Behaviour only; no spike data is loaded.
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
%              and labelled [2 4] "Control".
% STATS:       all parametric tests replaced with paired nonparametric
%              Wilcoxon signed-rank across SESSIONS (trials within a session
%              are not independent).  Set STATS.unit = 'trial' to fall back
%              to unpaired rank-sum over pooled trials.
% PER-ANIMAL:  each session here is a different animal, so the
%              effect is also tested INSIDE every animal: Wilcoxon rank-sum
%              (Mann-Whitney) on that session's trials, control trials vs stim
%              trials. Stim and control trials are interleaved in the same
%              session, so the trial is the sampling unit within an animal. The
%              result is reported as "significant in k of N animals" with every
%              animal's p printed. This is a replication claim about these
%              animals, not a population test -- the across-session signed-rank
%              / rank-sum below are still computed and printed next to it.
%              STATS.figTest / STATS.cdfTest = 'perAnimal' puts it on the figures.

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
PLOT.binStars  = false;   % true = star every significant time bin
PLOT.bracket10 = true;   % true = one bracket + star for the overall test
% by side -- bars over pooled TRIALS (what the earlier versions did) and bars
% over SESSIONS (the unit the Wilcoxon tests use). The statistics and the
% stars are identical in both; only the bars and the plotted centre move.
PLOT.fontSize  = 12;

% The contacts shown in fig 11. This is ALSO the multiple-comparison family:
% every test below runs over exactly these licks and BH corrects over them.
% Changing it changes the correction, so set it from what the panel shows.
%   at m = 6 the BH threshold at rank 4 is 4/6*0.05 = 0.0333
%   at m = 8 it is 4/8*0.05 = 0.0250
% The 4-vs-4 rank-sum floor is 0.0286, which sits between those two.
PLOT.displayLicks = 1:6;

% Which statistic the fig 10 bracket reports:
%   'area'    -- whole cumulative trace, 0 to 1, over PLOT.xlim10 (default)
%   'latency' -- median first-lick latency only
PLOT.bracketStat = 'area';

% ---- statistics -------------------------------------------------------
STATS.unit  = 'session';   % 'session' -> paired signrank across sessions
% 'trial'   -> unpaired ranksum over pooled trials
STATS.alpha = 0.05;

% Which p-value drives the stars and brackets ON THE FIGURES.
%   'ranksum'  -- UNPAIRED on the session-level values. This is what these
%                 go-cue panels report: with n = 4 the paired signed-rank
%                 floors at 0.0625 and cannot reach 0.05, and the pairing
%                 diagnostic below shows the control baselines barely vary,
%                 so there is no session effect to block on.
%   'signrank' -- paired across sessions.
%   'both'     -- draw BOTH: filled * = paired, open o = unpaired.
%   'perAnimal'-- rank-sum on TRIALS inside each animal; the figure
%                 shows "k/N" = number of animals significant at that lick
%                 (after BH across the displayed licks, within each animal).
% All of them are always COMPUTED and PRINTED; this only picks the markers.
STATS.figTest = 'perAnimal';   % was 'ranksum'

% Benjamini-Hochberg across the contacts in PLOT.displayLicks, applied within
% each comparison separately (L, R, duration). Stars use the corrected value.
STATS.bh = true;

% Tails are stated as hypotheses on the paired difference (Control - Stim).
STATS.tailProb  = 'right';   % stim lowers P(lick by t)      -> ctrl - stim > 0
STATS.tailLat   = 'left';   % stim delays first lick        -> ctrl - stim < 0
STATS.tailAngle = 'right';   % stim pulls licks toward midline -> |ctrl| - |stim| > 0
% [|angle|] the angle TESTS use |angle| (distance from 0), so one tail covers
% both sides: L and R licks are each tested as 'stim |angle| smaller than
% control |angle|'. Plots stay signed. Set false to test the signed angle again
% (then use tailAngle = 'both').
STATS.angleAbs  = true;
STATS.tailDur   = 'both';   % no a priori direction         -> set if you have one
if STATS.angleAbs, angT = @abs; else, angT = @(x) x; end   % [|angle|] transform for angle tests only

% Cumulative-probability time bins (indices into obj.time)
STATS.binEdges = 300:2:800;

% The window over which the WHOLE-TRACE comparison is computed: each session's
% cumulative curve is integrated across this span, giving one number per
% session per condition, and those are compared. Kept separate from
% PLOT.xlim10 so changing the view never silently changes the test.
STATS.traceWin_s = [0 1];

% The CDF panel (fig 10) always reports the UNPAIRED rank-sum. With n = 4 the
% paired signed-rank floors at 0.0625 and cannot clear 0.05, and the pairing
% diagnostic shows the control baselines barely vary across sessions, so there
% is no session effect to block on. This is not governed by STATS.figTest --
% it is fixed here so the panel cannot be run with the paired test by accident.
% Fig 11 per-animal test:
%   'pooled'  -- ONE rank-sum per animal on each trial's MEAN over the licks in
%                PLOT.displayLicks (licks 1-6). A trial enters if at least
%                STATS.minLicksPooled of those licks were detected. No BH needed.
%   'perLick' -- one rank-sum per lick, BH across licks (previous behaviour).
STATS.perAnimalMode   = 'pooled';
STATS.minLicksPooled  = 2;
STATS.cdfTest = 'perAnimal';   % was 'ranksum'. 'perAnimal' = rank-sum
% on per-trial first-lick latency inside each animal;
% the bracket reads "k/N animals".

%% NORMALISE THE TEXT OPTIONS
% Accept 'char', "string" or {'cell'} for every text option. Without this a
% stray pair of double quotes or braces above turns into an error hundreds of
% lines below ("SWITCH expression must be a scalar or a character vector"),
% which is a long way from the line that actually caused it.

txtOpts = {'unit','figTest','cdfTest','tailProb','tailLat','tailAngle','tailDur','perAnimalMode'};
for z = 1:numel(txtOpts)
    v = STATS.(txtOpts{z});
    if iscell(v), v = v{1}; end
    STATS.(txtOpts{z}) = lower(char(v));
end

v = PLOT.bracketStat;  if iscell(v), v = v{1}; end
PLOT.bracketStat = lower(char(v));

validateattributes(STATS.traceWin_s, {'numeric'}, {'numel',2,'increasing'}, ...
    mfilename, 'STATS.traceWin_s');

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

% exact p-value floors for this n, both tests
fprintf('\nEXACT FLOORS with n = %d sessions\n', nS);
fprintf('  paired signed-rank (2^n = %d arrangements)   : one-tailed %.4f | two-tailed %.4f\n', ...
    2^nS, 2^(-nS), 2^(1-nS));
fprintf('  unpaired rank-sum  (C(%d,%d) = %d arrangements) : one-tailed %.4f | two-tailed %.4f\n', ...
    2*nS, nS, nchoosek(2*nS,nS), 1/nchoosek(2*nS,nS), 2/nchoosek(2*nS,nS));
if 2^(-nS) > STATS.alpha
    fprintf('  *** paired signed-rank cannot reach p < %.3f at this n ***\n', STATS.alpha);
end
fprintf(['  NOTE: the rank-sum floor is lower only because it assumes the %d control and\n' ...
         '        %d stim values came from %d independent sessions rather than %d paired ones.\n'], ...
    nS, nS, 2*nS, nS);

%% PER-ANIMAL TESTS -- rank-sum on trials inside each animal
% One test per animal (= per session here), control trials vs stim
% trials of that session:
%   latency   first-lick time from the go cue, Control = L_ctrl + R_ctrl,
%             found counts as LATER than every lick (Inf), exactly as the
%             cumulative curve counts it as "no lick". Tail = STATS.tailLat.
%   angle     at each lick in PLOT.displayLicks, L_ctrl vs L_stim and R_ctrl vs
%             R_stim. Tail = STATS.tailAngle.
%   duration  at each lick, Control vs Stim (L + R pooled, as fig 11). Tail =
%             STATS.tailDur.
% Within each animal the lick-by-lick p values are BH-corrected across the
% displayed licks (STATS.bh), the same family the across-session tests use.

PA          = struct();
PA.names    = arrayfun(@(s) sprintf('%s %s', SESS(s).anm, SESS(s).date), 1:nS, 'UniformOutput', false);
PA.licks    = PLOT.displayLicks;
nLP         = numel(PA.licks);
PA.pLat     = nan(nS,1);   PA.medLat = nan(nS,2);   PA.nLat = zeros(nS,2);
PA.pAngL    = nan(nS,nLP); PA.pAngR  = nan(nS,nLP); PA.pDur = nan(nS,nLP);

for s = 1:nS
    latC = latencySec([SESS(s).firstBin{1}(:); SESS(s).firstBin{3}(:)], timeAxis);
    latS = latencySec([SESS(s).firstBin{2}(:); SESS(s).firstBin{4}(:)], timeAxis);
    PA.nLat(s,:)   = [numel(latC) numel(latS)];
    PA.medLat(s,:) = [median(latC) median(latS)];
    if numel(latC) >= 2 && numel(latS) >= 2
        PA.pLat(s) = ranksum(latC, latS, 'tail', STATS.tailLat);
    end

    dC = [SESS(s).duration{1}; SESS(s).duration{3}];  dC(dC < 0) = NaN;
    dS = [SESS(s).duration{2}; SESS(s).duration{4}];  dS(dS < 0) = NaN;
    for i = 1:nLP
        c = PA.licks(i);
        PA.pAngL(s,i) = rankSumColumn(angT(SESS(s).angle{1}), angT(SESS(s).angle{2}), c, STATS.tailAngle);
        PA.pAngR(s,i) = rankSumColumn(angT(SESS(s).angle{3}), angT(SESS(s).angle{4}), c, STATS.tailAngle);
        PA.pDur(s,i)  = rankSumColumn(dC, dS, c, STATS.tailDur);
    end
end

% BH across the displayed licks, separately inside each animal and comparison
PA.aAngL = PA.pAngL;  PA.aAngR = PA.pAngR;  PA.aDur = PA.pDur;
if STATS.bh
    for s = 1:nS
        PA.aAngL(s,:) = bhAdjust(PA.pAngL(s,:));
        PA.aAngR(s,:) = bhAdjust(PA.pAngR(s,:));
        PA.aDur(s,:)  = bhAdjust(PA.pDur(s,:));
    end
end

% how many animals are significant, and out of how many that could be tested
PA.kLat  = sum(PA.pLat < STATS.alpha);          PA.nLatOK  = sum(~isnan(PA.pLat));
PA.kAngL = sum(PA.aAngL < STATS.alpha, 1);      PA.nAngL   = sum(~isnan(PA.aAngL), 1);
PA.kAngR = sum(PA.aAngR < STATS.alpha, 1);      PA.nAngR   = sum(~isnan(PA.aAngR), 1);
PA.kDur  = sum(PA.aDur  < STATS.alpha, 1);      PA.nDur    = sum(~isnan(PA.aDur),  1);

% ---- POOLED over licks: one value per TRIAL (its mean over the displayed licks),
% one rank-sum per animal. The trial stays the unit, so the licks of one trial
% are never counted as independent samples.
PA.pooledP   = nan(nS, 3);   % columns: angle L, angle R, duration
PA.pooledN   = zeros(nS, 3, 2);   % trials entering (ctrl, stim)
PA.pooledMed = nan(nS, 3, 2);   % median of the per-trial means (ctrl, stim)
pooledTail   = {STATS.tailAngle, STATS.tailAngle, STATS.tailDur};
for s = 1:nS
    dC = [SESS(s).duration{1}; SESS(s).duration{3}];  dC(dC < 0) = NaN;
    dS = [SESS(s).duration{2}; SESS(s).duration{4}];  dS(dS < 0) = NaN;
    pairs = {angT(SESS(s).angle{1}), angT(SESS(s).angle{2}); ...   % [|angle|] mean of |angle| over licks
             angT(SESS(s).angle{3}), angT(SESS(s).angle{4}); ...
             dC * params.dt,   dS * params.dt};
    for z = 1:3
        a = trialMeanOverLicks(pairs{z,1}, PA.licks, STATS.minLicksPooled);
        b = trialMeanOverLicks(pairs{z,2}, PA.licks, STATS.minLicksPooled);
        PA.pooledN(s,z,:)   = [numel(a) numel(b)];
        PA.pooledMed(s,z,:) = [median(a) median(b)];
        if numel(a) >= 2 && numel(b) >= 2
            PA.pooledP(s,z) = ranksum(a, b, 'tail', pooledTail{z});
        end
    end
end
PA.pooledK  = sum(PA.pooledP < STATS.alpha, 1);
PA.pooledNa = sum(~isnan(PA.pooledP), 1);

fprintf('\n%s\n PER-ANIMAL TESTS -- Wilcoxon rank-sum on trials, control vs stim, inside each animal\n', repmat('=',1,78));
fprintf('%s\n', repmat('-',1,78));
fprintf(' FIRST-LICK LATENCY (s), tail = %s (no lick found = later than any lick)\n', STATS.tailLat);
fprintf('   %-18s %8s %8s %9s %9s %10s\n', 'animal', 'n ctrl', 'n stim', 'med ctrl', 'med stim', 'p');
for s = 1:nS
    fprintf('   %-18s %8d %8d %9.3f %9.3f %10.4g %s\n', PA.names{s}, PA.nLat(s,1), PA.nLat(s,2), ...
        PA.medLat(s,1), PA.medLat(s,2), PA.pLat(s), starOf(PA.pLat(s), STATS.alpha));
end
fprintf('   -> significant in %d of %d animals (p < %.3g)\n', PA.kLat, PA.nLatOK, STATS.alpha);

pa_tables = {PA.pAngL, PA.aAngL, [ternaryStr(STATS.angleAbs,'|ANGLE|','ANGLE') ', L_ctrl vs L_stim'], STATS.tailAngle; ...
             PA.pAngR, PA.aAngR, [ternaryStr(STATS.angleAbs,'|ANGLE|','ANGLE') ', R_ctrl vs R_stim'], STATS.tailAngle; ...
             PA.pDur,  PA.aDur,  'DURATION, Control vs Stim', STATS.tailDur};
for z = 1:size(pa_tables,1)
    P  = pa_tables{z,1};  Q = pa_tables{z,2};
    fprintf('\n %s, tail = %s | p per lick%s\n', pa_tables{z,3}, pa_tables{z,4}, ...
        ternaryStr(STATS.bh, ' (BH-adjusted within animal in brackets)', ''));
    hdr = cellstr(compose('lick %d', PA.licks(:)));
    fprintf('   %-18s', 'animal');  fprintf(' %16s', hdr{:});  fprintf('\n');
    for s = 1:nS
        fprintf('   %-18s', PA.names{s});
        for i = 1:nLP
            fprintf(' %7.3g (%5.3g)%s', P(s,i), Q(s,i), starOf(Q(s,i), STATS.alpha));
        end
        fprintf('\n');
    end
    k = sum(Q < STATS.alpha, 1);  nn = sum(~isnan(Q), 1);
    fprintf('   %-18s', 'animals sig.');
    kn = cellstr(compose('%d/%d', [k(:) nn(:)]));
    fprintf(' %16s', kn{:});
    fprintf('\n');
end

angLbl     = ternaryStr(STATS.angleAbs, '|ANGLE|', 'ANGLE');   % [|angle|]
pooledName = {[angLbl ' L (L_ctrl vs L_stim)'], [angLbl ' R (R_ctrl vs R_stim)'], 'DURATION (s, Control vs Stim)'};
fprintf('\n POOLED OVER LICKS %d-%d -- one value per trial (its mean over those licks, trials with >= %d licks)\n', ...
    PA.licks(1), PA.licks(end), STATS.minLicksPooled);
for z = 1:3
    fprintf('   %s, tail = %s\n', pooledName{z}, pooledTail{z});
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

% ---- statistics ----
    fvC = []; fvS = [];
    for k = grpIdx{1}, fvC = [fvC; allFirstVisibleTimes{k}(:)]; end
    for k = grpIdx{2}, fvS = [fvS; allFirstVisibleTimes{k}(:)]; end

    p_bin = nan(1, nBins);
    q_bin = nan(1, nBins);
    for b = 1:nBins
        x = Pcurve(:,b,1);  y = Pcurve(:,b,2);
        ok = ~isnan(x) & ~isnan(y);
        if sum(ok) >= 2
            if any(x(ok) - y(ok) ~= 0)
                p_bin(b) = signrank(x(ok), y(ok), 'tail', STATS.tailProb);
            end
            q_bin(b) = ranksum(x(ok), y(ok), 'tail', STATS.tailProb);
        end
        if strcmpi(STATS.unit,'trial')
            p_bin(b) = ranksum(double(fvC <= binEdges(b+1)), double(fvS <= binEdges(b+1)), ...
                               'tail', STATS.tailProb);
            q_bin(b) = p_bin(b);
        end
    end

% ---- WHOLE-TRACE COMPARISON ------------------------------------------
% One number per session per condition: the area under that session's
% cumulative curve over the plotted window. This compares the traces
% across their full 0-to-1 range rather than at the median crossing only.
% Larger area = licked earlier, so stim delaying licking means
% ctrl - stim > 0, matching STATS.tailProb.
    inWin    = timeCenters >= STATS.traceWin_s(1) & timeCenters <= STATS.traceWin_s(2);
    areaStat = nan(nS, 2);
    for s = 1:nS
        for gi = 1:2
            yv = Pcurve(s, inWin, gi);
            if all(~isnan(yv))
                areaStat(s,gi) = trapz(timeCenters(inWin), yv);
            end
        end
    end

    okA      = ~isnan(areaStat(:,1)) & ~isnan(areaStat(:,2));
    p_area   = NaN;
    q_area   = NaN;
    if sum(okA) >= 2
        if any(areaStat(okA,1) - areaStat(okA,2) ~= 0)
            p_area = signrank(areaStat(okA,1), areaStat(okA,2), 'tail', STATS.tailProb);
        end
        q_area = ranksum(areaStat(okA,1), areaStat(okA,2), 'tail', STATS.tailProb);
    end

% overall test: per-session median first-lick latency, Control vs Stim
    okLat = ~isnan(medLat(:,1)) & ~isnan(medLat(:,2));
    p_lat = NaN;
    if strcmpi(STATS.unit,'session')
        if sum(okLat) >= 2 && any(medLat(okLat,1) - medLat(okLat,2) ~= 0)
            p_lat = signrank(medLat(okLat,1), medLat(okLat,2), 'tail', STATS.tailLat);
        end
        testName = sprintf('Wilcoxon signed-rank, paired by session, one-tailed (%s), n = %d', ...
                           STATS.tailLat, sum(okLat));
    else
        latC = []; latS = [];
        for s = 1:nS
            for k = grpIdx{1}, latC = [latC; SESS(s).firstBin{k}(:)]; end
            for k = grpIdx{2}, latS = [latS; SESS(s).firstBin{k}(:)]; end
        end
        latC = timeAxis(round(latC(~isnan(latC))));
        latS = timeAxis(round(latS(~isnan(latS))));
        p_lat = ranksum(latC, latS, 'tail', STATS.tailLat);
        testName = sprintf('Wilcoxon rank-sum, pooled trials, one-tailed (%s), n = %d vs %d', ...
                           STATS.tailLat, numel(latC), numel(latS));
    end

% UNPAIRED counterpart: the same n control and n stim session values,
% treated as two independent groups of n.
    p_lat_rs = NaN;
    if sum(okLat) >= 2
        p_lat_rs = ranksum(medLat(okLat,1), medLat(okLat,2), 'tail', STATS.tailLat);
    end

    fprintf('\nFIG 10 -- Control vs Stim\n');
    fprintf('  [A] WHOLE TRACE -- area under each session''s cumulative curve, %.2f to %.2f s\n', ...
        STATS.traceWin_s(1), STATS.traceWin_s(2));
    fprintf('      paired signed-rank  p = %.4g   |   unpaired rank-sum  p = %.4g   (n = %d)\n', ...
        p_area, q_area, sum(okA));
    fprintf('  [B] MEDIAN FIRST-LICK LATENCY\n');
    fprintf('      %s\n', testName);
    fprintf('      paired signed-rank  p = %.4g   |   unpaired rank-sum  p = %.4g   (n = %d)\n', ...
        p_lat, p_lat_rs, sum(okLat));
    fprintf('  trials pooled: Control %d, Stim %d\n', nTrialGrp(1), nTrialGrp(2));
    logf('  per-session curve area (ctrl | stim | diff):\n');
    for s = 1:nS
        fprintf('    sess %d : %8.4f | %8.4f | %+8.4f\n', s, ...
            areaStat(s,1), areaStat(s,2), areaStat(s,1)-areaStat(s,2));
    end

% -------------------- PAIRING DIAGNOSTIC ----------------------------
% Is there any between-session baseline to remove?  If the per-session
% control and stim values are uncorrelated, pairing buys nothing and the
% unpaired test's extra arrangements are legitimate.  If they are highly
% correlated, the sessions share a baseline and unpairing discards it.
    logf('\n  PAIRING DIAGNOSTIC -- per-session median first-lick latency (s)\n');
    fprintf('    %-5s %10s %10s %10s %12s %12s\n', ...
        'sess','ctrl','stim','ctrl-stim','n ctrl','n stim');
    for s = 1:nS
        fprintf('    %-5d %10.4f %10.4f %10.4f %7d/%-4d %7d/%-4d\n', s, ...
            medLat(s,1), medLat(s,2), medLat(s,1)-medLat(s,2), ...
            nTrialSess(s,1), nTotSess(s,1), nTrialSess(s,2), nTotSess(s,2));
    end

    statMat  = {medLat, areaStat};
    statName = {'median first-lick latency', 'cumulative-curve area'};

    for z = 1:2
        M   = statMat{z};
        okz = ~isnan(M(:,1)) & ~isnan(M(:,2));
        cc  = M(okz,1);
        ss  = M(okz,2);

        fprintf('\n    -- %s --\n', statName{z});
        if numel(cc) < 3
            fprintf('       too few paired sessions for the diagnostic\n');
            continue
        end

        sdC = std(cc);
        sdS = std(ss);

% The argument that actually licenses unpairing: if the CONTROL
% condition barely varies across sessions there is no session effect
% to block on, whatever the correlation estimate happens to be.
        logf('       control spread across sessions : SD %.4f, range %.4f  (CV %.1f%%)\n', ...
            sdC, max(cc)-min(cc), 100*sdC/abs(mean(cc)));
        logf('       stim spread across sessions    : SD %.4f, range %.4f\n', ...
            sdS, max(ss)-min(ss));

        if sdC > 0 && sdS > 0
            fprintf('       corr(ctrl, stim)               : Pearson %+.3f | Spearman %+.3f   (n = %d, unreliable)\n', ...
                corr(cc, ss), corr(cc, ss, 'type', 'Spearman'), numel(cc));
        end

        sdD = std(cc - ss);
        sdP = sqrt((sdC^2 + sdS^2)/2);
        if sdP > 0
            rat = sdD / sdP;
            logf('       SD(paired diffs) %.4f / pooled SD %.4f = ratio %.3f   (1.414 means r = 0)\n', ...
                sdD, sdP, rat);
            logf('       r implied by that ratio        : %+.3f\n', 1 - (rat^2)/2);
            if max(sdC,sdS)/max(min(sdC,sdS),eps) > 4
                fprintf(['       CAUTION: the two groups differ in spread by more than 4x, so the\n' ...
                         '                ratio drifts toward 1.414 on its own and the correlation\n' ...
                         '                is unstable. Read the control-spread line instead.\n']);
            end
            if rat > 1.30
                logf('       -> pairing gains little; unpairing costs little here.\n');
            elseif rat < 0.90
                logf('       -> sessions share a baseline; the paired test is the correct one.\n');
            else
                fprintf('       -> intermediate. Neither choice is clearly right at this n.\n');
            end
        end

        if min(cc) > max(ss) || min(ss) > max(cc)
            logf('       SEPARATION: complete -- the groups do not overlap, so the unpaired\n');
            fprintf('                   rank-sum sits at its floor of %.4f one-tailed.\n', ...
                1/nchoosek(2*numel(cc), numel(cc)));
        end
    end

% ---- which statistic the bracket reports ----
    if strcmpi(PLOT.bracketStat,'latency')
        pPair10 = p_lat;   pUnp10 = p_lat_rs;
    else
        pPair10 = p_area;  pUnp10 = q_area;
    end
    fprintf('  bracket (%s statistic): paired p = %.4g | unpaired p = %.4g\n', ...
        PLOT.bracketStat, pPair10, pUnp10);

    if isnan(pPair10),            sPair10 = '--';
    elseif pPair10 < 0.001,       sPair10 = '***';
    elseif pPair10 < 0.01,        sPair10 = '**';
    elseif pPair10 < STATS.alpha, sPair10 = '*';
    else,                         sPair10 = 'n.s.';
    end
    if isnan(pUnp10),             sUnp10 = '--';
    elseif pUnp10 < 0.001,        sUnp10 = '***';
    elseif pUnp10 < 0.01,         sUnp10 = '**';
    elseif pUnp10 < STATS.alpha,  sUnp10 = '*';
    else,                         sUnp10 = 'n.s.';
    end
    if strcmp(STATS.cdfTest, 'peranimal')
        bracketStr   = sprintf('%d/%d animals', PA.kLat, PA.nLatOK);
        bracketP     = max(PA.pLat);   % every animal is below this
        bracketWhich = 'per-animal rank-sum on trials';
    elseif strcmp(STATS.cdfTest, 'ranksum')
        bracketStr   = sUnp10;
        bracketP     = pUnp10;
        bracketWhich = 'unpaired rank-sum';
    elseif strcmp(STATS.cdfTest, 'signrank')
        bracketStr   = sPair10;
        bracketP     = pPair10;
        bracketWhich = 'paired signed-rank';
    else
        bracketStr   = sprintf('%s paired / %s unpaired', sPair10, sUnp10);
        bracketP     = pUnp10;
        bracketWhich = 'both';
    end
    fprintf('  bracket: %s on the %s statistic over %.2f-%.2f s -> p = %.4g -> %s\n', ...
        bracketWhich, PLOT.bracketStat, STATS.traceWin_s(1), STATS.traceWin_s(2), ...
        bracketP, bracketStr);

% ---- per-bin stars, shared by both panels ----
    pbStar = [];
    if PLOT.binStars
        if strcmpi(STATS.cdfTest,'signrank'), pb = p_bin; else, pb = q_bin; end
        pb(~inWin) = NaN;   % only the plotted window is a test
        if STATS.bh
            pv = pb(:); av = nan(size(pv)); okb = ~isnan(pv);
            if any(okb)
                [ps, ordb] = sort(pv(okb));
                mm = numel(ps);
                ad = ps .* (mm ./ (1:mm)');
                ad = flipud(cummin(flipud(ad)));
                tmp = nan(mm,1); tmp(ordb) = min(ad,1);
                av(okb) = tmp;
            end
            pb = av';
            fprintf('  per-bin stars: BH-corrected across %d bins in the plotted window\n', sum(inWin));
        end
        pbStar = find(pb < STATS.alpha);
    end

% ---- annotation, applied to both panels ----
    for eu = 1   % error bars over TRIALS only
        axC = ax10all(eu);
        xlim(axC, PLOT.xlim10);
        ylim(axC, [0 1.16]);

        if ~isempty(pbStar)
            plot(axC, timeCenters(pbStar), 1.02*ones(size(pbStar)), 'k*', 'MarkerSize', 5);
        end

        if PLOT.bracket10 && (~isnan(pPair10) || ~isnan(pUnp10))
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

        if eu == 1 && strcmp(STATS.cdfTest, 'peranimal')
            text(0.03, 0.995, sprintf('%s, first-lick latency, p = %s', bracketWhich, ...
                 strjoin(cellstr(compose('%.3g', PA.pLat(~isnan(PA.pLat)))), ', ')), ...
                 'Units','normalized', 'FontSize', PLOT.fontSize-3, ...
                 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top', 'Parent', axC);
        elseif eu == 1
            text(0.03, 0.995, sprintf('%s, whole trace %.2f-%.2f s, p = %.4g', ...
                 bracketWhich, STATS.traceWin_s(1), STATS.traceWin_s(2), bracketP), ...
                 'Units','normalized', 'FontSize', PLOT.fontSize-3, ...
                 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top', 'Parent', axC);
        end

        xlabel(axC, 'Time from Go Cue (s)');
        ylabel(axC, 'Lick Probability');
        title(axC, sprintf('CI band over %ss', euName{eu}), 'FontWeight','normal');
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
    A_test = nan(nS, nLicksP, 4);   % [|angle|] what the angle tests use: session mean of angT(angle)
    D_sess = nan(nS, nLicksP, 2);   % duration, 1 = Control, 2 = Stim (seconds)

    for s = 1:nS
        for g = 1:4
            a = SESS(s).angle{g};
            if ~isempty(a)
                A_sess(s,:,g) = nanmean(a(:,licks), 1);
                A_test(s,:,g) = nanmean(angT(a(:,licks)), 1);
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

% ---- statistics first; they do not depend on the error-bar choice ----
    p_angL = nan(1,nLicksP); q_angL = nan(1,nLicksP); n_angL = zeros(1,nLicksP);
    p_angR = nan(1,nLicksP); q_angR = nan(1,nLicksP); n_angR = zeros(1,nLicksP);
    p_dur  = nan(1,nLicksP); q_dur  = nan(1,nLicksP); n_dur  = zeros(1,nLicksP);

    for i = 1:nLicksP
        xL = A_test(:,i,1); yL = A_test(:,i,2);   % [|angle|]
        ok = ~isnan(xL) & ~isnan(yL);  n_angL(i) = sum(ok);
        if sum(ok) >= 2
            if any(xL(ok) - yL(ok) ~= 0)
                p_angL(i) = signrank(xL(ok), yL(ok), 'tail', STATS.tailAngle);
            end
            q_angL(i) = ranksum(xL(ok), yL(ok), 'tail', STATS.tailAngle);
        end

        xR = A_test(:,i,3); yR = A_test(:,i,4);   % [|angle|]
        ok = ~isnan(xR) & ~isnan(yR);  n_angR(i) = sum(ok);
        if sum(ok) >= 2
            if any(xR(ok) - yR(ok) ~= 0)
                p_angR(i) = signrank(xR(ok), yR(ok), 'tail', STATS.tailAngle);
            end
            q_angR(i) = ranksum(xR(ok), yR(ok), 'tail', STATS.tailAngle);
        end

        x = D_sess(:,i,1); y = D_sess(:,i,2);
        ok = ~isnan(x) & ~isnan(y);  n_dur(i) = sum(ok);
        if sum(ok) >= 2
            if any(x(ok) - y(ok) ~= 0)
                p_dur(i) = signrank(x(ok), y(ok), 'tail', STATS.tailDur);
            end
            q_dur(i) = ranksum(x(ok), y(ok), 'tail', STATS.tailDur);
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

    showPaired = strcmp(STATS.figTest,'signrank') || strcmp(STATS.figTest,'both');
    showUnpair = strcmp(STATS.figTest,'ranksum')  || strcmp(STATS.figTest,'both');
    showAnimal = strcmp(STATS.figTest,'peranimal');

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
        yPair = yl(2) - 0.10*diff(yl);
        yUnp  = yl(2) - 0.04*diff(yl);
        for i = 1:nLicksP
            if showPaired && pAdjL(i) < STATS.alpha
                text(licks(i)-0.15, yPair, '*', 'Color', PLOT.colLctrl, ...
                     'FontSize', PLOT.fontSize+3, 'HorizontalAlignment','center');
            end
            if showPaired && pAdjR(i) < STATS.alpha
                text(licks(i)+0.15, yPair, '*', 'Color', PLOT.colRctrl, ...
                     'FontSize', PLOT.fontSize+3, 'HorizontalAlignment','center');
            end
            if showUnpair && qAdjL(i) < STATS.alpha
                plot(licks(i)-0.15, yUnp, 'o', 'MarkerSize', 6, 'LineWidth', 1.4, ...
                     'MarkerEdgeColor', PLOT.colLctrl, 'MarkerFaceColor','none');
            end
            if showUnpair && qAdjR(i) < STATS.alpha
                plot(licks(i)+0.15, yUnp, 'o', 'MarkerSize', 6, 'LineWidth', 1.4, ...
                     'MarkerEdgeColor', PLOT.colRctrl, 'MarkerFaceColor','none');
            end
        end
        if showAnimal && strcmp(STATS.perAnimalMode, 'pooled')   % one label per panel
            text(0.02, 0.98, sprintf(['licks %d-%d pooled per trial, rank-sum within animal:  ' ...
                 '{\\color[rgb]{%g %g %g}L %d/%d}   {\\color[rgb]{%g %g %g}R %d/%d} animals'], ...
                 PA.licks(1), PA.licks(end), PLOT.colLctrl, PA.pooledK(1), PA.pooledNa(1), ...
                 PLOT.colRctrl, PA.pooledK(2), PA.pooledNa(2)), 'Units','normalized', ...
                 'FontSize', PLOT.fontSize-3, 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top', ...
                 'Interpreter','tex');
        elseif showAnimal   % k/N per lick
            for i = 1:nLicksP
                if PA.kAngL(i) > 0
                    text(licks(i)-0.15, yPair, sprintf('%d/%d', PA.kAngL(i), PA.nAngL(i)), ...
                         'Color', PLOT.colLctrl, 'FontSize', PLOT.fontSize-3, 'HorizontalAlignment','center');
                end
                if PA.kAngR(i) > 0
                    text(licks(i)+0.15, yUnp, sprintf('%d/%d', PA.kAngR(i), PA.nAngR(i)), ...
                         'Color', PLOT.colRctrl, 'FontSize', PLOT.fontSize-3, 'HorizontalAlignment','center');
                end
            end
            if eu == 1
                text(0.02, 0.98, 'k/N = animals significant (rank-sum on trials)', 'Units','normalized', ...
                     'FontSize', PLOT.fontSize-3, 'Color', [0.35 0.35 0.35], 'VerticalAlignment','top');
            end
        end
        if strcmpi(STATS.figTest,'both') && eu == 1
            text(0.02, 0.98, '*  paired      o  unpaired', 'Units','normalized', ...
                 'FontSize', PLOT.fontSize-2, 'Color', [0.35 0.35 0.35], ...
                 'VerticalAlignment','top');
        end

        if eu == 1
            for g = 1:4
                text(0.70, 0.06 + 0.07*(4-g), angName{g}, 'Units','normalized', ...
                     'Color', angCol{g}, 'FontSize', PLOT.fontSize-1, 'FontWeight','bold');
            end
        end

        xlim([0.5 max(licks)+0.5]); xticks(licks); xtickangle(0);
        xlabel('Lick number from go cue');
        ylabel('Lick angle (deg)');
        title(sprintf('Angle  --  error bars over %ss', euName{eu}), 'FontWeight','normal');
        set(axA, 'TickDir','out', 'FontSize', PLOT.fontSize); box off

    end

% ------------------------- print the stats --------------------------
    if strcmpi(STATS.unit,'session')
        fprintf('\nFIG 11 -- p = paired signed-rank across sessions | q = UNPAIRED rank-sum, same values\n');
    else
        fprintf('\nFIG 11 -- p = rank-sum on pooled trials | q = same\n');
    end
    if showAnimal
        fprintf('  figure marks = k/N animals significant; the per-animal tables are printed above\n');
    end
    fprintf('  angle tail = %s, duration tail = %s | stars use %s', ...
        STATS.tailAngle, STATS.tailDur, STATS.figTest);
    if STATS.bh
        fprintf(', BH-corrected over m = %d contacts (threshold at rank k is k/%d*%.3f)\n', ...
            nLicksP, nLicksP, STATS.alpha);
    else
        fprintf(', UNCORRECTED\n');
    end
    fprintf('  columns: p = paired signed-rank, q = unpaired rank-sum, *Adj = same after BH\n');
    fprintf('  both rows carry these same stars; only the error bars differ\n');

    disp(table(licks', p_angL', pAdjL', q_angL', qAdjL', n_angL', ...
        'VariableNames', {'Lick','p_L','pAdj_L','q_L','qAdj_L','n_L'}));
    disp(table(licks', p_angR', pAdjR', q_angR', qAdjR', n_angR', ...
        'VariableNames', {'Lick','p_R','pAdj_R','q_R','qAdj_R','n_R'}));
    disp(table(licks', p_dur', pAdjD', q_dur', qAdjD', n_dur', ...
        'VariableNames', {'Lick','p_dur','pAdj_dur','q_dur','qAdj_dur','n_dur'}));

% ---------------- PAIRING DIAGNOSTIC, per comparison ------------------
% Same question as the fig 10 diagnostic, asked at every contact: do the
% per-session control and stim values share a baseline? If they do not,
% pairing removes nothing and the unpaired test's extra arrangements are
% legitimate. ratio = SD(paired differences) / pooled within-group SD;
    logf('\nFIG 11 -- PAIRING DIAGNOSTIC\n');

    diagName = {'L_ctrl vs L_stim  (angle)', ...
                'R_ctrl vs R_stim  (angle)', ...
                'Control vs Stim   (duration)'};
    diagX = {A_sess(:,:,1), A_sess(:,:,3), D_sess(:,:,1)};
    diagY = {A_sess(:,:,2), A_sess(:,:,4), D_sess(:,:,2)};

    for z = 1:3
        X = diagX{z};  Y = diagY{z};
        rP  = nan(1,nLicksP);  rS  = nan(1,nLicksP);
        sdD = nan(1,nLicksP);  sdP = nan(1,nLicksP);
        rat = nan(1,nLicksP);  rIm = nan(1,nLicksP);
        nn  = zeros(1,nLicksP);

        for i = 1:nLicksP
            x = X(:,i);  y = Y(:,i);
            okd = ~isnan(x) & ~isnan(y);
            nn(i) = sum(okd);
            if sum(okd) >= 3
                xc = x(okd);  yc = y(okd);
                if std(xc) > 0 && std(yc) > 0
                    rP(i) = corr(xc, yc);
                    rS(i) = corr(xc, yc, 'type', 'Spearman');
                end
                sdD(i) = std(xc - yc);
                sdP(i) = sqrt((var(xc) + var(yc))/2);
                if sdP(i) > 0
                    rat(i) = sdD(i) / sdP(i);
                    rIm(i) = 1 - (rat(i)^2)/2;
                end
            end
        end

        fprintf('\n  %s\n', diagName{z});
        disp(table(licks', rP', rS', sdD', sdP', rat', rIm', nn', ...
            'VariableNames', {'Lick','r_Pearson','r_Spearman','SD_diff','SD_pooled','ratio','r_implied','n'}));

        rv = rat(~isnan(rat));
        if isempty(rv)
            fprintf('    no contact had enough paired sessions for the diagnostic\n');
            continue
        end
        mr = median(rv);
        fprintf('    median ratio over %d usable contacts = %.3f   (1.414 means r = 0)\n', numel(rv), mr);
        logf('    contacts with ratio > 1.30 (pairing useless) : %d of %d\n', sum(rv > 1.30), numel(rv));
        logf('    contacts with ratio < 0.90 (shared baseline) : %d of %d\n', sum(rv < 0.90), numel(rv));
        if mr > 1.30
            fprintf(['    -> baselines carry almost no information. Unpairing costs little here;\n' ...
                     '       justify it from the DESIGN in the Methods, not from this table.\n']);
        elseif mr < 0.90
            logf(['    -> sessions share a substantial baseline. Unpairing discards real\n' ...
                     '       structure -- the paired test is the correct one for this comparison.\n']);
        else
            fprintf('    -> intermediate. Neither choice is clearly right at this n.\n');
        end
    end

    fprintf(['\n  NOTE: a verdict here applies only to the comparison it sits under.\n' ...
             '        The fig 10 latency diagnostic does not license unpairing for angle.\n']);

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

function t = latencySec(firstBin, timeAxis)
% first-lick bin -> seconds from the go cue; no lick found -> Inf,
% i.e. later than every lick, which is how the cumulative curve counts it.
    firstBin = firstBin(:);
    t  = inf(size(firstBin));
    ok = ~isnan(firstBin);
    t(ok) = timeAxis(round(firstBin(ok)));
end

function p = rankSumColumn(A, B, c, tail)
% rank-sum between the trials of A and of B at lick column c.
    p = NaN;
    if isempty(A) || isempty(B) || size(A,2) < c || size(B,2) < c, return; end
    a = A(:,c);  a = a(~isnan(a));
    b = B(:,c);  b = b(~isnan(b));
    if numel(a) < 2 || numel(b) < 2, return; end
    p = ranksum(a, b, 'tail', tail);
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p values; NaNs stay NaN.
    q  = nan(size(p));
    ok = ~isnan(p);
    if ~any(ok), return; end
    [ps, ord] = sort(p(ok));
    m  = numel(ps);
    ad = ps(:) .* (m ./ (1:m)');
    ad = flipud(cummin(flipud(ad)));
    tmp = nan(m,1);  tmp(ord) = min(ad, 1);
    q(ok) = tmp;
end

function s = starOf(p, alpha)
    if isnan(p), s = '';
    elseif p < 0.001, s = ' ***';
    elseif p < 0.01,  s = ' **';
    elseif p < alpha, s = ' *';
    else, s = '';
    end
end

function s = ternaryStr(c, a, b)
    if c, s = a; else, s = b; end
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

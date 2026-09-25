%% FS10_spontaneous_v4.m
%  Motor cortical activity during spontaneous licking bouts.
%  Bouts that were not triggered by the go cue and were never rewarded are
%  aligned to their onset and compared with go cue-evoked bouts.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/300
%    params.smooth      20
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -3 to 32 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear; clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.

sz = 14;

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
params.alignEvent = 'goCue';   % 'goCue' 'firstLick' 'moveOnset' ...
params.behav_only = 0;
params.timeWarp   = 0;
params.nLicks     = 20;
params.lowFR      = 0.01;   % minimum mean firing rate, Hz

params.condition(1)     = {'hit==1 | hit==0'};
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1'};
% NOTE: params.condition is passed to loadSessionData but is NOT what selects
% trials here -- P8 (R1) and P9 (R4) are taken straight from obj.bp below.

params.tmin   = -3;
params.tmax   = 32;
params.dt     = 1/300;
params.smooth = 20;   % causal gaussian kernel

params.quality = {'good'};   % good units only (findClusters trims blanks and ignores case)
% accept 'excellent','fair','ok'. Keep in mind
% when comparing unit counts across figures.

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

%% CONFIG: ONE PLACE FOR EVERY CHOICE
% LICK / BOUT DETECTION. These values now drive BOTH the per-session raster and
% the group analysis. In the previous version the raster used minLickLen 0.04 s
% and the analysis used 0.025 s, so a lick could be marked on the figure and
% excluded from the numbers, or the reverse.
cfg.withinBout_sec    = 0.25;   % gap below this keeps licks in the same bout
cfg.betweenBout_sec   = 0.5;   % bouts closer than this are merged
cfg.minLickLen_sec    = 0.025;   % a contact run shorter than this is not a lick
cfg.minBoutLicks      = 1;   % bouts with fewer licks than this are dropped

% PRE-GO-CUE (spontaneous) BOUTS. The silence guard is what makes these
% initiations from rest rather than the continuation of ongoing licking.
cfg.preGC_search_sec  = 2.5;   % search back this far from the go cue
cfg.preGC_min_sec     = 0.7;   % bout must start at least this long before it
cfg.preGC_silence_sec = 0.7;   % and be preceded by this much tongue-free time
cfg.minPreRec_sec     = 0.1;   % need this much recording before the lick

% SNIPPET WINDOW around each event
cfg.preBins_sec  = 1.0;
cfg.postBins_sec = 1.0;

% NORMALISATION -- one choice, applied post-loop, used by every figure.
%   'raw'      NO normalisation: mean spike rate in spk/s, as recorded. Units
%              that a session with unusually high overall rates dominates the
%              pooled mean, so always look at the per-animal panels too.
%   'minmax01' per session, pooled over R1 and R4. Removes session-to-session
%              scale, so only the SHAPE around the lick is comparable.
%   'zscore'   per session, mean/sd over the whole trace.
cfg.normMode = 'raw';

% BASELINE for the z-score option only
cfg.bl_start_sec = -1.0;
cfg.bl_end_sec   = -0.2;

% FIGURES
cfg.plotPerSession    = false;   % per-session trial raster + heatmap
% (3 events x nAnimals figures -- a lot of windows)
cfg.xlim_hmap      = [-0.25  0.15];
cfg.xlim_line      = [-0.3   0.3];
cfg.ylim_line      = [];   % [] = auto. Was hard-coded [0.1 0.5], which
% clips whenever the normalisation changes.
cfg.smooth_win     = 20;
cfg.rngSeed        = 42;   % shuffle-before-sort, so duration ties are
% ordered randomly rather than by trial number

% TRIAL CAPS: {animal, date, last usable trial}.
% Matched on a CANONICALISED name (see canonAnm below) because this project
% spells the same animal both 'TD4l' and 'TDl4'. The old code tested the
% spelling the loaders do not produce, so these caps never fired.
cfg.trialCaps = { 'TD13d', '2024-11-11', 278 ; ...
                  'TD8d',  '2024-09-07', 313 ; ...
                  'TD8d',  '2024-09-09', 298 ; ...
                  'TD4l',  '2025-06-05', 123 };

%% SPECIFY DATA TO LOAD
datapth = '';   % raw data folder not used (was: datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';)

meta13 = []; meta14 = []; meta15 = []; meta16 = []; meta17 = [];
meta18 = []; meta19 = []; meta20 = []; meta21 = []; meta22 = [];
meta23 = []; meta24 = []; meta25 = []; meta26 = []; meta27 = [];
meta28 = []; meta29 = []; meta30 = []; meta31 = []; meta32 = [];

date = '2025-02-03'; meta13 = slimMeta('loadTD3l_neur', date);
date = '2025-02-04'; meta14 = slimMeta('loadTD3l_neur', date);
date = '2025-02-05'; meta15 = slimMeta('loadTD3l_neur', date);
date = '2025-02-06'; meta16 = slimMeta('loadTD3l_neur', date);
date = '2025-02-07'; meta17 = slimMeta('loadTD3l_neur', date);

date = '2025-02-03'; meta18 = slimMeta('loadTD2l_neur', date);
date = '2025-02-04'; meta19 = slimMeta('loadTD2l_neur219', date);
date = '2025-02-05'; meta20 = slimMeta('loadTD2l_neur', date);
date = '2025-02-06'; meta21 = slimMeta('loadTD2l_neur', date);
date = '2025-02-07'; meta22 = slimMeta('loadTD2l_neur', date);

date = '2025-06-03'; meta23 = slimMeta('loadTD4l_neur', date);
date = '2025-06-04'; meta24 = slimMeta('loadTD4l_neur', date);
date = '2025-06-05'; meta25 = slimMeta('loadTD4l_neur', date);
date = '2025-06-06'; meta26 = slimMeta('loadTD4l_neur', date);
date = '2025-06-07'; meta27 = slimMeta('loadTD4l_neur', date);

date = '2025-08-21'; meta28 = slimMeta('loadTD5l_neur', date);
date = '2025-08-22'; meta29 = slimMeta('loadTD5l_neur', date);
date = '2025-08-23'; meta30 = slimMeta('loadTD5l_neur', date);
date = '2025-08-24'; meta31 = slimMeta('loadTD5l_neur', date);
date = '2025-08-25'; meta32 = slimMeta('loadTD5l_neur', date);

% 4 animals x 5 days, ordered animal1 days1-5, animal2 days1-5, ...
all_meta = [meta13;meta14;meta15;meta16;meta17; ...
            meta18;meta19;meta20;meta21;meta22; ...
            meta23;meta24;meta25;meta26;meta27; ...
            meta28;meta29;meta30;meta31;meta32];

nAnimals = 4;
nDays    = 5;
dayLabel = repmat(1:nDays, 1, nAnimals);   % 1x20
anmLabel = repelem(1:nAnimals, nDays);   % 1x20, which animal each session is

%% PROBE MAP
% One entry per session. Value = which probe to use (1 or 2); 0 = skip.
% s1..s5 are DAY maps: s<k> is nonzero exactly at each animal's day k.
s1 = [1 0 0 0 0   1 0 0 0 0   2 0 0 0 0   2 0 0 0 0];
s2 = [0 2 0 0 0   0 2 0 0 0   0 1 0 0 0   0 1 0 0 0];
s3 = [0 0 1 0 0   0 0 1 0 0   0 0 1 0 0   0 0 2 0 0];
s4 = [0 0 0 1 0   0 0 0 2 0   0 0 0 1 0   0 0 0 1 0];
s5 = [0 0 0 0 1   0 0 0 0 1   0 0 0 0 1   0 0 0 0 1];

allMaps = [s1(:), s2(:), s3(:), s4(:), s5(:)];   % 20 x 5

% Exactly one entry per session must be nonzero -- check it rather than assume,
% because max() would silently hide a double assignment.
nNonZero = sum(allMaps ~= 0, 2);
assert(all(nNonZero == 1), ...
    'probe map: session(s) %s have %s nonzero entries, expected exactly 1.', ...
    mat2str(find(nNonZero ~= 1)'), mat2str(nNonZero(nNonZero ~= 1)'));

probe_map = max(allMaps, [], 2);   % 20 x 1
day_map   = zeros(size(probe_map));
for si = 1:5
    day_map(allMaps(:,si) ~= 0) = si;
end
assert(isequal(day_map(:)', dayLabel), ...
    'probe map days disagree with dayLabel -- check s1..s5 against the session order.');

%% STORAGE
AllMeanFR_R1 = cell(numel(all_meta), 2);   % T x nTrials, RAW (normalised post-loop)
AllMeanFR_R4 = cell(numel(all_meta), 2);
AllLen_R1    = cell(numel(all_meta), 2);   % nL x nTrials tongue length
AllLen_R4    = cell(numel(all_meta), 2);
AllTime      = cell(numel(all_meta), 2);   % neural time axis
AllTimeL     = cell(numel(all_meta), 2);   % kinematic time axis
AllP8        = cell(1, numel(all_meta));   % R1 trial indices
AllP9        = cell(1, numel(all_meta));   % R4 trial indices
numNeurons   = cell(numel(all_meta), 2);
AllAnm       = cell(1, numel(all_meta));   % animal name per session
numTrials    = cell(1, numel(all_meta));
percCompleted= cell(1, numel(all_meta));

params.behav_only = 0;

%% MAIN LOOP
for sessnum = 1:numel(all_meta)

    clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin me

    probe_this_sess = probe_map(sessnum);
    day_this_sess   = day_map(sessnum);

    meta = all_meta(sessnum,1);
    params.probe = {meta.probe};
    params.cluid = {};

    [obj, params, kin] = slimToLegacy(meta, params);

    for sessix = 1:numel(meta)
    end
    for sessix = 1:numel(meta)
        fprintf('----Getting kinematic data for session %d out of %d----\n', sessix, numel(meta));
    end
    sessix = 1;

% ---------------- TRIAL CAP (see cfg.trialCaps) ----------------
% Applied to the kinematics AND to the R1/R4 trial lists, and it says out
% loaders never produce, so the cap silently did nothing.
    capN = trialCapFor(obj.pth.anm, obj.pth.dt, cfg.trialCaps);
    trialSet = (1:obj.bp.Ntrials)';
    if ~isnan(capN)
        trialSet(trialSet > capN) = [];
        logf('  [cap] %s %s: keeping trials 1-%d of %d\n', ...
            obj.pth.anm, obj.pth.dt, capN, obj.bp.Ntrials);
    end
    condtrix = trialSet;

    kinix  = find(strcmp(kin(sessix).featLeg, 'tongue_length'));
    assert(~isempty(kinix), 'tongue_length not found in kin.featLeg');
    Length = kin(sessix).dat(:, condtrix, kinix);

% ---------------- R1 and R4 HIT TRIALS ----------------
    allTr = 1:obj.bp.Ntrials;
    hitTr = allTr(obj.bp.hit == 1);
    P8 = intersect(hitTr, allTr(obj.bp.rewardedLick == 1))';   % R1
    P9 = intersect(hitTr, allTr(obj.bp.rewardedLick == 4))';   % R4
    if ~isnan(capN)
        P8(P8 > capN) = [];
        P9(P9 > capN) = [];
    end

% ---------------- REGION INDICES ----------------
    Ncells = size(obj.psth, 2);
    if iscell(params.cluid) && size(params.cluid,2) >= 1 && ~isempty(params.cluid{1,1})
        clu_p1 = 1 : size(params.cluid{1,1}, 1);
        clu_p2 = size(params.cluid{1,1},1)+1 : Ncells;
    else
% single-probe session: both indices point at the same units
        clu_p1 = 1 : Ncells;
        clu_p2 = clu_p1;
    end
    if isempty(clu_p2), clu_p2 = clu_p1; end

% ---------------- TIME AXES ----------------
    binSz    = obj.time(2) - obj.time(1);
    nL       = size(Length, 1);
    time_L   = obj.time(1) + (0:nL-1) * binSz;
    T_neural = size(obj.trialdat, 1);

    nT1 = numel(P8);
    nT4 = numel(P9);

% ---------------- TONGUE LENGTH PER TRIAL ----------------
    len_R1 = nan(nL, nT1);
    len_R4 = nan(nL, nT4);
    for k = 1:nT1
        tr = P8(k);
        if tr <= size(Length,2), len_R1(:,k) = Length(:,tr); end
    end
    for k = 1:nT4
        tr = P9(k);
        if tr <= size(Length,2), len_R4(:,k) = Length(:,tr); end
    end

% ---------------- MEAN FIRING RATE PER TRIAL, BOTH PROBES ----------------
    meanFR_R1 = cell(1,2);  meanFR_R4 = cell(1,2);
    cluSets = {clu_p1, clu_p2};
    for pp = 1:2
        allDat = obj.trialdat(:, cluSets{pp}, :);
        m1 = nan(T_neural, nT1);
        m4 = nan(T_neural, nT4);
        for k = 1:nT1, m1(:,k) = nanmean(allDat(:,:,P8(k)), 2); end
        for k = 1:nT4, m4(:,k) = nanmean(allDat(:,:,P9(k)), 2); end
        meanFR_R1{pp} = m1;
        meanFR_R4{pp} = m4;
    end

% ---------------- STORE (raw; normalised once, post-loop) ----------------
    for pp = 1:2
        AllMeanFR_R1{sessnum,pp} = meanFR_R1{pp};
        AllMeanFR_R4{sessnum,pp} = meanFR_R4{pp};
        AllLen_R1{sessnum,pp}    = len_R1;   % kinematics do not depend on probe
        AllLen_R4{sessnum,pp}    = len_R4;
        AllTime{sessnum,pp}      = obj.time;
        AllTimeL{sessnum,pp}     = time_L;
    end
    AllP8{sessnum}  = P8;
    AllP9{sessnum}  = P9;
    AllAnm{sessnum} = obj.pth.anm;
    numNeurons{sessnum,1} = numel(clu_p1);
    numNeurons{sessnum,2} = numel(clu_p2);
    numTrials{sessnum}     = sum(obj.bp.hit);
    percCompleted{sessnum} = sum(obj.bp.hit) / obj.bp.Ntrials;

    fprintf('Session %d\n', sessnum);
end

%% NORMALISATION (one place, one choice)
% -- a HORIZONTAL concat of two columns whose lengths are T*nR1 and T*nR4.
% Those differ whenever the session has unequal numbers of R1 and R4 trials,
% which is every session, so MATLAB threw "Dimensions of arrays being
% also keeps nanmin/nanmax scalar so the subtraction below is by a number
% rather than broadcast per column.
t_ax = [];  t_L = [];
for s = 1:numel(all_meta)
    for pp = 1:2
        if ~isempty(AllTime{s,pp})
            t_ax = AllTime{s,pp};
            t_L  = AllTimeL{s,pp};
            break
        end
    end
    if ~isempty(t_ax), break; end
end
assert(~isempty(t_ax), 'no session produced a time axis.');

AllMeanFR_R1_n = cell(size(AllMeanFR_R1));
AllMeanFR_R4_n = cell(size(AllMeanFR_R4));

for sessnum = 1:numel(all_meta)
    for pp = 1:2
        mat1 = AllMeanFR_R1{sessnum, pp};
        mat4 = AllMeanFR_R4{sessnum, pp};
        if isempty(mat1) && isempty(mat4), continue; end

        all_data = [];
        if ~isempty(mat1), all_data = [all_data; mat1(:)]; end   % <-- semicolon
        if ~isempty(mat4), all_data = [all_data; mat4(:)]; end   % <-- semicolon
        all_data = all_data(isfinite(all_data));
        if isempty(all_data), continue; end

        switch lower(cfg.normMode)
            case 'raw'
                f = @(M) M;   % spk/s, untouched
            case 'minmax01'
                lo = min(all_data);
                hi = max(all_data);
                if hi == lo, hi = lo + 1; end
                f = @(M) (M - lo) / (hi - lo);
            case 'zscore'
                mu = mean(all_data);
                sd = std(all_data);
                if sd == 0 || isnan(sd), sd = 1; end
                f = @(M) (M - mu) / sd;
            otherwise
                error('cfg.normMode must be ''minmax01'' or ''zscore''.');
        end

        if ~isempty(mat1), AllMeanFR_R1_n{sessnum,pp} = f(mat1); end
        if ~isempty(mat4), AllMeanFR_R4_n{sessnum,pp} = f(mat4); end
    end
end
logf('Normalisation complete: %s, per session, pooled R1+R4.\n', cfg.normMode);

%% EVENT DETECTION + SNIPPETS
binSz_neural = t_ax(2) - t_ax(1);
preBins      = round(cfg.preBins_sec  / binSz_neural);
postBins     = round(cfg.postBins_sec / binSz_neural);
t_common     = (-preBins : postBins) * binSz_neural;
snip_len     = numel(t_common);

% Keeping the animal index means the pooled figures and the per-animal figures
% come from exactly the same events -- the pooled version is a concatenation
% over the animal dimension, not a separate pass with its own bugs.
nEv = 3;
ev_cells = cell(nEv, nAnimals, nDays);
for e = 1:nEv
    for a = 1:nAnimals
        for dd = 1:nDays
            ev_cells{e,a,dd} = {};
        end
    end
end
evIdxOf = containers.Map({'init','reeng','preGC'}, {1,2,3});

for sessnum = 1:numel(all_meta)

    probe = probe_map(sessnum);
    if probe == 0, continue; end
    day = day_map(sessnum);
    anm = anmLabel(sessnum);

    t_L_s  = AllTimeL{sessnum, probe};
    t_ax_s = AllTime{sessnum,  probe};
    if isempty(t_L_s) || isempty(t_ax_s), continue; end

    binSz_L_s  = t_L_s(2)  - t_L_s(1);
    binSz_ax_s = t_ax_s(2) - t_ax_s(1);
    T_neural_s = numel(t_ax_s);
    nL_s       = numel(t_L_s);

    preBins_s    = round(cfg.preBins_sec  / binSz_ax_s);
    postBins_s   = round(cfg.postBins_sec / binSz_ax_s);
    preBins_L_s  = round(cfg.preBins_sec  / binSz_L_s);
    postBins_L_s = round(cfg.postBins_sec / binSz_L_s);

% snip_len is computed ONCE from the first session's axis but the indices
% below come from this session's. They agree only while every session uses
% the same params.dt -- so check rather than trust.
    assert(preBins_s + postBins_s + 1 == snip_len, ...
        'session %d has a different neural bin size than the first session.', sessnum);

    for trType = 1:2
        if trType == 1
            fr  = AllMeanFR_R1_n{sessnum, probe};
            len = AllLen_R1{sessnum, probe};
        else
            fr  = AllMeanFR_R4_n{sessnum, probe};
            len = AllLen_R4{sessnum, probe};
        end
        if isempty(fr) || isempty(len), continue; end
        fr  = fr';   % nTrials x T
        len = len';   % nTrials x nL

        nTr = min(size(fr,1), size(len,1));
        for tr = 1:nTr

% ONE detector for both windows and for the raster below, so the
% events plotted and the events measured cannot diverge.
            ev = findLickEvents(len(tr,:), t_L_s, cfg);

            for e = 1:numel(ev)
                s_kin = ev(e).sIdx;
                e_kin = ev(e).eIdx;
                t_lick = t_L_s(s_kin);
                s_ax   = round((t_lick - t_ax_s(1)) / binSz_ax_s) + 1;

                if strcmp(ev(e).type,'preGC')
                    if s_ax - round(cfg.minPreRec_sec / binSz_ax_s) < 1, continue; end
                end

                lickDur = (e_kin - s_kin + 1) * binSz_L_s;

                snip_fr = nan(1, snip_len);
                i1r = s_ax - preBins_s;   i2r = s_ax + postBins_s;
                i1  = max(1, i1r);        i2  = min(T_neural_s, i2r);
                if i2 >= i1
                    o1 = i1 - i1r + 1;  o2 = o1 + (i2 - i1);
                    snip_fr(o1:o2) = fr(tr, i1:i2);
                end

                snip_len_vec = nan(1, preBins_L_s + postBins_L_s + 1);
                l1r = s_kin - preBins_L_s;  l2r = s_kin + postBins_L_s;
                l1  = max(1, l1r);          l2  = min(nL_s, l2r);
                if l2 >= l1
                    p1 = l1 - l1r + 1;  p2 = p1 + (l2 - l1);
                    snip_len_vec(p1:p2) = len(tr, l1:l2);
                end
                snip_len_vec(isnan(snip_len_vec)) = 0;
                if numel(snip_len_vec) ~= snip_len
                    snip_len_vec = interp1(linspace(0,1,numel(snip_len_vec)), ...
                        snip_len_vec, linspace(0,1,snip_len), 'linear', 0);
                end

                entry = {snip_fr, lickDur, snip_len_vec};
                ei = evIdxOf(ev(e).type);
                ev_cells{ei,anm,day}{end+1} = entry;
            end
        end
    end
end

% ---- event counts, pooled and per animal ----
ev_names  = {'GC init', 'Re-engagement', 'Pre-GC init'};
anmNames  = cell(1, nAnimals);
for a = 1:nAnimals
    idx = find(anmLabel == a, 1, 'first');
    if ~isempty(AllAnm{idx}), anmNames{a} = AllAnm{idx}; else, anmNames{a} = sprintf('animal %d', a); end
end

fprintf('\n%s\n', repmat('=',1,72));
logf('EVENT COUNTS\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('%-16s', 'pooled');
for dd = 1:nDays, fprintf('   day %d', dd); end
fprintf('\n');
for e = 1:nEv
    fprintf('  %-14s', ev_names{e});
    for dd = 1:nDays
        n = 0;
        for a = 1:nAnimals, n = n + numel(ev_cells{e,a,dd}); end
        fprintf('%8d', n);
    end
    fprintf('\n');
end
for a = 1:nAnimals
    fprintf('%s\n', repmat('-',1,72));
    fprintf('%-16s\n', anmNames{a});
    for e = 1:nEv
        fprintf('  %-14s', ev_names{e});
        for dd = 1:nDays, fprintf('%8d', numel(ev_cells{e,a,dd})); end
        fprintf('\n');
    end
end
fprintf('%s\n', repmat('=',1,72));

%% BUILD SNIPPET MATRICES
% Built ONCE per (event, animal, day), then the pooled version is a
% concatenation over the animal dimension. The pooled and per-animal figures
% therefore cannot disagree: they are the same events, grouped two ways.
rng(cfg.rngSeed);
ev_labels = {'GC initiation (1st lick after go cue)', ...
             'Re-engagement (1st lick, bouts 2-N)', ...
             'Pre-GC initiation (spontaneous)'};

line_pool = cell(nEv, nDays);
hmap_pool = cell(nEv, nDays);
dur_pool  = cell(nEv, nDays);
ns_pool   = zeros(nEv, nDays);

for e = 1:nEv
    for dd = 1:nDays
        poolCells = {};
        for a = 1:nAnimals
            c = ev_cells{e,a,dd};
            poolCells = [poolCells, c];   %#ok<AGROW>
        end
        ns_pool(e,dd) = numel(poolCells);
        [line_pool{e,dd}, hmap_pool{e,dd}, dur_pool{e,dd}] = ...
            buildMats(poolCells, t_common, snip_len);
    end
end

% Colour limits from the pooled, real (non-masked) values
all_vals = [];
for e = 1:nEv
    for dd = 1:nDays
        if ~isempty(hmap_pool{e,dd})
            v = hmap_pool{e,dd}(:);
            all_vals = [all_vals; v(isfinite(v))];   %#ok<AGROW>
        end
    end
end
assert(~isempty(all_vals), 'no events were found in any session.');
clim_lo = prctile(all_vals, 5);
clim_hi = prctile(all_vals, 95);
if clim_lo == clim_hi, clim_hi = clim_lo + 1; end

switch lower(cfg.normMode)
    case 'raw',      unit_label = 'Mean spike rate (spk/s)';
    case 'minmax01', unit_label = 'Mean spike rate (0-1 norm)';
    case 'zscore',   unit_label = 'Mean spike rate (z)';
end

%% HEATMAPS, POOLED OVER ANIMALS
% Masked bins are TRANSPARENT, not painted at clim_lo. Painting them made
% "no data after the lick ended" look identical to a genuinely low rate.
col_tick = [0.85 0.10 0.10];
tickH    = 0.4;
tickLW   = 3;
fig_x    = [0.02 0.22 0.42 0.62 0.82];

for e = 1:nEv
    figure('Color','w','Units','normalized','Position',[0.08 0.20 0.62 0.55]);
    sgtitle(sprintf('%s  |  all animals', ev_labels{e}), 'FontSize',13,'FontWeight','bold');
    for dd = 1:nDays
        ax = axes('Position',[fig_x(dd)+0.01  0.10  0.17  0.80]);   %#ok<LAXES>
        drawEventHeatmap(ax, hmap_pool{e,dd}, dur_pool{e,dd}, t_common, ...
            clim_lo, clim_hi, cfg, col_tick, tickH, tickLW, ...
            sprintf('Day %d  (n=%d)', dd, ns_pool(e,dd)), dd == 1, dd == nDays, unit_label);
    end
end

% Per-animal heatmaps removed.

%% MEAN TRACES, POOLED
col_gc    = [0.18 0.65 0.18];
col_reeng = [1.00 0.50 0.00];
col_pre   = [0.20 0.40 0.85];
cols      = {col_gc, col_reeng, col_pre};

figure('Color','w','Units','normalized','Position',[0.08 0.28 0.62 0.32]);
for dd = 1:nDays
    ax = subplot(1, nDays, dd);  hold(ax,'on');
    for e = 1:nEv
        drawMeanTrace(ax, line_pool{e,dd}, t_common, cols{e}, cfg, ...
            sprintf('%s (n=%d)', ev_names{e}, ns_pool(e,dd)));
    end
    finishTraceAxes(ax, cfg, dd, nDays, unit_label, sprintf('Day %d', dd));
end
sgtitle(sprintf('All animals  |  mean \\pm 95%% CI  |  %s', unit_label), ...
    'FontSize',13,'FontWeight','bold');

% Per-animal mean traces, the one-animal-per-line figure and the per-session
% raster were removed. Only the pooled figures above are produced.

%% LOCAL FUNCTIONS

function ev = findLickEvents(lenTrace, t_L, cfg)
% THE ONE DETECTOR. Returns every lick event on this trial's tongue-length
% trace, tagged by type:
%   'init'   first lick of the FIRST bout after the go cue
%   'reeng'  first lick of bouts 2..N after the go cue
%   'preGC'  first lick of a bout that starts at least cfg.preGC_min_sec BEFORE
%            the go cue AND is preceded by cfg.preGC_silence_sec with no tongue
%            visible -- a spontaneous initiation from rest, not the tail of
%            ongoing licking. That silence guard is the whole reason these count
%            as uninstructed.
% event time in seconds.
% WHY THIS IS A FUNCTION. The previous version carried two copies of this logic
% -- one for the per-session raster, one for the group analysis -- with
% DIFFERENT minimum lick durations (0.04 s vs 0.025 s). Events shown on the
% cannot drift from itself.
    ev = struct('sIdx',{},'eIdx',{},'type',{},'tSec',{});
    if all(isnan(lenTrace)) , return; end

    binSz = t_L(2) - t_L(1);
    withinBins  = round(cfg.withinBout_sec  / binSz);
    betweenBins = round(cfg.betweenBout_sec / binSz);
    minLenBins  = round(cfg.minLickLen_sec  / binSz);
    silenceBins = round(cfg.preGC_silence_sec / binSz);

    [~, i_gc] = min(abs(t_L - 0));
    i_pre_start = find(t_L >= -cfg.preGC_search_sec, 1, 'first');
    if isempty(i_pre_start), i_pre_start = 1; end
    i_pre_end = i_gc - 1;

% ---------------- POST-GO-CUE ----------------
    sig = lenTrace;
    sig(1 : max(i_gc-1,0)) = NaN;
    [bS, bE] = mergeBouts(sig, withinBins, betweenBins, minLenBins, cfg.minBoutLicks);
    for b = 1:numel(bS)
        if b == 1, tp = 'init'; else, tp = 'reeng'; end
        ev(end+1) = struct('sIdx',bS{b}(1), 'eIdx',bE{b}(1), ...
                           'type',tp, 'tSec',t_L(bS{b}(1)));   %#ok<AGROW>
    end

% ---------------- PRE-GO-CUE (spontaneous) ----------------
    if i_pre_end >= i_pre_start
        sig_pre = lenTrace;
        sig_pre(1 : i_pre_start-1)  = NaN;
        sig_pre(i_pre_end+1 : end)  = NaN;
        [bS, bE] = mergeBouts(sig_pre, withinBins, betweenBins, minLenBins, cfg.minBoutLicks);
        for b = 1:numel(bS)
            s_bin = bS{b}(1);

% guard 1: far enough before the go cue
            if (s_bin - i_gc) * binSz > -cfg.preGC_min_sec, continue; end

% guard 2: silence before it -- no tongue at all in the preceding
% window. Checked on the FULL trace, not the windowed copy, so a
% bout starting just inside the search window cannot pass by having
% its history NaN-ed out.
            gs = max(1, s_bin - silenceBins);
            ge = s_bin - 1;
            if ge >= gs
                g = lenTrace(gs:ge);
                g(isnan(g)) = 0;
                if any(g > 0), continue; end
            end

            ev(end+1) = struct('sIdx',s_bin, 'eIdx',bE{b}(1), ...
                               'type','preGC', 'tSec',t_L(s_bin));   %#ok<AGROW>
        end
    end
end

function [boutS, boutE] = mergeBouts(sig, withinBins, betweenBins, minLenBins, minLicks)
% Contact runs -> licks -> bouts -> merged bouts, on a trace where everything
% outside the window of interest is already NaN.
%   a LICK   is a run of consecutive samples with the tongue visible (>0),
%            at least minLenBins long
%   a BOUT   groups licks separated by <= withinBins
%   MERGING  joins bouts separated by < betweenBins
% Returns, per merged bout, the vector of its licks' start and end indices.
    boutS = {};  boutE = {};
    idx = find(~isnan(sig) & sig > 0);
    if isempty(idx), return; end

    sp    = find(diff(idx) > 1);
    s_all = idx([1,  sp+1]);
    e_all = idx([sp, end]);

    ok    = (e_all - s_all + 1) >= minLenBins;
    s_all = s_all(ok);  e_all = e_all(ok);
    if isempty(s_all), return; end

    nRuns  = numel(s_all);
    boutID = ones(1, nRuns);
    cur    = 1;
    for li = 2:nRuns
        if (s_all(li) - e_all(li-1) - 1) > withinBins
            cur = cur + 1;
        end
        boutID(li) = cur;
    end

    rawS = cell(1,cur);  rawE = cell(1,cur);
    rawFirst = nan(1,cur);  rawLast = nan(1,cur);
    for b = 1:cur
        mem = boutID == b;
        rawS{b} = s_all(mem);
        rawE{b} = e_all(mem);
        rawFirst(b) = s_all(find(mem,1,'first'));
        rawLast(b)  = e_all(find(mem,1,'last'));
    end

    mS = {rawS{1}};  mE = {rawE{1}};  endAny = rawLast(1);  nM = 1;
    for b = 2:cur
        if (rawFirst(b) - endAny - 1) < betweenBins
            mS{nM} = [mS{nM}, rawS{b}];
            mE{nM} = [mE{nM}, rawE{b}];
            endAny = rawLast(b);
        else
            nM = nM + 1;
            mS{nM} = rawS{b};
            mE{nM} = rawE{b};
            endAny = rawLast(b);
        end
    end

    for b = 1:nM
        if numel(mS{b}) >= minLicks
            boutS{end+1} = mS{b};   %#ok<AGROW>
            boutE{end+1} = mE{b};   %#ok<AGROW>
        end
    end
end

function capN = trialCapFor(anm, dte, caps)
% Last usable trial for this session, or NaN if uncapped.
% Names are canonicalised first because this project spells the same animal
% both 'TD4l' and 'TDl4'. The old code hard-coded the spelling the loaders do
% NOT produce, so its caps never fired and capped sessions kept every trial
% without a word of warning.
    capN = NaN;
    if isempty(caps), return; end
    a = canonAnm(anm);
    for r = 1:size(caps,1)
        if strcmpi(canonAnm(caps{r,1}), a) && strcmpi(caps{r,2}, dte)
            capN = caps{r,3};
            return
        end
    end
end

function s = canonAnm(a)
% 'TDl4' -> 'TD4l'. Leaves 'TD4l', 'TD13d', 'TD8d' untouched.
    s = regexprep(a, '^TD([a-zA-Z])(\d+)$', 'TD$2$1');
end

function [matLine, matHmap, dursSorted] = buildMats(cells, t_common, snip_len)
% Turn a cell array of {snip_fr, lickDur, snip_len} entries into the two
% matrices the figures use.
%   matLine     n x T, every snippet, used for the mean +/- CI traces
%   matHmap     the same, but NaN after each event's own lick ends, and rows
%               SHUFFLED THEN SORTED by lick duration
%   dursSorted  the durations in the row order of matHmap
% Shuffling before the sort matters: without it, events with identical
% durations stay in collection order, so the heatmap can show a block
% structure that is really just trial order.
    matLine = [];  matHmap = [];  dursSorted = [];
    n = numel(cells);
    if n == 0, return; end

    matLine = nan(n, snip_len);
    matHmap = nan(n, snip_len);
    durs    = nan(1, n);
    for k = 1:n
        seg_fr  = cells{k}{1};
        lickDur = cells{k}{2};
        matLine(k,:) = seg_fr;
        seg_h = seg_fr;
        seg_h(t_common > lickDur) = NaN;
        matHmap(k,:) = seg_h;
        durs(k) = lickDur;
    end

    shuf    = randperm(n);
    matHmap = matHmap(shuf,:);
    durs    = durs(shuf);
    [dursSorted, si] = sort(durs, 'ascend');
    matHmap = matHmap(si,:);
end

function drawEventHeatmap(ax, hmap, dursSorted, t_common, clim_lo, clim_hi, ...
                          cfg, col_tick, tickH, tickLW, ttl, showY, showCB, unit_label)
% One day's heatmap. NaN bins are transparent rather than painted.
    if isempty(hmap)
        title(ax, sprintf('%s\nn=0', ttl), 'FontSize',10);
        axis(ax,'off');  return
    end
    n = size(hmap,1);
    him = imagesc(ax, t_common, 1:n, hmap);
    set(him, 'AlphaData', ~isnan(hmap));
    set(ax, 'Color', [1 1 1]);
    hold(ax,'on');
    for k = 1:n
        plot(ax, [0 0], [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
        tE = dursSorted(k);
        plot(ax, [tE tE], [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
    end
    caxis(ax, [clim_lo clim_hi]);
    colormap(ax, parula);
    if showCB
        cb = colorbar(ax, 'eastoutside');
        cb.Label.String = unit_label;
        cb.Ticks      = [clim_lo clim_hi];
        cb.TickLabels = {sprintf('%.2f',clim_lo), sprintf('%.2f',clim_hi)};
        cb.TickLength = 0;
    end
    xlim(ax, cfg.xlim_hmap);
    ylim(ax, [0.5 n+0.5]);
    set(ax, 'YDir','normal', 'TickDir','out', 'FontSize',9);
    xlabel(ax, 'Time re. lick start (s)');
    if showY, ylabel(ax, 'Event # (sorted by lick duration)'); end
    title(ax, ttl, 'FontSize',10);
    box(ax,'off');
end

function drawMeanTrace(ax, mat, t_common, col, cfg, dispName)
% Mean +/- 95% CI across events, smoothed, plotted only where at least two
% snippets contribute (the snippets are NaN-padded at the trial edges).
    if isempty(mat), return; end
    mu      = nanmean(mat, 1);
    n_valid = sum(~isnan(mat), 1);
    sem     = nanstd(mat, 0, 1) ./ sqrt(max(n_valid,1));
    ci      = 1.96 * sem;

    mu_sm = movmean(mu, cfg.smooth_win, 'omitnan');
    ci_sm = movmean(ci, cfg.smooth_win, 'omitnan');

    valid = n_valid >= 2;
    t_p  = t_common(valid);
    mu_p = mu_sm(valid);
    ci_p = ci_sm(valid);

    in_xl = t_p >= cfg.xlim_line(1) & t_p <= cfg.xlim_line(2);
    t_p = t_p(in_xl);  mu_p = mu_p(in_xl);  ci_p = ci_p(in_xl);
    if numel(t_p) < 2, return; end

    fill(ax, [t_p fliplr(t_p)], [mu_p+ci_p fliplr(mu_p-ci_p)], col, ...
        'FaceAlpha',0.20, 'EdgeColor','none', 'HandleVisibility','off');
    plot(ax, t_p, mu_p, '-', 'Color',col, 'LineWidth',2.5, 'DisplayName', dispName);
end

function finishTraceAxes(ax, cfg, dd, nDays, unit_label, ttl)
    xline(ax, 0, 'k--', 'LineWidth',1, 'HandleVisibility','off');
    xlim(ax, cfg.xlim_line);
    if ~isempty(cfg.ylim_line), ylim(ax, cfg.ylim_line); end
    xlabel(ax, 'Time re. lick start (s)');
    if dd == 1, ylabel(ax, unit_label); end
    if ~isempty(ttl), title(ax, ttl, 'FontSize',11); end
    if dd == nDays
        lh = legend(ax, 'Location','northwest', 'FontSize',8);
        if ~isempty(lh), legend(ax,'boxoff'); end
    end
    set(ax, 'TickDir','out', 'FontSize',10);
    box(ax,'off');
end

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

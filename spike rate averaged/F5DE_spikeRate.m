%% F5DE_spikeRate.m
%  Session-normalized mean spike rate across units, Simple Reward Task, first five days of training.
%  Rates are averaged within a session and normalized, then averaged across
%  sessions and plotted against mean tongue length. Produces Fig. 5D-E.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; Data = dataRoot in setPaths.m
%  ANALYSIS SETTINGS
%    params.alignEvent  'firstLick'
%    params.dt          1/200
%    params.smooth      30
%    params.quality     {'good','excellent'}
%    params.lowFR       0.01
%    params.window      -1.5 to 4 s
%  Run the whole file.

clear; clc;

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.

%% PATHS

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
params.alignEvent          = 'firstLick';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

params.behav_only = 0;
params.timeWarp            = 0;   % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20;   % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.01;   % minimum mean firing rate, Hz

params.condition(1) = {'hit==1 | hit==0' };
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 4'};
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};
params.condition(end+1) = {'hit==1' };


params.tmin = -1.5;
params.tmax = 4;
params.dt     = 1/200;   % 5 ms bins

% smooth with causal gaussian kernel
params.smooth = 30;   % causal gaussian, in bins (neural AND tongue)

% cluster qualities to use
params.quality = {'good','excellent',' good','good '};   % good + excellent units; exactly the export's list, so slimToLegacy uses the exported units directly

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};
params.feat_varToExplain = 80;   % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80;   % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;

% Params for finding kinematic modes
params.fcut = 10;   % smoothing cutoff frequency
params.cond = 5;   % which conditions to use to find mode
params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
params.fa = false;   % if true, reduces neural dimensions to 10 with factor analysis
params.bctype = 'reflect';   % options are : reflect  zeropad  none

%% EPOCH / CONTACT / TEST CONFIG
%                            CONTROL PANEL
cfg = struct();

% --- the epochs ----------------------------------------------------------
cfg.epochContacts = [1 3; 3 5; 5 8];   % [firstContact lastContact] per segment
cfg.epochLabels   = {'Rest','Licks 1-3','Licks 3-5','Licks 5-8'};

% --- trials --------------------------------------------------------------
cfg.trialOutcome    = 'all';   % 'all' = hit==1 | hit==0   'hit' = hits only
cfg.rewardedLick    = 1;   % which rewardedLick type is the main trial set
cfg.minLicks        = 8;   % need >= this many contacts ...   (0 = no filter)
cfg.minLickWindow_s = 1.5;   % ... within this many seconds of the go cue

% --- what the figure shows -----------------------------------------------
cfg.usePop  = 1;   % 1 = raw population mean
% 2 = each neuron min-maxed to 0-1, THEN averaged
cfg.useRest = 3;   % 1 = Rest [-1.5 -1.0] measured from the first contact
% 2 = Rest [-1.0 -0.5] measured from that trial's go cue
% 3 = Rest [-0.3  0 ] measured from that trial's go cue,
%     i.e. the 300 ms immediately before the cue
cfg.sessionNorm = 'minmax4';
% how each session is put on a common scale:
%   'traceMinMax' normalize the whole SESSION TRACE to [0 1] first, then read
%                 the four epoch values off that scale. None of the four is
%                 pinned, so both the spacing and the overall level survive.
%   'minmax4'     min-max the FOUR VALUES alone:
%                     (E - min(E)) / (max(E) - min(E))
%                 smallest -> 0, largest -> 1. BOTH ends are pinned, and the
%                 zero is whichever of the four happens to be smallest.
%   'peak4'       divide the four values by their own peak:  E / max(E)
%                 Largest -> 1, NO floor subtracted, so Rest keeps its real
%                 size relative to the peak and the points stay proportional
%                 to each other. Only one end is pinned.
%   'restToPeak'  Rest -> 0, biggest segment -> 1
%   'none'        leave the raw values alone
cfg.traceNormPct = [0 100];
% percentiles of the trace used as the [0 1] anchors for 'traceMinMax'; [0 100]
% is a straight min-max, [1 99] limits the effect of a single noisy sample.
cfg.showLegend = true;

% --- neurons -------------------------------------------------------------
cfg.minNeuronRate = 0;   % Hz. Drop units quieter than this (0 = keep all).

% --- statistics ----------------------------------------------------------
cfg.testPoints = [3 4];   % which of the 4 points to test (3 and 4 = the segments)
cfg.testDays   = [1 5];   % which two day-groups to compare

% --- Day A vs Day B INSIDE each animal, trials as units --------
% Each trial's four values are put on its SESSION's scale (the same
% cfg.sessionNorm map that turns the session's four means into the plotted
% point, so the trial values average to exactly the plotted value). Then, per
% animal and per point, Wilcoxon rank-sum Day A trials vs Day B trials,
% one-tailed, H1: Day B LOWER than Day A. Stars go on the figure.
cfg.paTest      = true;   % false = skip it
cfg.paSets      = {[3 4]};   % one test per entry; each trial's value = MEAN of these points
% {[3 4]} = Licks 3-5 and 5-8 combined (ONE test per animal)
% {2, 3, 4} = the three segments tested separately
cfg.paBH        = true;   % BH across cfg.paSets inside each animal (no effect with one set)
cfg.paStarRule  = 'all';   % star when significant in: 'all' tested animals | 'majority'
cfg.paShowCount = true;   % print k/N above every tested set
cfg.paAlpha     = 0.05;

% Definitions behind the control panel
cfg.popModes  = {'raw','minmax01'};
cfg.popLabels = {'Raw population mean', 'Each neuron 0-1, then averaged'};

cfg.restWindows = {[-1.50 -1.00], [-1.00 -0.50], [-0.30 0.00]};
cfg.restRefs    = {'align',       'goCue',       'goCue'};
cfg.restLabels  = {'Rest [-1.5 -1.0] re: first contact', ...
                   'Rest [-1.0 -0.5] re: go cue', ...
                   'Rest [-0.3 0] re: go cue'};
cfg.restWindow_s = cfg.restWindows{cfg.useRest};

cfg.lickFields = {'lickL'};   % {'lickL','lickR'} to pool both spouts
cfg.alignEvent = params.alignEvent;

% NOTES
% epochContacts: a segment runs from the time of its first contact to the time
%   of its last contact ON THAT TRIAL, so the bins follow the animal's own lick
%   rate rather than assuming ~6.7 Hz. Contact 1 is the first contact at or
%   after the go cue. A trial needs max(epochContacts(:)) contacts to be scored.
% usePop: 'raw' is the plain mean over neurons, so high-rate units dominate it.
%   'minmax01' min-maxes each neuron's trial-averaged trace to [0 1] and applies
%   that same map to all its trials before averaging, so every neuron
%   contributes on the same bounded scale. Both are computed in the main loop
%   and stored, so switching usePop needs only a re-run of the figure cell.
% useRest: windows 2 and 3 are measured from each trial's OWN go cue, so they
%   track reaction time instead of assuming one. All three are computed and
%   stored, so switching useRest needs only a re-run of the figure cell.
%   Window 3 (300 ms before the cue) is closest to the licking but may include
%   anticipatory activity, which can grow with training.
% sessionNorm: 'minmax4' pins the SMALLEST of the four values to 0 -- which is
%   not necessarily Rest. 'restToPeak' pins Rest to 0 by construction and lets
%   later points go negative if they fall below baseline.
% minNeuronRate: optional rate floor on top of params.lowFR. Near-silent units
%   are harmless under 'raw' and 'minmax01'.
% --- run control -------------------------------------------------------
cfg.traceDays   = 1:5;   % days to draw the rate-vs-tongue figure for
cfg.traceStyles = {'A'};   % only style A is drawn (auto-scaled axes, mean +/- SEM)
cfg.showMethodComparison = true;
cfg.useCache    = true;
cfg.cacheDir    = fullfile(cfgPaths.cacheRoot, 'spikeRate_learning_cache');

%% SESSIONS (20 sessions: 4 animals x 5 training days)
datapth = '';   % raw data folder not used
spec.name = 'Learning';
spec.sessionDates = { ...
    '2025-02-03','2025-02-04','2025-02-05','2025-02-06','2025-02-07', ...   % TD3l
    '2025-02-03','2025-02-04','2025-02-05','2025-02-06','2025-02-07', ...   % TD2l
    '2025-06-03','2025-06-04','2025-06-05','2025-06-06','2025-06-07', ...   % TD4l
    '2025-08-21','2025-08-22','2025-08-23','2025-08-24','2025-08-25' };   % TD5l
spec.sessionLoaders = { ...
    @loadTD3l_neur, @loadTD3l_neur,    @loadTD3l_neur, @loadTD3l_neur, @loadTD3l_neur, ...
    @loadTD2l_neur, @loadTD2l_neur219, @loadTD2l_neur, @loadTD2l_neur, @loadTD2l_neur, ...
    @loadTD4l_neur, @loadTD4l_neur,    @loadTD4l_neur, @loadTD4l_neur, @loadTD4l_neur, ...
    @loadTD5l_neur, @loadTD5l_neur,    @loadTD5l_neur, @loadTD5l_neur, @loadTD5l_neur };

% Day maps: probe number per session for each training day (0 = unused)
s1 = [1 0 0 0 0   1 0 0 0 0   2 0 0 0 0   2 0 0 0 0 ];
s2 = [0 2 0 0 0   0 2 0 0 0   0 1 0 0 0   0 1 0 0 0 ];
s3 = [0 0 1 0 0   0 0 1 0 0   0 0 1 0 0   0 0 2 0 0 ];
s4 = [0 0 0 1 0   0 0 0 2 0   0 0 0 1 0   0 0 0 1 0 ];
s5 = [0 0 0 0 1   0 0 0 0 1   0 0 0 0 1   0 0 0 0 1 ];
dayNames = {'Day 1','Day 2','Day 3','Day 4','Day 5'};
dayCols  = [linspace(0,1,5)' zeros(5,2)];   % day 1 black -> day 5 red

nSessAll = numel(spec.sessionDates);
assert(numel(spec.sessionLoaders) == nSessAll && all(cellfun(@numel, {s1,s2,s3,s4,s5}) == nSessAll), ...
    'session table and day maps disagree in length');

% per-session accumulators
allP1         = cell(1, nSessAll);
allP2         = cell(1, nSessAll);
allmoveP1     = cell(1, nSessAll);
allmoveP4     = cell(1, nSessAll);
numNeurons    = cell(nSessAll, 2);
numTrials     = cell(1, nSessAll);
percCompleted = cell(1, nSessAll);
allE1    = cell(1, nSessAll);
allE2    = cell(1, nSessAll);
epInfo1  = cell(1, nSessAll);
epInfo2  = cell(1, nSessAll);
sessAnm  = cell(1, nSessAll);
sessDate = cell(1, nSessAll);

%% MAIN LOOP
for sessnum = 1:nSessAll
    fprintf('Session %d\n', sessnum);

[obj, pp, tongueAll] = srLoadSession(spec.sessionLoaders{sessnum}, spec.sessionDates{sessnum}, ...
                                     spec, params, cfg, datapth);
params.cluid   = pp.cluid;
params.trialid = pp.trialid;

trialSet = [1:obj.bp.Ntrials]';

condtrix = trialSet;
if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    condtrix(condtrix > 278) = [];
elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
    condtrix(condtrix > 313) = [];
elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
    condtrix(condtrix > 298) = [];
elseif strcmp(obj.pth.dt,'2025-06-05') && strcmp(obj.pth.anm,'TDl4')
    condtrix(condtrix > 123) = [];
end

Length = tongueAll(:, condtrix);   % = kin.dat(:,condtrix,kinix) for 'tongue_length'

all11 = [1:obj.bp.Ntrials];

switch lower(cfg.trialOutcome)
    case 'hit', outcomeMask = (obj.bp.hit == 1);
    case 'all', outcomeMask = true(size(obj.bp.hit));
    otherwise,  error('cfg.trialOutcome must be ''all'' or ''hit''.');
end
hit11 = all11(outcomeMask);
r111 = all11((obj.bp.rewardedLick == 1));
r444 = all11((obj.bp.rewardedLick == 4));

allr11 = intersect(hit11,r111)';
allr44 = intersect(hit11,r444)';

nHitSel    = sum(obj.bp.hit(allr11) == 1);
nNonHitSel = numel(allr11) - nHitSel;

P8 = allr11';
P9 = allr44';

if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
P9(P9 > 278) = [];
elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
P8(P8 > 313) = [];
P9(P9 > 313) = [];
elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
P8(P8 > 298) = [];
P9(P9 > 298) = [];
elseif strcmp(obj.pth.dt,'2025-06-05') && strcmp(obj.pth.anm,'TDl4')
P8(P8 > 123) = [];
P9(P9 > 123) = [];
end

[P8, nDrop8, nKeep8] = applyMinLickFilter(obj, P8, cfg);
[P9, nDrop9, nKeep9] = applyMinLickFilter(obj, P9, cfg);

if isempty(P8)
    warning('Session %d (%s %s): no trials survive the min-lick filter -- skipping session.', ...
        sessnum, obj.pth.anm, obj.pth.dt);
    continue
end

neurall = {};
if strcmp(obj.pth.dt,'2024-07-13') && strcmp(obj.pth.anm,'TD10si') || strcmp(obj.pth.dt,'2024-07-09') && strcmp(obj.pth.anm,'TD9si') ...
        || strcmp(obj.pth.dt,'2025-02-22') && strcmp(obj.pth.anm,'TDl3')  || strcmp(obj.pth.dt,'2025-02-21') && strcmp(obj.pth.anm,'TDl2') ...
        || strcmp(obj.pth.dt,'2025-02-19') && strcmp(obj.pth.anm,'TDl2')
    Ncells = size(obj.psth, 2);
    clu_m1TJ = 1:size(params.cluid,1);
    clu_ALM  = [];
    reg  = clu_m1TJ;
    reg1 = clu_m1TJ;
else
    Ncells = size(obj.psth, 2);
    clu_m1TJ = 1:size(params.cluid{1, 1},1);
    clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;
    reg  = clu_m1TJ;
    reg1 = clu_ALM;
end
neurall{1} = [reg];
neurall{2} = [reg1];

rang = [1:numel(obj.time)];

% tongue trace (identical for both probes, so computed once)
traj_1 = Length(:,P8);
traj_1 = movmean(traj_1,1);
traj_1(isnan(traj_1)) = 0;
traj_1 = mySmooth(traj_1, params.smooth, params.bctype);   % same causal kernel as the spikes, per trial
traj_1 = mean(traj_1,2); traj_1 = traj_1(rang);
traj_1 = normalize(traj_1,'range',[0 1]);
all_traj = [traj_1];
all_norm = rescale(all_traj, 0, 1);
n1 = size(traj_1,1);
traj_1 = all_norm(1:n1, :);

for k = 1:2
    neurons = neurall{k};
    aa = obj.trialdat(:,neurons,P8);

    mean_spike_rate_1 = squeeze(mean(mean(aa, 2), 3));
    mean_spike_rate_1 = mean_spike_rate_1(rang);
    mean_spike_rate_1 = movmean(mean_spike_rate_1,1);

    all_rates = [mean_spike_rate_1];
    mn = min(all_rates);
    mx = max(all_rates);
    if mx == mn
        error('All values are identical – cannot normalize to [0,1].');
    end
    all_norm = (all_rates - mn) / (mx - mn);
    mean_spike_rate_1 = all_norm(1:numel(mean_spike_rate_1));

    if k == 1
        allP1{sessnum} = [mean_spike_rate_1];
    else
        allP2{sessnum} = [mean_spike_rate_1];
    end

    [rawEp, epInfo] = epochMeansBothWays(obj, neurons, P8, cfg);
    if any(epInfo.nRestOutOfRange > 0)
        logf('             !! Rest window off the end of the trace on');
        for w = 1:numel(epInfo.nRestOutOfRange)
            if epInfo.nRestOutOfRange(w) > 0
                logf(' %d trial(s) for window %d', epInfo.nRestOutOfRange(w), w);
            end
        end
        logf(' -- params.tmin may be too late.\n');
    end
    if k == 1
        allE1{sessnum}   = rawEp;
        epInfo1{sessnum} = epInfo;
    else
        allE2{sessnum}   = rawEp;
        epInfo2{sessnum} = epInfo;
    end
end

allmoveP1{sessnum} = [traj_1'];
allmoveP4{sessnum} = [traj_1'];

numNeurons{sessnum,1} = [numel(clu_m1TJ)];
numNeurons{sessnum,2} = [numel(clu_ALM)];

sessAnm{sessnum}  = obj.pth.anm;
sessDate{sessnum} = obj.pth.dt;

numTrials{sessnum} = sum( obj.bp.hit);
percCompleted{sessnum} = sum( obj.bp.hit) / obj.bp.Ntrials;

clear obj tongueAll Length aa
end

allmoveP1 = allmoveP1(:);

%% RATE vs TONGUE, per day
if any(strcmpi(cfg.traceStyles, 'A'))
for pos = cfg.traceDays
    allpos = [s1;s2;s3;s4;s5];
    lw = 3.5;
    p = allpos(pos,:);
    validRow1    = [];
    validRow1mov = [];
    for i = 1:numel(p)
        idx = p(i);
        if ~idx, continue; end
        SD = eval(sprintf('allP%d', idx));
        SD = SD{i};
        if isempty(SD), continue; end   % session produced no data
        validRow1(end+1,:)    = SD(:,1);
        validRow1mov(end+1,:) = cell2mat(allmoveP1(i,:));
    end

    n1         = size(validRow1, 1);
    winSR1     = 19;
    winMov1    = 19;

    meanRow1 = movmean(mean(validRow1, 1), winSR1);

    if exist('tinv','file')
        tcrit = tinv(0.975, n1-1);
    else
        tcrit = 1.96;
    end
    semRow1    = (std(validRow1, 0, 1) ./ sqrt(n1)) * tcrit;
    semRow1    = movmean(semRow1, winSR1);   % smooth with same window as mean

    a_smooth_1 = movmean(mean(validRow1mov, 1), winMov1);

    t = (1:numel(meanRow1)) * 0.005 - 1.5;

    pad    = 0.05;
    mn_r   = min(meanRow1 - semRow1); mx_r = max(meanRow1 + semRow1);
    mn_k   = min(a_smooth_1);         mx_k = max(a_smooth_1);
    rng_r  = mx_r - mn_r;             rng_k = mx_k - mn_k;
    ylim_l = [mn_r - pad*rng_r,  mx_r + pad*rng_r];
    ylim_r = [mn_k - pad*rng_k,  mx_k + pad*rng_k];

    figure; hold on

    yyaxis left
    dcol = dayCols(pos,:);   % day 1 black -> day 5 red
    fill([t, fliplr(t)], [meanRow1+semRow1, fliplr(meanRow1-semRow1)], ...
        dcol, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(t, meanRow1, '-', 'Color', dcol, 'LineWidth', lw);
    ylabel('Mean Spike Rate');
    ylim(ylim_l);
    ax = gca; ax.YAxis(1).Color = dcol;

    yyaxis right
    plot(t, a_smooth_1, '-', 'Color', [0 0 0], 'LineWidth', lw);
    ylabel('Mean Tongue Length / Motion Energy');
    ylim(ylim_r);
    ax.YAxis(2).Color = [0 0 0];

    xline(0, 'k--', 'LineWidth', 1);
    xlim([-0.6 1]);
    xlabel('Time re: first contact (s)');
    box off
    set(gca, 'FontSize', 11)
    px = 50; py = 50; width = 350; height = 400;
    set(gcf, 'Position', [px, py, width, height]);
    set(gcf, 'Name', sprintf('%s - style A', dayNames{pos}));   % window name only
end
end

%% THE FIGURE
% One panel. What it shows is set by cfg.usePop, cfg.useRest and
% cfg.sessionNorm in the control panel at the top of the file -- change those
% and re-run this cell; no need to re-run the main loop.
% Pipeline, per session and probe:
%   neurons' activity -> per-trial segment bounds from that trial's own contact
%   times -> average across the whole population inside each interval ->
%   average over trials -> four raw values -> normalize those four
%   (cfg.sessionNorm) -> mean +/- SEM across sessions within a day.

region_codes = {s1,s2,s3,s4,s5};
region_names = {'Day 1','Day 2','Day 3','Day 4','Day 5'};
colors = dayCols;   % day 1 black -> day 5 red

x = 1:numel(cfg.epochLabels);
[GRP, GANM] = collectGridGroups(region_codes, allE1, allE2, sessAnm, ...
                                cfg.usePop, cfg.useRest, cfg.sessionNorm);

% ---- Day A vs Day B inside each animal, trials as units ----
PA = struct('anm',{{}}, 'p',[], 'q',[], 'nA',[], 'nB',[], 'medA',[], 'medB',[], 'k',[], 'n',[]);
if cfg.paTest
    gA = cfg.testDays(1);  gB = cfg.testDays(2);
    [TA, AA] = collectTrialGroups(region_codes{gA}, allE1, allE2, sessAnm, cfg.usePop, cfg.useRest, cfg.sessionNorm);
    [TB, AB] = collectTrialGroups(region_codes{gB}, allE1, allE2, sessAnm, cfg.usePop, cfg.useRest, cfg.sessionNorm);
    PA.anm = unique([AA(:); AB(:)], 'stable');
    nAn = numel(PA.anm);  nPt = numel(cfg.paSets);
    PA.p = nan(nAn,nPt);  PA.nA = zeros(nAn,nPt);  PA.nB = zeros(nAn,nPt);
    PA.medA = nan(nAn,nPt);  PA.medB = nan(nAn,nPt);
    for ai = 1:nAn
        VA = [TA{strcmp(AA, PA.anm{ai})}];   % all Day-A trials of this animal
        VB = [TB{strcmp(AB, PA.anm{ai})}];
        for ii = 1:nPt
            pt = cfg.paSets{ii};
            a = zeros(0,1);  b = zeros(0,1);
% each trial -> mean of its values at the points of this set (all must exist)
            if ~isempty(VA), a = mean(VA(pt,:), 1)'; a = a(isfinite(a)); end
            if ~isempty(VB), b = mean(VB(pt,:), 1)'; b = b(isfinite(b)); end
            PA.nA(ai,ii) = numel(a);  PA.nB(ai,ii) = numel(b);
            if ~isempty(a), PA.medA(ai,ii) = median(a); end
            if ~isempty(b), PA.medB(ai,ii) = median(b); end
            if numel(a) >= 2 && numel(b) >= 2
                PA.p(ai,ii) = ranksum(a, b, 'tail', 'right');   % Day A > Day B
            end
        end
    end
    PA.q = PA.p;
    if cfg.paBH
        for ai = 1:nAn, PA.q(ai,:) = bhAdjustPA(PA.p(ai,:)); end
    end
    PA.k = sum(PA.q < cfg.paAlpha, 1);
    PA.n = sum(isfinite(PA.q), 1);

    fprintf('\n=================== DAY %d vs DAY %d INSIDE EACH ANIMAL ===================\n', gA, gB);
    fprintf('%s | %s | sessionNorm %s (each trial on its session''s scale)\n', ...
        cfg.popLabels{cfg.usePop}, cfg.restLabels{cfg.useRest}, cfg.sessionNorm);
    fprintf('Wilcoxon rank-sum, Day %d trials vs Day %d trials, one-tailed H1: Day %d LOWER%s\n', ...
        gA, gB, gB, ternaryPA(cfg.paBH && nPt > 1, ' | BH across sets within animal', ''));
    for ii = 1:nPt
        pt = cfg.paSets{ii};
        fprintf('  --- per-trial mean of point(s) %s (%s) ---\n', mat2str(pt), strjoin(cfg.epochLabels(pt), ' + '));
        fprintf('    %-10s %7s %7s %10s %10s %10s %10s\n', 'animal', ...
            sprintf('n d%d',gA), sprintf('n d%d',gB), sprintf('med d%d',gA), sprintf('med d%d',gB), 'p', 'q');
        for ai = 1:nAn
            st = '';  if isfinite(PA.q(ai,ii)) && PA.q(ai,ii) < cfg.paAlpha, st = '  *'; end
            fprintf('    %-10s %7d %7d %10.3f %10.3f %10.4g %10.4g%s\n', PA.anm{ai}, PA.nA(ai,ii), ...
                PA.nB(ai,ii), PA.medA(ai,ii), PA.medB(ai,ii), PA.p(ai,ii), PA.q(ai,ii), st);
        end
        fprintf('    -> significant in %d of %d animals\n', PA.k(ii), PA.n(ii));
    end
    fprintf('=====================================================================================\n');
end

figure('Color','w','Position',[120 120 540 470]);
hold on;

hL = [];  hN = {};
for r = 1:numel(region_codes)
    R = GRP{r};
    if isempty(R), continue; end
    m   = mean(R, 2, 'omitnan');
    sem = std(R, 0, 2, 'omitnan') ./ sqrt(max(sum(~isnan(R),2),1));

    fill([x, fliplr(x)], [(m+sem)', fliplr((m-sem)')], colors(r,:), ...
         'EdgeColor','none', 'FaceAlpha',0.15);
    hL(end+1) = plot(x, m, '-o', 'Color', colors(r,:), ...   %#ok<SAGROW>
        'MarkerFaceColor', colors(r,:), 'LineWidth', 2, 'MarkerSize', 6);
    hN{end+1} = sprintf('%s (n=%d)', region_names{r}, size(R,2));   %#ok<SAGROW>
end

% ---- stars above the tested points ----
if cfg.paTest && ~isempty(PA.n)
    yTop = -inf(1, numel(x));
    for r = 1:numel(region_codes)
        R = GRP{r};
        if isempty(R), continue; end
        m   = mean(R, 2, 'omitnan');
        sem = std(R, 0, 2, 'omitnan') ./ sqrt(max(sum(~isnan(R),2),1));
        top = m + sem;  top(~isfinite(top)) = m(~isfinite(top));
        yTop = max(yTop, top(:)');
    end
    yl  = ylim;  pad = 0.05 * diff(yl);
    for ii = 1:numel(cfg.paSets)
        pt = cfg.paSets{ii};
        if PA.n(ii) == 0 || ~all(isfinite(yTop(pt))), continue; end
        switch lower(cfg.paStarRule)
            case 'majority', isSig = PA.k(ii) > PA.n(ii)/2;
            otherwise,       isSig = PA.k(ii) == PA.n(ii);
        end
        yS = max(yTop(pt)) + pad;
        xc = mean(pt);
        if numel(pt) > 1   % bracket over the combined points
            plot([min(pt) min(pt) max(pt) max(pt)], yS + [-0.4 0 0 -0.4]*pad, 'k-', ...
                'LineWidth', 1, 'HandleVisibility','off');
        end
        if isSig
            text(xc, yS, '*', 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', 18, 'FontWeight','bold', 'Color','k');
        end
        if cfg.paShowCount
            text(xc, yS + ternaryPA(isSig, 1.6, 0)*pad, sprintf('%d/%d', PA.k(ii), PA.n(ii)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', 8, 'Color', [0.35 0.35 0.35]);
        end
    end
    yl2 = ylim;
    ylim([yl2(1) max(yl2(2), max(yTop([cfg.paSets{:}])) + 4*pad)]);
end

xlim([0.5 numel(x)+0.5]);
set(gca, 'XTick', x, 'XTickLabel', cfg.epochLabels, ...
         'XTickLabelRotation', 45, 'FontSize', 11);
xlabel('Epoch');
ylabel(sprintf('Activity (%s)', cfg.sessionNorm));
nAnmFig = numel(unique(sessAnm(~cellfun(@isempty, sessAnm))));
title(sprintf('tjM1  |  %d sessions, %d animals', nSessAll, nAnmFig), ...
    'FontSize', 10, 'FontWeight', 'normal');
if cfg.showLegend && ~isempty(hL)
    legend(hL, hN, 'Location','northoutside','Orientation','horizontal', ...
           'FontSize', 8, 'Box','off');
end
box off;
hold off;

% ---- the same numbers ----------------------------------------------------
logf('\n=================== THE FIGURE, as numbers ===================\n');
logf('%s | %s | sessionNorm %s\n', cfg.popLabels{cfg.usePop}, ...
    cfg.restLabels{cfg.useRest}, cfg.sessionNorm);
logf('  %-6s', 'day');
logf('%14s', cfg.epochLabels{:});
logf('%8s\n', 'n');
for r = 1:numel(region_codes)
    R = GRP{r};
    if isempty(R), continue; end
    logf('  %-6d', r);
    logf('%14.3f', mean(R, 2, 'omitnan'));
    logf('%8d\n', size(R,2));
end
logf('==============================================================\n');

%% ROBUSTNESS SWEEP
% The day-1 minus day-5 gap at each test point, under EVERY plot-time
% combination of binning method, floor and scale.
% A gap with the same sign and similar size in every row does not depend on
% the analysis choice.
% The population mode and Rest window are not swept here.

region_codes = {s1,s2,s3,s4,s5};
methods = {'perTrial','traceThenBin','traceThenBinAll'};
floors  = {'rest','traceMin'};
scales  = {'peakSegment','fixedSegment'};
gA = cfg.testDays(1);  gB = cfg.testDays(2);

% Local copy; every field below is swept.
dcfg = cfg;
dcfg.epochMethod  = 'perTrial';
dcfg.floorMode    = 'rest';
dcfg.scaleMode    = 'peakSegment';
dcfg.scaleSegment = 1;

logf('\n================= ROBUSTNESS SWEEP =================\n');
logf('day %d minus day %d, per-animal paired, popMode = %s\n\n', gA, gB, cfg.popModes{1});
logf('%-18s %-10s %-14s', 'method', 'floor', 'scale');
for pt = cfg.testPoints
    logf('%16s', cfg.epochLabels{pt});
end
logf('%8s\n', 'n');

for mi = 1:numel(methods)
  for fi = 1:numel(floors)
    for si = 1:numel(scales)
        c2 = dcfg;  c2.scaleMode = scales{si};
        try
            [G, GA] = collectEpochGroups(region_codes, allE1, allE2, sessAnm, c2, ...
                                         methods{mi}, floors{fi});
        catch
            continue
        end
        if isempty(G{gA}) || isempty(G{gB}), continue; end
        logf('%-18s %-10s %-14s', methods{mi}, floors{fi}, scales{si});
        nShow = 0;
        for pt = cfg.testPoints
            [aN, aV] = perAnimalMean(GA{gA}, G{gA}(pt,:));
            [bN, bV] = perAnimalMean(GA{gB}, G{gB}(pt,:));
            [~, ia, ib] = intersect(aN, bN, 'stable');
            d  = aV(ia) - bV(ib);
            d  = d(isfinite(d));
            nShow = numel(d);
            if isempty(d), logf('%16s', 'n/a');
            else,          logf('%16.3f', mean(d));
            end
        end
        logf('%8d\n', nShow);
    end
  end
end
logf('\n  A positive number means day %d sits BELOW day %d at that point.\n', gB, gA);
logf('====================================================\n');

%% EPOCH DIAGNOSTICS (no tests; Fig. 5E is the within-animal rank-sum above)

region_codes = {s1,s2,s3,s4,s5};
[GRP, GANM]  = collectGridGroups(region_codes, allE1, allE2, sessAnm, ...
                                 cfg.usePop, cfg.useRest, cfg.sessionNorm);

gA = cfg.testDays(1);
gB = cfg.testDays(2);

logf('\n=================== EPOCH DIAGNOSTICS ===================\n');
logf('trials        : outcome = %s, rewardedLick == %d\n', cfg.trialOutcome, cfg.rewardedLick);
logf('population    : %s\n', cfg.popLabels{cfg.usePop});
logf('rest window   : %s\n', cfg.restLabels{cfg.useRest});
logf('sessionNorm   : %s\n', cfg.sessionNorm);
logf('segments      : contacts %s\n', mat2str(cfg.epochContacts));
logf('trial filter  : >= %d contacts within %.2f s of the go cue\n', cfg.minLicks, cfg.minLickWindow_s);
logf('sessions/day  : ');
logf('%d ', cellfun(@(g) size(g,2), GRP));
logf('\n\n');

% ---- lick-rate check: is the confound present? --------------------------
logf('--- median segment duration by day (s) -- checks the lick-rate confound ---\n');
nSeg = size(cfg.epochContacts,1);
for g = 1:numel(region_codes)
    D = segDurByGroup(region_codes{g}, epInfo1, epInfo2);
    if isempty(D), continue; end
    md = mean(D, 2, 'omitnan');
    nContact = cfg.epochContacts(:,2) - cfg.epochContacts(:,1);
    logf('  day %d :', g);
    for e = 1:nSeg
        logf('  %s = %.3f s (ICI %.3f)', cfg.epochLabels{e+1}, md(e), md(e)/max(nContact(e),1));
    end
    logf('\n');
end
logf('  -> if inter-contact interval (ICI) shifts across days, fixed-time bins\n');
logf('     would confound lick rate with disengagement; contact-based bins do not.\n\n');

% ---- where each session's peak segment sits ------------------------------
% Under 'restToPeak' or 'minmax4' one point is pinned to 1 in every session, and
% the group mean reaches exactly 1 only if EVERY session peaks in the SAME
% segment; when some peak later that column averages below 1. Under
% 'traceMinMax' nothing is pinned, so no column need be 0 or 1 at all. Either
% way this census says which segment carries each session's maximum.
logf('--- which segment holds each session''s peak ---\n');
logf('  %-6s', 'day');
logf('%14s', cfg.epochLabels{2:end});
logf('\n');
for g = 1:numel(region_codes)
    R = GRP{g};
    if isempty(R), continue; end
    [~, pk] = max(R(2:end,:), [], 1, 'omitnan');
    logf('  %-6d', g);
    for e = 1:(size(R,1)-1)
        logf('%14s', sprintf('%d/%d', sum(pk==e), size(R,2)));
    end
    logf('\n');
end
logf('  (a column reads 1.000 in the figure only when every session peaks there)\n\n');

%% LOCAL FUNCTIONS

% CONTACT-BASED EPOCH MACHINERY

function lk = gatherLicks(obj, tr, cfg)
% All contact times for one trial, in session time, sorted.
    lk = [];
    for f = 1:numel(cfg.lickFields)
        fn = cfg.lickFields{f};
        if isfield(obj.bp.ev, fn)
            v = obj.bp.ev.(fn);
            if iscell(v) && numel(v) >= tr && ~isempty(v{tr})
                lk = [lk; v{tr}(:)];   %#ok<AGROW>
            end
        end
    end
    lk = sort(lk);
end

function t0 = alignTime(obj, tr, cfg)
% Session time of the alignment event for one trial, so contact times can be
% put on the same clock as obj.time. Reads obj.bp.ev.(params.alignEvent) when
% that field exists; otherwise falls back to the first contact at or after the
% go cue, which is what 'firstLick' means.
    t0 = nan;
    if isfield(obj.bp.ev, cfg.alignEvent)
        v = obj.bp.ev.(cfg.alignEvent);
        if iscell(v)
            if numel(v) >= tr && ~isempty(v{tr}), t0 = v{tr}(1); end
        elseif isnumeric(v)
            if numel(v) >= tr, t0 = v(tr); end
        end
    end
    if ~isfinite(t0)
        lk = gatherLicks(obj, tr, cfg);
        gc = obj.bp.ev.goCue(tr);
        lk = lk(lk >= gc - 1e-9);
        if ~isempty(lk), t0 = lk(1); end
    end
end

function [keepIdx, nDrop, nKeep] = applyMinLickFilter(obj, trials, cfg)
% Keep trials with at least cfg.minLicks contacts within cfg.minLickWindow_s
% of the go cue. cfg.minLicks = 0 keeps everything.
    trials = trials(:);
    if cfg.minLicks <= 0
        keepIdx = trials;  nDrop = 0;  nKeep = numel(trials);
        return
    end
    keep = false(numel(trials),1);
    for j = 1:numel(trials)
        tr = trials(j);
        lk = gatherLicks(obj, tr, cfg);
        gc = obj.bp.ev.goCue(tr);
        keep(j) = sum(lk >= gc - 1e-9 & lk <= gc + cfg.minLickWindow_s + 1e-9) >= cfg.minLicks;
    end
    keepIdx = trials(keep);
    nKeep   = sum(keep);
    nDrop   = numel(trials) - nKeep;
end

function [out, info] = epochMeansBothWays(obj, neurons, trials, cfg)
% RAW epoch values for one session/probe, computed for EVERY combination of
% population mode x rest window, plus the fields the robustness sweep uses.
%   out.grid{p,w}   (nSeg+1) x 1 raw epoch means for cfg.popModes{p} and
%                   cfg.restWindows{w}. Row 1 is Rest, rows 2..end the segments.
%   1. take the neurons' activity
%   2. per trial, work out the time bounds of each segment from THAT trial's
%      own contact times (and the Rest window from its own go cue, if that
%      window is go-cue referenced)
%   3. average across the whole population inside those intervals
%   4. average the per-trial values -> four numbers for this session
% Normalizing those four numbers happens at PLOT time (cfg.sessionNorm).
    nSeg   = size(cfg.epochContacts,1);
    nEp    = nSeg + 1;
    nNeed  = max(cfg.epochContacts(:));
    nPop   = numel(cfg.popModes);
    nRest  = numel(cfg.restWindows);
    trials = trials(:);

    out = struct();
    out.grid            = cell(nPop, nRest);
    out.traceRange      = nan(nPop, 2);   % [lo hi] anchors of the session trace
    out.perTrial        = nan(nEp,1);
    out.traceThenBin    = nan(nEp,1);
    out.traceThenBinAll = nan(nEp,1);
    out.traceMinMatched = nan;
    out.traceMinAll     = nan;
    out.trialS          = cell(1, nPop);   % nSeg  x nTrial, per pop mode
    out.trialR          = cell(1, nPop);   % nRest x nTrial, per pop mode

    info = struct('nNeuronIn', numel(neurons), 'nNeuronKept', 0, ...
                  'rateMin', nan, 'rateMed', nan, ...
                  'nTrialIn', numel(trials), 'nTrialWithContacts', 0, ...
                  'nContactNeeded', nNeed, 'medSegDur', nan(nSeg,1), ...
                  'meanContactT', nan(nNeed,1), 'medRT', nan, ...
                  'nRestOutOfRange', zeros(1,nRest));

    if isempty(trials) || isempty(neurons), return; end

    tt = obj.time(:);
    A0 = obj.trialdat(:, neurons, trials);   % time x nNeuron x nTrial

% ---- neuron rate floor, before any per-neuron scaling ----------------
    rate0 = mean(A0, [1 3], 'omitnan');
    rate0 = rate0(:)';
    info.rateMin = min(rate0);
    info.rateMed = median(rate0, 'omitnan');
    if cfg.minNeuronRate > 0
        A0 = A0(:, rate0 >= cfg.minNeuronRate, :);
    end
    info.nNeuronKept = size(A0,2);
    if size(A0,2) == 0
        warning('All %d neurons fell below cfg.minNeuronRate = %g Hz.', ...
            numel(neurons), cfg.minNeuronRate);
        return
    end

% ---- per-trial geometry: identical for every population mode ---------
    segLo = nan(numel(trials), nSeg);   segHi = segLo;
    resLo = nan(numel(trials), nRest);  resHi = resLo;
    CT    = nan(nNeed, numel(trials));
    dur   = nan(nSeg,  numel(trials));
    rt    = nan(1, numel(trials));
    good  = false(numel(trials),1);

    for j = 1:numel(trials)
        tr = trials(j);
        lk = gatherLicks(obj, tr, cfg);
        gc = obj.bp.ev.goCue(tr);
        lk = lk(lk >= gc - 1e-9);   % contacts at/after the go cue
        t0 = alignTime(obj, tr, cfg);
        if ~isfinite(t0) || numel(lk) < nNeed
            continue   % too few contacts to score
        end
        good(j)  = true;
        lkRel    = lk - t0;
        CT(:,j)  = lkRel(1:nNeed);
        rt(j)    = t0 - gc;   % reaction time

        for e = 1:nSeg
            segLo(j,e) = lkRel(cfg.epochContacts(e,1));
            segHi(j,e) = lkRel(cfg.epochContacts(e,2));
            dur(e,j)   = segHi(j,e) - segLo(j,e);
        end

        gcRel = gc - t0;   % go cue on the trace clock
        for w = 1:nRest
            W = cfg.restWindows{w};
            switch lower(cfg.restRefs{w})
                case 'align', ref = 0;
                case 'gocue', ref = gcRel;
                otherwise, error('cfg.restRefs{%d} must be ''align'' or ''goCue''.', w);
            end
            resLo(j,w) = ref + W(1);
            resHi(j,w) = ref + W(2);
            if resHi(j,w) < tt(1) || resLo(j,w) > tt(end)
                info.nRestOutOfRange(w) = info.nRestOutOfRange(w) + 1;
            end
        end
    end

    info.nTrialWithContacts = sum(good);
    info.medSegDur    = median(dur, 2, 'omitnan');
    info.meanContactT = mean(CT, 2, 'omitnan');
    info.medRT        = median(rt, 'omitnan');

% ---- one pass per population mode ------------------------------------
    for p = 1:nPop
        A   = applyPopMode(A0, cfg.popModes{p});
        pop = mean(A, 2, 'omitnan');
        pop = reshape(pop, size(A,1), numel(trials));   % time x nTrial

        S = nan(nSeg,  numel(trials));
        R = nan(nRest, numel(trials));
        for j = 1:numel(trials)
            if ~good(j), continue; end
            for e = 1:nSeg
                if ~(segHi(j,e) > segLo(j,e)), continue; end
                m = tt >= segLo(j,e) & tt <= segHi(j,e);
                if any(m), S(e,j) = mean(pop(m,j), 'omitnan'); end
            end
            for w = 1:nRest
                m = tt >= resLo(j,w) & tt <= resHi(j,w);
                if any(m), R(w,j) = mean(pop(m,j), 'omitnan'); end
            end
        end

        segMean = mean(S, 2, 'omitnan');
        out.trialS{p} = S;   out.trialR{p} = R;
        for w = 1:nRest
            out.grid{p,w} = [mean(R(w,:), 'omitnan'); segMean];
        end

% ---- anchors for cfg.sessionNorm = 'traceMinMax' ------------------
% The whole trial-averaged trace, over the SAME trials that produced
% the segment values. Because the map (x - lo)/(hi - lo) is affine and
% averaging is linear, normalizing the trace here and reading the epoch
% values off it is identical to normalizing every trial before binning.
        if any(good)
            traceM = mean(pop(:, good), 2, 'omitnan');
            out.traceRange(p,:) = [pctl(traceM, cfg.traceNormPct(1)), ...
                                   pctl(traceM, cfg.traceNormPct(2))];
        end

% ---- fields for the robustness sweep (raw pass, first rest window) ---
        if strcmpi(cfg.popModes{p}, 'raw')
            out.perTrial = out.grid{p,1};
            m1 = tt >= cfg.restWindows{1}(1) & tt <= cfg.restWindows{1}(2);
            if any(good)
                traceM = mean(pop(:, good), 2, 'omitnan');
                out.traceMinMatched = min(traceM);
                out.traceThenBin    = binTraceAtTimes(traceM, tt, m1, info.meanContactT, cfg);
            end
            traceA = mean(pop, 2, 'omitnan');
            out.traceMinAll = min(traceA);
            if any(isfinite(info.meanContactT))
                out.traceThenBinAll = binTraceAtTimes(traceA, tt, m1, info.meanContactT, cfg);
            end
        end
    end
end

function A = applyPopMode(A, mode)
% Transform the time x nNeuron x nTrial array BEFORE the population average.
    switch lower(mode)
        case 'raw'
% nothing -- plain mean over neurons, high-rate units dominate

        case 'minmax01'
% Min-max each neuron's TRIAL-AVERAGED trace to [0 1], then apply
% that same affine map to all of its trials. Bounded, so no neuron
% can be multiplied up into noise the way 1/mean does.
            tr1 = mean(A, 3, 'omitnan');   % time x nNeuron
            lo  = min(tr1, [], 1);
            hi  = max(tr1, [], 1);
            rgn = hi - lo;
            rgn(~isfinite(rgn) | rgn <= 0) = NaN;   % flat neuron -> dropped
            A = (A - reshape(lo,1,[],1)) ./ reshape(rgn,1,[],1);

        case 'rate'
            mu = mean(A, [1 3], 'omitnan');
            mu(mu <= 0) = NaN;
            A = A ./ reshape(mu,1,[],1);

        otherwise
            error('Unknown population mode: %s', mode);
    end
end

function v = normSessionEpochs(E, mode, traceRange)
% Put ONE session's four raw epoch values on a common scale.
%   traceRange gives the [lo hi] anchors used by 'traceMinMax'.
    if nargin < 3, traceRange = [nan nan]; end
    v = E(:);
    if all(~isfinite(v)), v = nan(size(v)); return; end
    switch lower(mode)
        case 'none'
% leave alone
        case 'traceminmax'
            lo = traceRange(1);  hi = traceRange(2);
            if isfinite(hi-lo) && (hi-lo) > 0, v = (v-lo)./(hi-lo); else, v = nan(size(v)); end
        case 'minmax4'
            lo = min(v); hi = max(v);
            if isfinite(hi-lo) && (hi-lo) > 0, v = (v-lo)./(hi-lo); else, v = nan(size(v)); end
        case 'peak4'
% Divide by the largest of the four. No floor is subtracted, so the
% four values keep their ratios to one another and Rest lands at its
% true fraction of the peak instead of being forced to zero.
            d = max(v);
            if isfinite(d) && d ~= 0, v = v ./ d; else, v = nan(size(v)); end
        case 'resttopeak'
            rest = v(1); den = max(v(2:end)) - rest;
            if isfinite(den) && den > 0, v = (v-rest)./den; else, v = nan(size(v)); end
        otherwise
            error('Unknown cfg.sessionNorm: %s', mode);
    end
end

function q = pctl(v, p)
% Linear-interpolated percentile. No Statistics Toolbox needed, and identical
% to min()/max() at p = 0 and p = 100.
    v = sort(v(isfinite(v)));
    if isempty(v), q = nan; return; end
    if numel(v) == 1, q = v; return; end
    idx = 1 + (p/100) * (numel(v) - 1);
    lo  = floor(idx);  hi = ceil(idx);
    if lo == hi
        q = v(lo);
    else
        q = v(lo) + (idx - lo) * (v(hi) - v(lo));
    end
end

function [GRP, GANM] = collectGridGroups(region_codes, allE1, allE2, sessAnm, p, w, mode)
% Gather every session's normalized epoch vector for one (popMode, restWindow).
    nG = numel(region_codes);
    GRP = cell(1,nG);  GANM = cell(1,nG);
    for g = 1:nG
        code = region_codes{g};
        R = [];  A = {};
        for sess = 1:numel(code)
            if sess > numel(allE1), break; end
            switch code(sess)
                case 1,    raw = allE1{sess};
                case 2,    raw = allE2{sess};
                otherwise, continue
            end
            if isempty(raw) || ~isfield(raw,'grid'), continue; end
            if p > size(raw.grid,1) || w > size(raw.grid,2), continue; end
            E = raw.grid{p,w};
            if isempty(E), continue; end
            rngP = [nan nan];
            if isfield(raw,'traceRange') && p <= size(raw.traceRange,1)
                rngP = raw.traceRange(p,:);
            end
            R(:,end+1) = normSessionEpochs(E, mode, rngP);   %#ok<AGROW>
            if sess <= numel(sessAnm) && ~isempty(sessAnm{sess})
                A{end+1} = char(sessAnm{sess});   %#ok<AGROW>
            else
                A{end+1} = sprintf('sess%02d', sess);   %#ok<AGROW>
            end
        end
        GRP{g} = R;  GANM{g} = A;
    end
end

function v = binTraceAtTimes(trace, tt, restM, contactT, cfg)
% Bin one already-averaged trace using a single set of contact times.
    nSeg = size(cfg.epochContacts,1);
    v    = nan(nSeg+1, 1);
    if any(restM), v(1) = mean(trace(restM), 'omitnan'); end
    for e = 1:nSeg
        lo = contactT(cfg.epochContacts(e,1));
        hi = contactT(cfg.epochContacts(e,2));
        if ~isfinite(lo) || ~isfinite(hi) || ~(hi > lo), continue; end
        m = tt >= lo & tt <= hi;
        if any(m), v(1+e) = mean(trace(m), 'omitnan'); end
    end
end

function [v, e] = normEpochs(S, cfg, methodOverride, floorOverride)
% Apply cfg.epochMethod / cfg.floorMode / cfg.scaleMode to one session's raw
% epoch struct. The two overrides let the comparison cell sweep the options
% without touching cfg.
    meth = cfg.epochMethod;   if nargin >= 3 && ~isempty(methodOverride), meth = methodOverride; end
    flm  = cfg.floorMode;     if nargin >= 4 && ~isempty(floorOverride),  flm  = floorOverride;  end

    switch lower(meth)
        case 'pertrial',        e = S.perTrial;        minField = 'traceMinMatched';
        case 'tracethenbin',    e = S.traceThenBin;    minField = 'traceMinMatched';
        case 'tracethenbinall', e = S.traceThenBinAll; minField = 'traceMinAll';
        otherwise, error('Unknown cfg.epochMethod: %s', meth);
    end

% The trace minimum is needed only for floorMode 'traceMin'.
    tmin = nan;
    if isfield(S, minField), tmin = S.(minField); end
    e = e(:);
    v = nan(size(e));
    if all(~isfinite(e)), return; end

    switch lower(flm)
        case 'rest',     fl = e(1);
        case 'tracemin'
            if ~isfield(S, minField)
                error(['cfg.floorMode = ''traceMin'' needs the trace minimum, which ' ...
                       'this stored result does not carry. Re-run the main loop once ' ...
                       '(the epoch struct saves it), or use cfg.floorMode = ''rest''.']);
            end
            fl = tmin;
        case 'none',     fl = 0;
        otherwise, error('Unknown cfg.floorMode: %s', flm);
    end
    if ~isfinite(fl), return; end

    switch lower(cfg.scaleMode)
        case 'peaksegment'
            den = max(e(2:end)) - fl;
            if isfinite(den) && den > 0, v = (e - fl) ./ den; end
        case 'fixedsegment'
            den = e(1 + cfg.scaleSegment) - fl;
            if isfinite(den) && den > 0, v = (e - fl) ./ den; end
        case 'none'
            v = e - fl;
        otherwise, error('Unknown cfg.scaleMode: %s', cfg.scaleMode);
    end
end

function [GRP, GANM] = collectEpochGroups(region_codes, allE1, allE2, sessAnm, cfg, methodOverride, floorOverride)
% For each group (day), the normalized epoch vector of every session in it,
% plus that session's animal name (used for the paired test).
    if nargin < 6, methodOverride = ''; end
    if nargin < 7, floorOverride  = ''; end
    nG   = numel(region_codes);
    GRP  = cell(1,nG);
    GANM = cell(1,nG);

    for g = 1:nG
        code = region_codes{g};
        R = [];  A = {};
        for sess = 1:numel(code)
            if sess > numel(allE1), break; end
            switch code(sess)
                case 1,    raw = allE1{sess};
                case 2,    raw = allE2{sess};
                otherwise, continue
            end
            if isempty(raw), continue; end

            R(:,end+1) = normEpochs(raw, cfg, methodOverride, floorOverride);   %#ok<AGROW>
            if sess <= numel(sessAnm) && ~isempty(sessAnm{sess})
                A{end+1} = char(sessAnm{sess});   %#ok<AGROW>
            else
                A{end+1} = sprintf('sess%02d', sess);   %#ok<AGROW>
            end
        end
        GRP{g}  = R;
        GANM{g} = A;
    end
end

function D = segDurByGroup(code, epInfo1, epInfo2)
% Median segment durations (s) for every session in one group, one column per session.
    D = [];
    for sess = 1:numel(code)
        if sess > numel(epInfo1), break; end
        switch code(sess)
            case 1,    inf_s = epInfo1{sess};
            case 2,    inf_s = epInfo2{sess};
            otherwise, continue
        end
        if isempty(inf_s) || ~isfield(inf_s, 'medSegDur'), continue; end
        D(:,end+1) = inf_s.medSegDur(:);   %#ok<AGROW>
    end
end

function [T, A] = collectTrialGroups(code, allE1, allE2, sessAnm, p, w, mode)
% per-trial epoch values for one day map: T{i} is nEp x nTrial for
% session i, each trial put on that session's cfg.sessionNorm scale (the same
% affine map normSessionEpochs applies to the session means). A{i} = animal.
    T = {};  A = {};
    for sess = 1:numel(code)
        if sess > numel(allE1), break; end
        switch code(sess)
            case 1,    raw = allE1{sess};
            case 2,    raw = allE2{sess};
            otherwise, continue
        end
        if isempty(raw) || ~isfield(raw,'grid') || ~isfield(raw,'trialS'), continue; end
        if p > size(raw.grid,1) || w > size(raw.grid,2) || p > numel(raw.trialS), continue; end
        E = raw.grid{p,w};  St = raw.trialS{p};  Rt = raw.trialR{p};
        if isempty(E) || isempty(St), continue; end
        rngP = [nan nan];
        if isfield(raw,'traceRange') && p <= size(raw.traceRange,1), rngP = raw.traceRange(p,:); end
        [off, scl] = normMapSession(E, mode, rngP);
        T{end+1} = ([Rt(w,:); St] - off) ./ scl;   %#ok<AGROW>
        if sess <= numel(sessAnm) && ~isempty(sessAnm{sess})
            A{end+1} = char(sessAnm{sess});   %#ok<AGROW>
        else
            A{end+1} = sprintf('sess%02d', sess);   %#ok<AGROW>
        end
    end
end

function [off, scl] = normMapSession(E, mode, traceRange)
% the affine map of normSessionEpochs as (E - off) ./ scl, so it can
% be applied to single trials. Mirrors normSessionEpochs case for case.
    v = E(:);  off = 0;  scl = 1;
    if all(~isfinite(v)), off = NaN; scl = NaN; return; end
    switch lower(mode)
        case 'none'
        case 'traceminmax', off = traceRange(1);  scl = traceRange(2) - traceRange(1);
        case 'minmax4',     off = min(v);         scl = max(v) - min(v);
        case 'peak4',       off = 0;              scl = max(v);
        case 'resttopeak',  off = v(1);           scl = max(v(2:end)) - v(1);
        otherwise, error('Unknown cfg.sessionNorm: %s', mode);
    end
    bad = ~isfinite(scl) || ~isfinite(off) || scl == 0;
    if strcmpi(mode,'minmax4') || strcmpi(mode,'resttopeak') || strcmpi(mode,'traceminmax')
        bad = bad || scl < 0;
    end
    if bad, off = NaN; scl = NaN; end
end

function q = bhAdjustPA(p)
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

function out = ternaryPA(c, a, b)
    if c, out = a; else, out = b; end
end

function [anm, val] = perAnimalMean(names, vals)
% Average a group's session-level values within each animal, so the day 1 vs
% day 5 comparison is one number per animal and the pairing is well defined.
    if isempty(names), anm = {}; val = []; return; end
    names = names(:);
    vals  = vals(:);
    [anm, ~, ic] = unique(names, 'stable');
    val = nan(numel(anm),1);
    for i = 1:numel(anm)
        val(i) = mean(vals(ic == i), 'omitnan');
    end
    anm = anm(:)';
    val = val(:)';
end

%% LOADING + CACHE

function [obj, pp, tongueAll] = srLoadSession(loader, dateStr, spec, params, cfg, datapth)
% Load one session and keep only what the analysis
% reads: obj.time / trialdat / psth / bp / pth, params.cluid / trialid, and the
% tongue_length column of kin.dat.
% Cache: one file per (loader, date); an unreadable file is deleted and the
% session reloaded. Saves go to a
% temporary file that is renamed only once complete.
    key = cacheKey(loader, dateStr, params);
    useCache = cfg.useCache && ~isempty(cfg.cacheDir);
    fname = '';
    if useCache
        if ~exist(cfg.cacheDir, 'dir'), mkdir(cfg.cacheDir); end
        fname = fullfile(cfg.cacheDir, sprintf('%s_%s_%s.mat', spec.name, ...
            func2str(loader), matlab.lang.makeValidName(dateStr)));
        if exist(fname, 'file')
            C = [];
            try
                C = load(fname, 'S', 'key');
            catch ME
                logf('  [cache] unreadable %s (%s) -- deleting and reloading\n', fname, ME.message);
                try, delete(fname); catch, end
            end
            if isstruct(C) && isfield(C,'key') && isfield(C,'S') && strcmp(C.key, key)
                [obj, pp, tongueAll] = unpackS(C.S);
                logf('  [cache] hit  %s %s\n', func2str(loader), dateStr);
                return
            elseif isstruct(C)
                logf('  [cache] stale %s %s (params changed) -- reloading\n', func2str(loader), dateStr);
            end
        end
    end

    tLoad = tic;
    meta = slimMeta(loader, dateStr);
    params.probe = {meta.probe};
    params.cluid = {};
    [o, prm, kin] = slimToLegacy(meta, params);
    o = o(1);  prm = prm(1);
    kinix = find(strcmp(kin.featLeg, 'tongue_length'), 1);
    assert(~isempty(kinix), 'tongue_length not found in kin.featLeg for %s', dateStr);

    S.time      = o.time;
    S.trialdat  = single(o.trialdat);   % single on disk; back to double on load
    S.psth      = o.psth;
    S.bp        = o.bp;
    S.pth       = o.pth;
    S.cluid     = prm.cluid;
    S.trialid   = prm.trialid;
    S.tongueAll = kin.dat(:,:,kinix);
    logf('  [load]  %s %s in %.1f s\n', func2str(loader), dateStr, toc(tLoad));

    if useCache
        tmp = [fname(1:end-4) '_tmp.mat'];
        try
            save(tmp, 'S', 'key', '-v7.3');
            movefile(tmp, fname, 'f');
        catch ME
            logf('  [cache] could not save %s (%s) -- continuing without it\n', fname, ME.message);
            try, delete(tmp); catch, end
        end
    end
% Return the freshly loaded arrays rather than the single-precision cached copy.
    obj = struct('time', o.time, 'trialdat', o.trialdat, 'psth', o.psth, 'bp', o.bp, 'pth', o.pth);
    pp  = struct('cluid', {prm.cluid}, 'trialid', {prm.trialid});
    tongueAll = S.tongueAll;
end

function [obj, pp, tongueAll] = unpackS(S)
    obj = struct('time', S.time, 'trialdat', double(S.trialdat), 'psth', S.psth, ...
                 'bp', S.bp, 'pth', S.pth);
    pp  = struct('cluid', {S.cluid}, 'trialid', {S.trialid});
    tongueAll = S.tongueAll;
end

function key = cacheKey(loader, dateStr, params)
% Everything that can change what gets loaded, and nothing that cannot.
    k.version = '1';
    k.loader  = func2str(loader);
    k.date    = dateStr;
    k.align   = params.alignEvent;
    k.behav   = params.behav_only;
    k.warp    = params.timeWarp;
    k.nLicks  = params.nLicks;
    k.lowFR   = params.lowFR;
    k.quality = params.quality;
    k.cond    = params.condition;
    k.tmin    = params.tmin;
    k.tmax    = params.tmax;
    k.dt      = params.dt;
    k.smooth  = params.smooth;
    k.traj    = params.traj_features;
    k.featVar = params.feat_varToExplain;
    k.NVar    = params.N_varToExplain;
    k.adv     = params.advance_movement;
    k.fcut    = params.fcut;
    k.condN   = params.cond;
    k.method  = params.method;
    k.fa      = params.fa;
    k.bctype  = params.bctype;
    k.source  = 'slimExport';   % marks entries built from the exported data
    key = jsonencode(k);
end

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

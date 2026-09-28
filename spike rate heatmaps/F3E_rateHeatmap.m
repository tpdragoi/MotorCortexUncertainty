%% F3E_rateHeatmap.m
%  Per-session spike-rate heatmaps (Fig. 3E), Double Reward Task: single-trial rate averaged
%  over each probe's units, trials sorted by the time of the last port contact in the bout.
%  Reads Data\<task>\<ANM>_<DATE>_obj.mat and _kin.mat (Data = dataRoot in setPaths.m)
%  through shared\slimMeta and shared\slimToLegacy.
%  Settings: aligned to goCue, dt 1/200 s, smooth 25, good + excellent units, lowFR 0.01 Hz,
%  window -1.5 to 4 s.

clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.

sz = 14;

% Locate the repository root and read data locations from setPaths.m.
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % this file sits one folder below the repository root
if isempty(repoRoot) || ~exist(fullfile(repoRoot, 'setPaths.m'), 'file'), repoRoot = pwd; end
assert(exist(fullfile(repoRoot, 'setPaths.m'), 'file') == 2, ...
    'Cannot find setPaths.m. Run this script from its file, or cd to the repository root first.');
addpath(repoRoot);
cfgPaths = setPaths();   % data locations are set once, in setPaths.m

%% PARAMETERS
params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping applies to neural data only
params.behav_only = 0;
params.timeWarp            = 0;   % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20;   % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.01;   % minimum mean firing rate, Hz

params.condition(1) = {'hit==1 | hit==0' };   % 1: all trials
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};   % 2: hits, trial type 1, reward on lick 1
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};   % 3: hits, trial type 2, reward on lick 1
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};   % 4: hits, trial type 3, reward on lick 1
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 6'};   % 5: hits, trial type 1, reward on lick 6
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 6'};   % 6: hits, trial type 2, reward on lick 6
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 6'};   % 7: hits, trial type 3, reward on lick 6
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};   % 8: hits, reward on lick 1
params.condition(end+1) = {'hit==1 & rewardedLick == 6'};   % 9: hits, reward on lick 6
params.condition(end+1) = {'hit==1' };   % 10: all hits

params.tmin = -1.5;
params.tmax = 4;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 25;

% cluster qualities to use
params.quality = {'good','excellent',' good','good '};   % good + excellent units

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};
params.feat_varToExplain = 80;   % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80;   % keep num dims that explain this much variance in neural data
params.advance_movement = 0;

% Params for finding kinematic modes
params.fcut = 10;   % smoothing cutoff frequency
params.cond = 5;   % which conditions to use to find mode
params.method = 'xcorr';   % 'xcorr' | 'regress'
params.fa = false;   % if true, reduces neural dimensions to 10 with factor analysis
params.bctype = 'reflect';   % 'reflect' | 'zeropad' | 'none'

%% SPECIFY DATA TO LOAD

datapth = '';   % not used

% one empty placeholder per session slot; the ones a task uses are
% filled in below and the rest drop out of the all_meta concatenation
[meta, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11, ...
    meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19, meta20, meta21, ...
    meta22, meta23, meta24, meta25, meta26, meta27, meta28, meta29, meta30, meta31, ...
    meta32, meta33, meta34, meta35, meta36, meta37, meta38, meta39, meta40] = deal([]);
date = '2025-07-17';
meta39 = slimMeta('loadTD25_neural', date);

y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12 ...
    ;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24 ...
    ;meta25;meta26;meta27;meta28;meta29;meta30;meta31;meta32;meta33;meta34;meta35;meta36...
    ;meta37;meta38;meta39;meta40];

y1_all = [];
y2_all = [];
y3_all = [];
y4_all = [];

allmoveP1 = {};
allmoveP4 = {};

for sessnum = 1:length(all_meta)

clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin

meta = all_meta(sessnum,1);

params.probe = {meta.probe};

% LOAD DATA
clear obj

params.cluid = {};
[obj, params, kin] = slimToLegacy(meta, params);

trialSet = [1:obj.bp.Ntrials]';

fprintf('Session %d\n', sessnum);

conds2use = [1];   % index into params.condition
kinfeat = 'tongue_length';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
sessix = 1;

psthForProj = [];
for c = conds2use
    condtrix = trialSet;   % trials of this condition

    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    condtrix(condtrix > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 298) = [];
     elseif strcmp(obj.pth.dt,'2023-12-16') && strcmp(obj.pth.anm,'TD4f')
        condtrix(condtrix > 203) = [];
     elseif strcmp(obj.pth.dt,'2025-07-17') && strcmp(obj.pth.anm,'TD25d')
        condtrix(condtrix > 199) = [];
    end

end

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
Length = kin.dat(:,condtrix,kinix);

%% Calculate Last Lick

hitTrials = params.trialid{1,1};

LastL = [];
for i = 1:length(hitTrials)

    tr = hitTrials(i);
    lickL = obj.bp.ev.lickL{i};
    if isempty(lickL)
        LastL = [LastL 0];
    elseif ~isempty(lickL)
        lickL = lickL(lickL > obj.bp.ev.goCue(i));
        if ~isempty(lickL)
            LastL = [LastL lickL(end) - obj.bp.ev.goCue(i)];
        else
            LastL = [LastL 0];
        end
    end

end

all11 = [1:obj.bp.Ntrials];
hit11 = all11((obj.bp.hit == 1));
r111 = all11((obj.bp.rewardedLick == 1));
r444 = all11((obj.bp.rewardedLick == 6));

allr11 = intersect(hit11,r111)';
allr44 = intersect(hit11,r444)';

P8 = allr11;
P9 = allr44;

if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
P9(P9 > 278) = [];
elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
P8(P8 > 313) = [];
P9(P9 > 313) = [];
elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
P8(P8 > 298) = [];
P9(P9 > 298) = [];
elseif strcmp(obj.pth.dt,'2023-12-16') && strcmp(obj.pth.anm,'TD4f')
P8(P8 > 203) = [];
P9(P9 > 203) = [];
elseif strcmp(obj.pth.dt,'2025-07-17') && strcmp(obj.pth.anm,'TD25d')
P8(P8 > 199) = [];
P9(P9 > 199) = [];
end

%% TONGUE

task = 14;

conds2use = [1];   % index into params.condition
kinfeat = 'tongue_length';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1

sessix = 1;

Ncells = size(obj.psth, 2);

condpsth = obj.trialdat;

    if strcmp(obj.pth.dt,'2025-06-30') && strcmp(obj.pth.anm,'TD24d') || strcmp(obj.pth.dt,'2025-07-01') && strcmp(obj.pth.anm,'TD24d')

        Ncells = size(obj.psth, 2);
        clu_m1TJ = 1:size(params.cluid,1);
        clu_ALM = clu_m1TJ;

        reg = clu_m1TJ;
        reg1 = clu_ALM;
        neurall{1} = [reg];
        neurall{2} = [reg1];

    else

neurall = {};
Ncells = size(obj.psth, 2);
clu_m1TJ = 1:size(params.cluid{1, 1},1);
clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;

reg = clu_m1TJ;
reg1 = clu_ALM;
neurall{1} = [reg];
neurall{2} = [reg1];
    end

allreg = {};

allreg{1} = [reg];
allreg{2} = [reg1];

figure;

for kk = 1:2

numClu = allreg{kk};
% skip a probe with no units in this session
if isempty(numClu)
    logf('  %s %s: no units on probe %d -- panel skipped\n', obj.pth.anm, obj.pth.dt, kk);
    continue
end

% single-trial rates for this condition
psthForProj = [];
for c = conds2use
    condtrix = trialSet;   % trials of this condition

    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    condtrix(condtrix > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 298) = [];
     elseif strcmp(obj.pth.dt,'2023-12-16') && strcmp(obj.pth.anm,'TD4f')
        condtrix(condtrix > 203) = [];
     elseif strcmp(obj.pth.dt,'2025-07-17') && strcmp(obj.pth.anm,'TD25d')
        condtrix(condtrix > 199) = [];
    end
condpsth = obj(sessix).trialdat(:,:,condtrix);   % Take the single trial PSTHs for these trials

end

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));

% mean rate across this probe's units: [time x trials]

condpsth = condpsth(:,numClu,:);

tnsorp1 = squeeze(mean(condpsth, 2));

Kinematics1 = kin.dat(:,condtrix,kinix);

tnsorp_r1 = tnsorp1(:,P8);
tnsorp_r4 = tnsorp1(:,P9);

% last lick time: last sample where tongue length differs from its mode

clear mode

% R1 trials
Kinematics1 = kin.dat(:,condtrix,kinix);
Kinematics1 = Kinematics1(:, P8);

tol = 1e-9;
Kinematics1 = Kinematics1 - mode(Kinematics1);
last_nonzero_index = max((1:size(Kinematics1,1)).' .* (abs(Kinematics1) > tol), [], 1);
[~, idx1] = sort(last_nonzero_index);

last_nonzero_index = last_nonzero_index(idx1);
tnsorp_r1          = tnsorp_r1(:, idx1);
sorted_P8          = P8(idx1);   % keep P8 sorted alongside

lastlickr1 = zeros(1, size(tnsorp_r1, 2));
rewardr1   = nan(1,   size(tnsorp_r1, 2));

for i = 1:size(tnsorp_r1, 2)
    temp = sorted_P8(i);   % trial ID in sorted order
    if last_nonzero_index(i) ~= 0
        lastlickr1(i) = obj.time(last_nonzero_index(i));
        if task == 16
            if ~isnan(obj.bp.ev.lickL{temp,1})
                liks  = obj.bp.ev.lickL{temp,1} > obj.bp.ev.goCue(temp);
                licks = obj.bp.ev.lickL{temp,1}(liks);
                if ~isempty(licks) && ~any(isnan(licks))
                    rewardr1(i) = licks(1) - obj.bp.ev.goCue(temp);
                end
            end
        else
            if ~isnan(obj.bp.ev.lickL{temp,1})
                liks  = obj.bp.ev.lickL{temp,1} > obj.bp.ev.goCue(temp);
                licks = obj.bp.ev.lickL{temp,1}(liks);
                if ~isempty(licks) && ~any(isnan(licks))
                    rewardr1(i) = licks(1) - obj.bp.ev.goCue(temp);
                end
            end
        end
    else
        lastlickr1(i) = obj.time(600);
    end
end

% Sort small to large and remove trials with last lick < 0.1
[lastlickr1_sorted, sort_idx] = sort(lastlickr1, 'ascend');
keep      = lastlickr1_sorted >= 0.1;
lastlickr1 = lastlickr1_sorted(keep);
tnsorp_r1  = tnsorp_r1(:, sort_idx(keep));
sorted_P8  = sorted_P8(sort_idx(keep));
rewardr1   = rewardr1(sort_idx(keep));

% R6 trials (P9)
Kinematics1 = kin.dat(:,condtrix,kinix);
Kinematics1 = Kinematics1(:, P9);

tol = 1e-9;
Kinematics1 = Kinematics1 - mode(Kinematics1);
last_nonzero_index = max((1:size(Kinematics1,1)).' .* (abs(Kinematics1) > tol), [], 1);
[~, idx1] = sort(last_nonzero_index);

last_nonzero_index = last_nonzero_index(idx1);
tnsorp_r4          = tnsorp_r4(:, idx1);
sorted_P9          = P9(idx1);   % keep P9 sorted alongside

lastlickr4 = zeros(1, size(tnsorp_r4, 2));
rewardr4   = nan(1,   size(tnsorp_r4, 2));
rewardr41  = nan(1,   size(tnsorp_r4, 2));

for i = 1:size(tnsorp_r4, 2)
    temp = sorted_P9(i);   % trial ID in sorted order
    if last_nonzero_index(i) ~= 0
        lastlickr4(i) = obj.time(last_nonzero_index(i));
        if task == 14
            if ~isnan(obj.bp.ev.lickL{temp,1})
                liks  = obj.bp.ev.lickL{temp,1} > obj.bp.ev.goCue(temp);
                licks = obj.bp.ev.lickL{temp,1}(liks);
            end
            rewardr4(i) = obj.bp.ev.reward(temp) - obj.bp.ev.goCue(temp);
        else
            if ~isnan(obj.bp.ev.lickL{temp,1})
                liks  = obj.bp.ev.lickL{temp,1} > obj.bp.ev.goCue(temp);
                licks = obj.bp.ev.lickL{temp,1}(liks);
            end
            rewardr4(i) = obj.bp.ev.reward(temp) - obj.bp.ev.goCue(temp);
            if ~isnan(obj.bp.ev.lickL{temp,1})
                rewardr41(i) = obj.bp.ev.lickL{temp,1}(1) - obj.bp.ev.goCue(temp);
            end
        end
    else
        lastlickr4(i) = obj.time(600);
    end
end

% Sort small to large and remove trials with last lick < 0.1
[lastlickr4_sorted, sort_idx] = sort(lastlickr4, 'ascend');
keep      = lastlickr4_sorted >= 0.1;
lastlickr4 = lastlickr4_sorted(keep);
tnsorp_r4  = tnsorp_r4(:, sort_idx(keep));
sorted_P9  = sorted_P9(sort_idx(keep));
rewardr4   = rewardr4(sort_idx(keep));
rewardr41  = rewardr41(sort_idx(keep));

tnsorp_r1 = abs(tnsorp_r1);
tnsorp_r4 = abs(tnsorp_r4);

numColors = 100;
blue  = [0/255, 0/255, 153/255];
white = [1, 1, 1];
yy    = 1.2;
t     = linspace(0,1,numColors).^yy;
blueToWhite = [ ...
    blue(1) + (white(1)-blue(1)) * t', ...
    blue(2) + (white(2)-blue(2)) * t', ...
    blue(3) + (white(3)-blue(3)) * t'  ...
    ];

sz     = 10;
nR1    = size(tnsorp_r1, 2);
nR4    = size(tnsorp_r4, 2);
orange = [1 0.5 0.1];
pink   = [1 0 0.68];
allData = [tnsorp_r1(:); tnsorp_r4(:)];
% ignore NaNs and guard against a flat or empty range, which caxis rejects
finiteData = allData(isfinite(allData));
if isempty(finiteData), finiteData = 0; end
p_all   = prctile(finiteData, [2 98]);
rang    = [p_all(1) p_all(2)];
if ~all(isfinite(rang)) || rang(2) <= rang(1)
    rang = [min(finiteData) max(finiteData)];
end
if ~all(isfinite(rang)) || rang(2) <= rang(1)
    rang = [0 1];
end
cmp = blueToWhite;

% 1) R1 panel
subplot(2,2,kk)
imagesc(obj.time, 1:nR1, tnsorp_r1');
colormap(cmp); caxis(rang);
cb1 = colorbar;
cb1.Ticks      = [rang(1), rang(2)];
cb1.TickLabels = {sprintf('%.2f', rang(1)), sprintf('%.2f', rang(2))};
cb1.TickLength = 0;
axis tight;
xlabel(['time from ' num2str(params.alignEvent)])
ylabel('Trial # (R1)')
title(sprintf('Date %s  Animal %s  Probe %d (R1)', obj.pth.dt, obj.pth.anm, kk))
xlim([-0.3 3.2])
set(gca, 'FontSize', sz)
hold on
for i = 1:nR1
    line([lastlickr1(i) lastlickr1(i)], [i-0.5 i+0.5], 'Color', orange, 'LineWidth', 3);
end
hold off

% 2) R4 panel
subplot(2,2,kk+2)
imagesc(obj.time, 1:nR4, tnsorp_r4');
colormap(cmp); caxis(rang);
cb2 = colorbar;
cb2.Ticks      = [rang(1), rang(2)];
cb2.TickLabels = {sprintf('%.2f', rang(1)), sprintf('%.2f', rang(2))};
cb2.TickLength = 0;
axis tight;
xlabel(['time from ' num2str(params.alignEvent)])
ylabel('Trial # (R4)')
title(sprintf('Date %s  Animal %s  Probe %d (R4)', obj.pth.dt, obj.pth.anm, kk))
xlim([-0.3 3.2])
set(gca, 'FontSize', sz)
hold on
for i = 1:nR4
    line([lastlickr4(i) lastlickr4(i)], [i-0.5 i+0.5], 'Color', orange, 'LineWidth', 3);
end
hold off

set(gcf, 'Position', [50 100 600 800]);

end

end

a = 3;

a = 2;

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

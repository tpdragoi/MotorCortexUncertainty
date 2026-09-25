%% F1E_rateHeatmap.m
%  Per-session spike-rate heatmaps, Simple Reward Task.
%  Trials are sorted by the time of the last port contact in the bout.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/200
%    params.smooth      35
%    params.quality     {'good'}
%    params.lowFR       0.01
%    params.window      -1.5 to 5 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


sz = 14;

% addpath 'C:\Users\vdragoi\Desktop\uninstructedMovements_v2-main\ObjVis\warp'

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
params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping only operates on neural data for now.
params.behav_only = 0;
params.timeWarp            = 0;   % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20;   % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.01;   % minimum mean firing rate, Hz

params.condition(1) = {'hit==1 | hit==0' };   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 4'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 4'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 4'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1' };   % left to right         % right hits, no stim, aw off

% params.condition(1) = {'hit==1 | hit==0' };    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 6'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 6'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 6'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & rewardedLick == 6'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1' };    % left to right         % right hits, no stim, aw off

params.tmin = -1.5;
params.tmax = 5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 35;

% cluster qualities to use
params.quality = {'good'};   % good units only (findClusters trims blanks and ignores case)

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

%% SPECIFY DATA TO LOAD

datapth = '';   % raw data folder not used (was: datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';)

% one empty placeholder per session slot; the ones a task uses are
% filled in below and the rest drop out of the all_meta concatenation
[meta, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11, ...
    meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19, meta20, meta21, ...
    meta22, meta23, meta24, meta25, meta26, meta27, meta28, meta29, meta30, meta31, ...
    meta32, meta33, meta34, meta35, meta36] = deal([]);
date = '2024-07-11';
meta3 = slimMeta('loadTD10s_neur', date);
date = '2024-07-09';
meta11 = slimMeta('loadTD9s_neur709', date);

y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12 ...
    ;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24 ...
    ;meta25;meta26;meta27;meta28;meta29;meta30;meta31;meta32;meta33;meta34;meta35;meta36];
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

for sessix = 1:numel(meta)
end

% Get kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix),{' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
end

%% Define trial sets (P8 = rewardedLick==1 hits, P9 = rewardedLick==4 hits)

all11  = 1:obj.bp.Ntrials;
hit11  = all11(obj.bp.hit == 1);
r111   = all11(obj.bp.rewardedLick == 1);
r444   = all11(obj.bp.rewardedLick == 4);

P8 = intersect(hit11, r111)';
P9 = intersect(hit11, r444)';

if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    P9(P9 > 278) = [];
elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
    P8(P8 > 313) = []; P9(P9 > 313) = [];
elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
    P8(P8 > 298) = []; P9(P9 > 298) = [];
end

%% Define brain regions

task   = 14;
kinfeat = 'tongue_length';
sessix  = 1;
Ncells  = size(obj.psth, 2);

if strcmp(obj.pth.dt,'2024-07-13') && strcmp(obj.pth.anm,'TD10si') || strcmp(obj.pth.dt,'2024-07-09') && strcmp(obj.pth.anm,'TD9si') ...
        || strcmp(obj.pth.dt,'2025-02-22') && strcmp(obj.pth.anm,'TDl3')  || strcmp(obj.pth.dt,'2025-02-21') && strcmp(obj.pth.anm,'TDl2') ...
        || strcmp(obj.pth.dt,'2025-02-19') && strcmp(obj.pth.anm,'TDl2') || strcmp(obj.pth.dt,'2025-04-21') && strcmp(obj.pth.anm,'TD20d') ...
        || strcmp(obj.pth.dt,'2025-08-02') && strcmp(obj.pth.anm,'TD26d') || strcmp(obj.pth.dt,'2025-08-03') && strcmp(obj.pth.anm,'TD26d')
    clu_m1TJ = 1:size(params.cluid,1);
    reg  = clu_m1TJ;
    reg1 = clu_m1TJ;
else
    clu_m1TJ = 1:size(params.cluid{1,1},1);
    clu_ALM  = size(params.cluid{1,1},1)+1:Ncells;
    reg  = clu_m1TJ;
    reg1 = clu_ALM;
end

allreg{1} = reg;
allreg{2} = reg1;

%% Main loop over probes

figure;

for kk = 1:2

    numClu = allreg{kk};
% a probe the export did not keep (no group map used it) has no units
% here, so there is nothing to plot for that panel.
    if isempty(numClu)
        logf('  [clean] %s %s: no units on probe %d -- panel skipped\n', obj.pth.anm, obj.pth.dt, kk);
        continue
    end

% Get condtrix (trial indices for this condition)
    condtrix = params(sessix).trialid{1};
    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
        condtrix(condtrix > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 298) = [];
    end

% Get mean spike rate across neurons: [time x trials]
    condpsth = obj(sessix).trialdat(:, numClu, condtrix);
    tnsorp1  = squeeze(mean(condpsth, 2));   % [time x trials]

% kinematic feature
    kinix = find(strcmp(kin(sessix).featLeg, kinfeat));

% Map P8 and P9 (absolute trial numbers) to local indices within condtrix
    [~, P8_local] = ismember(P8, condtrix);
    P8_local = P8_local(P8_local > 0);
    P8_used  = P8(P8_local > 0);

    [~, P9_local] = ismember(P9, condtrix);
    P9_local = P9_local(P9_local > 0);
    P9_used  = P9(P9_local > 0);

% Slice neural data using local indices
    tnsorp_r1 = tnsorp1(:, P8_local);   % [time x nR1]
    tnsorp_r4 = tnsorp1(:, P9_local);   % [time x nR4]

% Slice kinematics using same local indices
    Kin_all = kin(sessix).dat(:, condtrix, kinix);   % [time x nCondTrials]
    Kin_r1  = Kin_all(:, P8_local);   % [time x nR1]
    Kin_r4  = Kin_all(:, P9_local);   % [time x nR4]

% Sort R1 trials by last lick time

    tol = 1e-9;
    Kin_r1_centered = Kin_r1 - mode(Kin_r1(:));
    last_nonzero_r1 = max((1:size(Kin_r1_centered,1)).' .* (abs(Kin_r1_centered) > tol), [], 1);
    [~, idx1] = sort(last_nonzero_r1);

    last_nonzero_r1    = last_nonzero_r1(idx1);
    tnsorp_r1          = tnsorp_r1(:, idx1);
    sorted_trialids_r1 = P8_used(idx1);

% Compute lastlickr1 and rewardr1
    lastlickr1 = zeros(1, size(tnsorp_r1, 2));
    rewardr1   = nan(1,   size(tnsorp_r1, 2));

    for i = 1:size(tnsorp_r1, 2)
        temp = sorted_trialids_r1(i);
        if last_nonzero_r1(i) ~= 0
            lastlickr1(i) = obj.time(last_nonzero_r1(i));
            if ~isnan(obj.bp.ev.lickL{temp,1})
                liks  = obj.bp.ev.lickL{temp,1} > obj.bp.ev.goCue(temp);
                licks = obj.bp.ev.lickL{temp,1}(liks);
                if ~isempty(licks) && ~any(isnan(licks))
                    rewardr1(i) = licks(1) - obj.bp.ev.goCue(temp);
                end
            end
        else
            lastlickr1(i) = obj.time(600);
        end
    end

% Sort lastlickr1 small to large and remove trials with last lick < 0.1
    [lastlickr1_sorted, sort_idx] = sort(lastlickr1, 'ascend');
    keep = lastlickr1_sorted >= 0.1;

    lastlickr1         = lastlickr1_sorted(keep);
    tnsorp_r1          = tnsorp_r1(:, sort_idx(keep));
    sorted_trialids_r1 = sorted_trialids_r1(sort_idx(keep));
    rewardr1           = rewardr1(sort_idx(keep));

% Sort R4 trials by last lick time

    Kin_r4_centered = Kin_r4 - mode(Kin_r4(:));
    last_nonzero_r4 = max((1:size(Kin_r4_centered,1)).' .* (abs(Kin_r4_centered) > tol), [], 1);
    [~, idx1] = sort(last_nonzero_r4);

    last_nonzero_r4    = last_nonzero_r4(idx1);
    tnsorp_r4          = tnsorp_r4(:, idx1);
    sorted_trialids_r4 = P9_used(idx1);

% Compute lastlickr4 and rewardr4
    lastlickr4 = zeros(1, size(tnsorp_r4, 2));
    rewardr4   = nan(1,   size(tnsorp_r4, 2));

    for i = 1:size(tnsorp_r4, 2)
        temp = sorted_trialids_r4(i);
        if last_nonzero_r4(i) ~= 0
            lastlickr4(i) = obj.time(last_nonzero_r4(i));
            rewardr4(i)   = obj.bp.ev.reward(temp) - obj.bp.ev.goCue(temp);
        else
            lastlickr4(i) = obj.time(600);
        end
    end

% Sort lastlickr4 small to large and remove trials with last lick < 0.1
    [lastlickr4_sorted, sort_idx] = sort(lastlickr4, 'ascend');
    keep = lastlickr4_sorted >= 0.1;

    lastlickr4         = lastlickr4_sorted(keep);
    tnsorp_r4          = tnsorp_r4(:, sort_idx(keep));
    sorted_trialids_r4 = sorted_trialids_r4(sort_idx(keep));
    rewardr4           = rewardr4(sort_idx(keep));

% Colormaps

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

    sz = 10;
    nR1 = size(tnsorp_r1, 2);
    nR4 = size(tnsorp_r4, 2);

    tnsorp_r1 = tnsorp_r1(1:size(tnsorp_r1, 1), :);
    tnsorp_r4 = tnsorp_r4(1:size(tnsorp_r1, 1), :);
    allData = [tnsorp_r1(:); tnsorp_r4(:)];
% ignore NaNs and guard against a flat or empty range, which caxis rejects
    finiteData = allData(isfinite(allData));
    if isempty(finiteData), finiteData = 0; end
    p_all   = prctile(finiteData, [1 99]);
    rang    = [p_all(1) p_all(2)];
    if ~all(isfinite(rang)) || rang(2) <= rang(1)
        rang = [min(finiteData) max(finiteData)];
    end
    if ~all(isfinite(rang)) || rang(2) <= rang(1)
        rang = [0 1];
    end
    combined_data = [tnsorp_r1'; tnsorp_r4'];

    subplot(1,2,kk)
    imagesc(obj.time(1:1300), 1:(nR1+nR4), combined_data);
    colormap(blueToWhite)
    caxis(rang);
    cb = colorbar;
    cb.Ticks      = [rang(1), rang(2)];
    cb.TickLabels = {sprintf('%.2f', rang(1)), sprintf('%.2f', rang(2))};
    cb.TickLength = 0;
    grid off; axis tight;
    box off;
    xlabel(['time from ', num2str(params.alignEvent)]);
    ylabel('Trial #');
    title(['Date ', num2str(obj.pth.dt), ' Animal ', num2str(obj.pth.anm), ' Probe ', num2str(kk)]);
    xlim([-0.3 3.75])
    set(gca, 'FontSize', sz)

    orange = [1 0.5 0.1];
    for i = 1:nR1
        line([lastlickr1(i) lastlickr1(i)], [i-0.5 i+0.5], 'Color', orange, 'LineWidth', 3);
    end
    for i = 1:nR4
        y = nR1 + i;
        line([lastlickr4(i) lastlickr4(i)], [y-0.5 y+0.5], 'Color', orange, 'LineWidth', 3);
    end

    set(gcf, 'Position', [50 100 800 400]);

end   % kk loop

end   % sessnum loop

a = 3;


% colorbar; grid off; axis tight;
%         'LineWidth', 3, 'Color', orange);
% colorbar; grid off; axis tight;
%         'LineWidth', 3, 'Color', orange);

a = 2;

function logf(varargin)
% Progress and diagnostic messages, silenced by default.
% Set verbose = true to print them.
verbose = false;
if verbose
    fprintf(varargin{:});
end
end

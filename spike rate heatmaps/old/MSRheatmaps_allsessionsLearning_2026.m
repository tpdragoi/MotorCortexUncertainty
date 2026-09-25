% Finding "Kinematic Modes"
clear,clc

sz = 14;

d = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main';

addpath(genpath(fullfile(d,'utils')))
addpath(genpath(fullfile(d,'DataLoadingScripts')))
addpath(genpath(fullfile(d,'funcs')))
rmpath(genpath(fullfile(d,'fig1')));

addpath 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\base code\functions_td'
addpath 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\ObjVis\warp'



%% PARAMETERS
params.alignEvent          = 'goCue'; % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping only operates on neural data for now.
params.behav_only = 0;
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.1; % remove clusters with firing rates across all trials less than this val

params.condition(1) = {'hit==1 | hit==0' };    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1' };    % left to right         % right hits, no stim, aw off

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
params.tmax = 4;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 25;

% cluster qualities to use
% params.quality = {'ok','good','mua','great'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'good','excellent',' good', 'good ','ok'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'good'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};
params.feat_varToExplain = 80;  % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80;     % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;

% Params for finding kinematic modes
params.fcut = 10;          % smoothing cutoff frequency
params.cond = 5;         % which conditions to use to find mode
params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
params.fa = false;         % if true, reduces neural dimensions to 10 with factor analysis
params.bctype = 'reflect'; % options are : reflect  zeropad  none

%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';


meta = [];
meta1 = [];
meta2 = [];
meta3 = [];
meta4 = [];
meta5 = [];
meta6 = [];
meta7 = [];
meta8 = [];
meta9 = [];
meta10 = [];
meta11 = [];
meta12 = [];
meta13 = [];
meta14 = [];
meta15 = [];
meta16 = [];
meta17 = [];
meta18 = [];
meta19 = [];
meta20 = [];
meta21 = [];
meta22 = [];
meta23 = [];
meta24 = [];
meta25 = [];
meta26 = [];
meta27 = [];
meta28 = [];
meta29 = [];
meta30 = [];
meta31 = [];
meta32 = [];
meta33 = [];
meta34 = [];
meta35 = [];
meta36 = [];

% meta1 = loadTD_inactivation1(meta1,datapth);
% meta2 = loadTD_inactivation2(meta2,datapth);
% meta3 = loadTD_inactivation3(meta3,datapth);
% meta4 = loadTD_inactivation4(meta4,datapth);
% meta5 = loadTD_inactivation5(meta5,datapth);


date = '2025-08-21';
meta1 = loadTD5l_neur(meta1,datapth,date);
date = '2025-08-22';
meta2 = loadTD5l_neur(meta2,datapth,date);
date = '2025-08-23';
meta3 = loadTD5l_neur(meta3,datapth,date);
date = '2025-08-24';
meta4 = loadTD5l_neur(meta4,datapth,date);
date = '2025-08-25';
meta5 = loadTD5l_neur(meta5,datapth,date);

date = '2025-06-03';
meta6 = loadTD4l_neur(meta6,datapth,date);
date = '2025-06-04';
meta7 = loadTD4l_neur(meta7,datapth,date);
date = '2025-06-05';
meta8 = loadTD4l_neur(meta8,datapth,date);
date = '2025-06-06';
meta9 = loadTD4l_neur(meta9,datapth,date);
date = '2025-06-07';
meta10 = loadTD4l_neur(meta10,datapth,date);

% TD15d R1
date = '2025-02-03';
meta11 = loadTD2l_neur(meta11,datapth,date);
date = '2025-02-04';
meta12 = loadTD2l_neur219(meta12,datapth,date);
date = '2025-02-05';
meta13 = loadTD2l_neur(meta13,datapth,date);
date = '2025-02-06';
meta14 = loadTD2l_neur(meta14,datapth,date);
date = '2025-02-07';
meta15 = loadTD2l_neur(meta15,datapth,date);
% TD8d R1

date = '2025-02-03';
meta16 = loadTD3l_neur(meta16,datapth,date);
date = '2025-02-04';
meta17 = loadTD3l_neur(meta17,datapth,date);
date = '2025-02-05';
meta18 = loadTD3l_neur(meta18,datapth,date);
date = '2025-02-06';
meta19 = loadTD3l_neur(meta19,datapth,date);
date = '2025-02-07';
meta20 = loadTD3l_neur(meta20,datapth,date);% 






y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

% all_meta = [meta1;meta2;meta3;meta4;meta5];
% all_meta = [meta2;meta1;meta3;meta1;meta1];
% all_meta = [meta1];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12 ...
    ;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24 ...
    ;meta25;meta26;meta27;meta28;meta29;meta30;meta31;meta32;meta33;meta34;meta35;meta36];
y1_all = [];
y2_all = [];
y3_all = [];
y4_all = [];

allmoveP1 = {};
allmoveP4 = {};



%% ============================================================
%  ONE FIGURE: 5 rows (sessions) x 4 columns (animals)
% =============================================================
nAnimals       = 4;
nSessPerAnimal = 5;
nRows          = nSessPerAnimal;
nCols          = nAnimals;

probeToUse = [2, 1, 2, 1, 1, ...   % Animal 1
              2, 1, 1, 1, 1, ...   % Animal 2
              1, 2, 1, 2, 1, ...   % Animal 3
              1, 2, 1, 1, 1];      % Animal 4

fig = figure('Position', [50 50 900 900]);

% ---- Layout parameters — larger gaps = smaller plots ----
lm = 0.11;     % left margin (wider to fit ylabels on every column)
rm = 0.04;
bm = 0.08;
tm = 0.05;
hg = 0.055;    % increased horizontal gap → narrower plots
vg = 0.042;    % increased vertical gap  → shorter plots

aw = (1 - lm - rm - (nCols-1)*hg) / nCols;
ah = (1 - bm - tm - (nRows-1)*vg) / nRows;

% ---- Pre-create ALL axes ----
axArray = gobjects(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        x0 = lm + (c-1)*(aw + hg);
        y0 = bm + (nRows - r)*(ah + vg);
        axArray(r,c) = axes(fig, 'Position', [x0 y0 aw ah]);
        axis(axArray(r,c), 'off');
    end
end

%% ============================================================
%  MAIN SESSION LOOP
% =============================================================
for sessnum = 1:length(all_meta)

    clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 ...
          Length angle obj aa aaa idxHit kin me ...
          lastlickr1 lastlickr4 rewardr4 ...
          tnsorp_r1 tnsorp_r4 P8 P9

    kk               = probeToUse(sessnum);
    animalIdx        = ceil(sessnum / nSessPerAnimal);
    sessWithinAnimal = mod(sessnum-1, nSessPerAnimal) + 1;

    ax = axArray(sessWithinAnimal, animalIdx);

    % ---- Load metadata & session data ----
    meta         = all_meta(sessnum, 1);
    params.probe = {meta.probe};
    clear obj
    params.cluid = {};
    [obj, params] = loadSessionData(meta, params, params.behav_only);

    trialSet = (1:obj.bp.Ntrials)';

    for sessix = 1:numel(meta)
        me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    end

    nSessions = numel(meta);
    for sessix = 1:numel(meta)
        fprintf('---- Kinematics: session %d / %d ----\n', sessix, nSessions);
        kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
    end

    conds2use = [1];
    kinfeat   = 'tongue_length';
    sessix    = 1;
    task      = 14;

    % ---- condtrix with trial limits ----
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

    kinix = find(strcmp(kin(sessix).featLeg, kinfeat));

    % ---- Trial type indices ----
    all11 = 1:obj.bp.Ntrials;
    hit11 = all11(obj.bp.hit == 1);
    r111  = all11(obj.bp.rewardedLick == 1);
    r444  = all11(obj.bp.rewardedLick == 4);

    P8 = intersect(hit11, r111)';
    P9 = intersect(hit11, r444)';

    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
        P8(P8 > 278) = []; P9(P9 > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        P8(P8 > 313) = []; P9(P9 > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        P8(P8 > 298) = []; P9(P9 > 298) = [];
    elseif strcmp(obj.pth.dt,'2025-06-05') && strcmp(obj.pth.anm,'TDl4')
        P8(P8 > 123) = []; P9(P9 > 123) = [];
    end

    % ---- Region/neuron assignment ----
    specialCase = (strcmp(obj.pth.dt,'2024-07-13') && strcmp(obj.pth.anm,'TD10si')) || ...
                  (strcmp(obj.pth.dt,'2024-07-09') && strcmp(obj.pth.anm,'TD9si'))  || ...
                  (strcmp(obj.pth.dt,'2025-02-22') && strcmp(obj.pth.anm,'TDl3'))   || ...
                  (strcmp(obj.pth.dt,'2025-02-21') && strcmp(obj.pth.anm,'TDl2'))   || ...
                  (strcmp(obj.pth.dt,'2025-02-19') && strcmp(obj.pth.anm,'TDl2'));

    Ncells = size(obj.psth, 2);
    if specialCase
        clu_m1TJ  = 1:size(params.cluid, 1);
        allreg{1} = clu_m1TJ;
        allreg{2} = clu_m1TJ;
    else
        clu_m1TJ  = 1:size(params.cluid{1,1}, 1);
        clu_ALM   = size(params.cluid{1,1}, 1)+1:Ncells;
        allreg{1} = clu_m1TJ;
        allreg{2} = clu_ALM;
    end

    numClu = allreg{kk};

    % ---- condtrix_kk with trial limits ----
    condtrix_kk = params(sessix).trialid{conds2use(1)};
    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
        condtrix_kk(condtrix_kk > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix_kk(condtrix_kk > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix_kk(condtrix_kk > 298) = [];
    elseif strcmp(obj.pth.dt,'2025-06-05') && strcmp(obj.pth.anm,'TDl4')
        condtrix_kk(condtrix_kk > 123) = [];
    end

% ---- Neural data ----
    condpsth = obj(sessix).trialdat(:, numClu, condtrix_kk);
    tnsorp1  = squeeze(mean(condpsth, 2));

    tnsorp_r1 = tnsorp1(:, P8);
    tnsorp_r4 = tnsorp1(:, P9);

    % ---- Sort R1 by last kinematic event ----
    tol = 1e-9;
    K1  = kin.dat(:, condtrix_kk, kinix);
    K1  = K1(:, P8);
    K1  = K1 - mode(K1(:));
    last_nz_r1 = max((1:size(K1,1)).' .* (abs(K1) > tol), [], 1);
    [~, idx1]  = sort(last_nz_r1);
    last_nz_r1 = last_nz_r1(idx1);
    tnsorp_r1  = tnsorp_r1(:, idx1);
    sorted_P8  = P8(idx1);

    lastlickr1 = zeros(1, size(tnsorp_r1, 2));
    for i = 1:size(tnsorp_r1, 2)
        if last_nz_r1(i) ~= 0
            lastlickr1(i) = obj.time(last_nz_r1(i));
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

    % ---- Sort R4 by last kinematic event ----
    K4  = kin.dat(:, condtrix_kk, kinix);
    K4  = K4(:, P9);
    K4  = K4 - mode(K4(:));
    last_nz_r4 = max((1:size(K4,1)).' .* (abs(K4) > tol), [], 1);
    [~, idx1]  = sort(last_nz_r4);
    last_nz_r4 = last_nz_r4(idx1);
    tnsorp_r4  = tnsorp_r4(:, idx1);
    sorted_P9  = P9(idx1);

    lastlickr4 = zeros(1, size(tnsorp_r4, 2));
    rewardr4   = zeros(1, size(tnsorp_r4, 2));
    for i = 1:size(tnsorp_r4, 2)
        temp = sorted_P9(i);
        if last_nz_r4(i) ~= 0
            lastlickr4(i) = obj.time(last_nz_r4(i));
            rewardr4(i)   = obj.bp.ev.reward(temp) - obj.bp.ev.goCue(temp);
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

    tnsorp_r1 = abs(tnsorp_r1);
    tnsorp_r4 = abs(tnsorp_r4);

    % ---- Colormap ----
    numColors   = 100;
    blue        = [0/255, 0/255, 153/255];
    white       = [1, 1, 1];
    t           = linspace(0, 1, numColors).^1.2;
    blueToWhite = [blue(1)+(white(1)-blue(1))*t', ...
                   blue(2)+(white(2)-blue(2))*t', ...
                   blue(3)+(white(3)-blue(3))*t'];

    % ---- Combine & per-trial baseline subtract [-1, 0) s ----
    combined_data = [tnsorp_r1'; tnsorp_r4'];
    nR1 = size(tnsorp_r1, 2);
    nR4 = size(tnsorp_r4, 2);

    baselineIdx = obj.time >= -1 & obj.time < 0;
    baselineMu  = mean(combined_data(:, baselineIdx), 2);
    dataZero    = combined_data - baselineMu;
    dataZero(dataZero < 0) = 0;

    posMax = prctile(dataZero(:), 99);
    if posMax == 0; posMax = 1; end

    % ============================================================
    %  PLOT
    % ============================================================
    imagesc(ax, obj.time, 1:(nR1+nR4), dataZero);
    colormap(ax, blueToWhite);
    caxis(ax, [0 posMax]);
    hold(ax, 'on');

    % White divider between R1 and R4
    yline(ax, nR1 + 0.5, 'w-', 'LineWidth', 1.5);

    % Red go-cue line
    xline(ax, 0, 'Color', [1 0 0 0.5], 'LineWidth', 1);

    % Orange last-lick lines R1
    orange = [1 0.5 0.1];
    for i = 1:nR1
        line(ax, [lastlickr1(i) lastlickr1(i)], [i-0.5 i+0.5], ...
            'Color', orange, 'LineWidth', 1.5);
    end

    % Orange last-lick lines R4
    for i = 1:nR4
        y_pos = nR1 + i;
        line(ax, [lastlickr4(i) lastlickr4(i)], [y_pos-0.5 y_pos+0.5], ...
            'Color', orange, 'LineWidth', 1.5);
    end

    % ---- Axis formatting ----
    axis(ax, 'tight');
    set(ax, 'Visible', 'on');
    xlim(ax, [-0.5 3]);
    box(ax, 'off');
    grid(ax, 'off');
    set(ax, 'FontSize', 9, 'TickDir', 'out');
    set(ax, 'YTick', []);

    % ---- x-label on bottom row only ----
    if sessWithinAnimal == nSessPerAnimal
        xlabel(ax, 'Time from Go Cue (s)', 'FontSize', 12);
    end

    % ---- y-label: trial count on EVERY plot ----
    ylabel(ax, sprintf('n = %d trials', nR1 + nR4), 'FontSize', 12);

    % ---- Per-plot colorbar spanning exact full axes height ----
    drawnow;
    axPos = ax.Position;
    cbW   = 0.008;
    cbGap = 0.003;
    cbPos = [axPos(1) + axPos(3) + cbGap, ...
             axPos(2), ...
             cbW, ...
             axPos(4)];
    cb = colorbar(ax, 'Position', cbPos);
    cb.Ticks      = [0, posMax];
    cb.TickLabels = {sprintf('%.1f Hz', 0), sprintf('%.1f Hz', posMax)};
    cb.FontSize   = 9;

    drawnow;
    
end  % sessnum loop
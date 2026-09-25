%% F3D_rasterPSTH.m
%  Rasters and trial-averaged spike rates for example single units, Double Reward Task.
%  One session, both probes. Trials are grouped by reward condition and the
%  raster is drawn from finely binned single-trial spike rates.
%  Produces Fig. 3D.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behaviour
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
clear,clc

sz = 14;

v2Root = fileparts(fileparts(mfilename('fullpath')));
if isempty(v2Root) || ~exist(fullfile(v2Root, 'shared', 'slimToLegacy.m'), 'file')
    v2Root = 'C:\Users\LabTech\Documents\Cortical Disengagement Figures\MATLAB Codes _ v2';
end
addpath(fullfile(v2Root, 'shared'));
addpath(fullfile(v2Root, 'shared', 'pipelineCopies'));

%% FIGURE OPTIONS
% Which unit to draw in the per-neuron raster / spike-rate figure.
% Leave empty to skip that figure. The overview grid drawn first labels every
% unit with its index, so run once, pick an index off the grid, and set it
% here. A vector draws one figure per entry.
plotNeuron = [];

%% PARAMETERS
params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping only operates on neural data for now.
params.behav_only = 0;
params.timeWarp            = 0;   % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20;   % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.01;   % remove clusters with firing rates across all trials less than this val

params.condition(1) = {'hit==1 | hit==0' };   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 6'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 6'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 6'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 1'};   % left to right         % right hits, no stim, aw off
params.condition(end+1) = {'hit==1 & rewardedLick == 6'};   % left to right         % right hits, no stim, aw off
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

params.tmin = -2.5;
params.tmax = 4;
params.dt = 1/300;   % finest rate in the exported files

% smooth with causal gaussian kernel
params.smooth = 1;

% cluster qualities to use
params.quality = {'good'};   % accepts any cell array of strings - special character 'all' returns clusters of any quality

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

meta = [];

animal = 'YH1';
date = '2023-05-10';
meta = slimMeta(animal, date);

params.probe = {meta.probe};

%% LOAD DATA

% obj, params and kin come from the exported session files, rebuilt for this
% script's own params by slimToLegacy. Replaces loadSessionData +
% loadMotionEnergy + getKinematics.
[obj, params, kin] = slimToLegacy(meta, params);

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

%% TONGUE
conds2use = [1];   % With reference to 'params.condition'
kinfeat = 'tongue_length';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
% tongue_length
% top_tongue_xdisp_view2
% jaw_ydisp_view1
sessix = 1;

psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};   % Get the trials from this condition
    condpsth = obj(sessix).trialdat(:,:,condtrix);   % Take the single trial PSTHs for these trials
end

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));

Ncells = size(obj.psth, 2);
clu_m1TJ = 1:size(params.cluid{1, 1},1);
clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;

% first lick mode
reg = clu_ALM;
reg1 = clu_m1TJ;

% if obj.pth.anm == 'TD25d'
% condtrix = condtrix(1:199);

Kinematics1 = kin.dat(:,condtrix,kinix);

conds2use = [1];   % With reference to 'params.condition'
kinfeat = 'tongue_length';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
sessix = 1;

psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};   % Get the trials from this condition
    condpsth = obj(sessix).trialdat(:,:,condtrix);   % Take the single trial PSTHs for these trials
end
% if obj.pth.anm == 'TD25d'
% condtrix = condtrix(1:199);

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
Length = kin.dat(:,condtrix,kinix);

kinfeat = 'tongue_angle';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1

psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};   % Get the trials from this condition
    condpsth = obj(sessix).trialdat(:,:,condtrix);   % Take the single trial PSTHs for these trials
end
% if obj.pth.anm == 'TD25d'
% condtrix = condtrix(1:199);

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
angle = kin.dat(:,condtrix,kinix);

kinfeat = 'top_tongue_xvel_view2';   % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1

psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};   % Get the trials from this condition
    condpsth = obj(sessix).trialdat(:,:,condtrix);   % Take the single trial PSTHs for these trials
end

% if obj.pth.anm == 'TD25d'
% condtrix = condtrix(1:199);

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
velocity = kin.dat(:,condtrix,kinix);


% C = Length';
% t = obj.time;            % 1×T
% nTrials = size(C,1);     % number of rows in C


numColors = 100;
blue = [0/255, 0/255, 153/255]; white = [255/255, 255/255, 255/255]; red = [255/255, 51/255, 0/255];
blueWhiteRed = [linspace(blue(1), white(1), numColors)', linspace(blue(2), white(2), numColors)', linspace(blue(3), white(3), numColors)';
                linspace(white(1), red(1), numColors)', linspace(white(2), red(2), numColors)', linspace(white(3), red(3), numColors)'];

numColors = 500; blue = [0/255, 0/255, 153/255]; white = [255/255, 255/255, 255/255];

blueToWhite = [linspace(blue(1), white(1), numColors)', linspace(blue(2), white(2), numColors)', linspace(blue(3), white(3), numColors)'];

% colormap(blueToWhite);colorbar;
% figure;imagesc(linspace(0, 1, numColors*2), [0 1], reshape(blueToWhite, [1, numColors, 3]));axis off;

numColors = 100; black = [0/255, 0/255, 0/255]; red = [255/255, 0/255, 0/255]; yellow = [255/255, 255/255, 0/255];

blackToRedToYellow = [linspace(black(1), red(1), numColors)', linspace(black(2), red(2), numColors)', linspace(black(3), red(3), numColors)';
                      linspace(red(1), yellow(1), numColors)', linspace(red(2), yellow(2), numColors)', linspace(red(3), yellow(3), numColors)'];

% colormap(blackToRedToYellow);colorbar;
% figure;imagesc(linspace(0, 1, numColors*2), [0 1], reshape(blackToRedToYellow, [1, numColors*2, 3]));axis off;

% ---- condition colours ------------------------------------------------
% Double Reward. The CI band is a pale version of the same hue so the
% mean line stays readable on top of it.
pale   = @(c) 1 - 0.55*(1 - c);
col_c1 = [0.90 0.45 0.70];   % R1 pink
col_c2 = [0.45 0.15 0.65];   % R16 purple

%% BIG PLOT %%

nNeurons = size(obj.trialdat,2);
if ~isempty(plotNeuron) && any(plotNeuron < 1 | plotNeuron > nNeurons | plotNeuron ~= round(plotNeuron))
    error('plotNeuron must be whole numbers between 1 and %d; this session has %d units.', nNeurons, nNeurons);
end

nCols = 15;
nRows = ceil(nNeurons / nCols);

% sample trials once
r1 = params.trialid{8};
r4 = params.trialid{9};

% r1(r1>278) = [];
% r4(r4>278) = [];

idx1 = randperm(numel(r1), 40);
idx4 = randperm(numel(r4), 40);
trials8 = r1(idx1);
trials9 = r4(idx4);

windowSize = 40;
time       = obj.time(:);

% make the figure and layout
figure('Units','normalized','Position',[0 0 1 1]);
t = tiledlayout(nRows, nCols, ...
    'TileSpacing','compact','Padding','compact');

for nn = 1:nNeurons
  ax = nexttile;   % get the next tile

  spkM8 = squeeze(obj.trialdat(:, nn, trials8));
  spkM9 = squeeze(obj.trialdat(:, nn, trials9));
  frM8 = movmean(spkM8, windowSize, 1);
  frM9 = movmean(spkM9, windowSize, 1);

% compute mean±95% CI
  mean8 = mean(frM8,2);
  err8  = 1.96 * std(frM8,0,2) / sqrt(size(frM8,2));
  mean9 = mean(frM9,2);
  err9  = 1.96 * std(frM9,0,2) / sqrt(size(frM9,2));

% plot
  hold(ax,'on');
  fill(ax, [time; flipud(time)], ...
           [mean8-err8; flipud(mean8+err8)], ...
           pale(col_c1), 'EdgeColor','none','FaceAlpha',0.9);
  plot(ax, time, mean8, '-', 'Color', col_c1, 'LineWidth',1);
  fill(ax, [time; flipud(time)], ...
           [mean9-err9; flipud(mean9+err9)], ...
           pale(col_c2), 'EdgeColor','none','FaceAlpha',0.9);
  plot(ax, time, mean9, '-', 'Color', col_c2, 'LineWidth',1);
  xline(ax, 0,'k--','LineWidth',0.5);
  xlim(ax,[-0.75 1.5]);

% tidy each small axes
  ax.XTick = [];
  ax.YTick = [];
  ax.Box   = 'off';
  title(ax, num2str(nn), 'FontSize',6);
end

% delete any completely unused tiles
totalTiles = nRows * nCols;
if totalTiles > nNeurons
% t.Children is in reverse plotting order
  for k = 1:(totalTiles - nNeurons)
    delete( t.Children(1) )
  end
end

% now set one global X‐ and Y‐label on the tiledlayout
t.XLabel.String = ['time from ' num2str(params.alignEvent)];
t.XLabel.FontSize = 12;
t.YLabel.String = 'firing rate (spikes/s)';
t.YLabel.FontSize = 12;

% 1) Select trials — SAME trials, SAME order, used for lick raster,
%    spike raster, AND PSTH for both R1 (blue) and R4 (red)
r1List = params.trialid{8};   % R1 trial IDs
r4List = params.trialid{9};   % R4 trial IDs

nDraw = 50;

idx1 = randperm(numel(r1List), nDraw);
idx4 = randperm(numel(r4List), nDraw);
sel1 = r1List(idx1);
sel4 = r4List(idx4);

nTrials = nDraw + nDraw;   % total rows in raster (R1 stacked below, R4 above)

% 2) Per-neuron, 4-row figure (3 data rows + 1 thin label row)
if isempty(plotNeuron)
    fprintf(['plotNeuron is empty, so the per-neuron figure was skipped. ' ...
             'Pick an index from the overview grid and set plotNeuron at the top.\n']);
end

for nn = plotNeuron(:)'

px     = 75;
py     = 75;
width  = 380;
height = 500;
ms     = 1.1;

neuron = nn;
time   = obj.time(:);

xLimVals    = [-0.2, 1.95];
sbX         = 1.78;
rasterScale = 50;   % SAME scale-bar length for both rasters
rasterYLim  = [0, nTrials+1];

rasterW = 0.60;
rasterH = 0.24;
rasterL = 0.18;

figure('Position',[px, py, width, height], 'Color','w');

% --- Row 1: Lick raster (R1 = blue, R4 = red, stacked) --------------------
ax1 = subplot('Position', [rasterL 0.70 rasterW rasterH]);
hold on
for i = 1:nDraw
    trial = sel1(i);
    lickContacts = obj.bp.ev.lickL{trial,1} - obj.bp.ev.goCue(trial);
    plot(lickContacts, i*ones(size(lickContacts)), '.', 'MarkerSize', ms, 'Color','b');
end
for i = 1:nDraw
    trial = sel4(i);
    lickContacts = obj.bp.ev.lickL{trial,1} - obj.bp.ev.goCue(trial);
    plot(lickContacts, (nDraw+i)*ones(size(lickContacts)), '.', 'MarkerSize', ms, 'Color','r');
end
xline(0, '--k', 'LineWidth', 1);
xlim(xLimVals);
ylim(rasterYLim);

text(-0.32, nTrials/2, 'Licks', 'Color','k', 'FontWeight','bold', ...
    'FontSize', 13, 'Rotation', 90, 'HorizontalAlignment','center')

% scale bar: rasterScale trials
plot([sbX sbX], [0, rasterScale], 'k-', 'LineWidth', 2.5)
text(sbX+0.05, rasterScale/2, [num2str(rasterScale) ' trials'], 'Color','k', ...
    'FontSize', 10, 'HorizontalAlignment','left', 'VerticalAlignment','middle')

axis off

% --- Row 2: PSTH (R1 = blue, R4 = red) — DRAGGABLE -------------------------
ax2 = subplot('Position', [rasterL 0.40 rasterW rasterH]);
hold on

windowSize = 30;

% USE THE SAME sel1/sel4 trials as the rasters (instead of all r1/r4 trials)
trials8 = sel1;
trials9 = sel4;

spkM8 = squeeze(obj.trialdat(:, neuron, trials8));
spkM9 = squeeze(obj.trialdat(:, neuron, trials9));

frM8 = movmean(spkM8, windowSize, 1);
frM9 = movmean(spkM9, windowSize, 1);

mean8 = mean(frM8, 2);
sem8  = std(frM8, 0, 2) ./ sqrt(size(frM8,2));
mean9 = mean(frM9, 2);
sem9  = std(frM9, 0, 2) ./ sqrt(size(frM9,2));

err8 = 1.96 * sem8;
err9 = 1.96 * sem9;

yMax = max([mean8+err8; mean9+err9]) * 1.15;

hFill8 = fill([time; flipud(time)], [mean8-err8; flipud(mean8+err8)], ...
    pale(col_c1), 'EdgeColor','none', 'FaceAlpha', 0.9);
hRate8 = plot(time, mean8, '-', 'Color', col_c1, 'LineWidth', 2);

hFill9 = fill([time; flipud(time)], [mean9-err9; flipud(mean9+err9)], ...
    pale(col_c2), 'EdgeColor','none', 'FaceAlpha', 0.9);
hRate9 = plot(time, mean9, '-', 'Color', col_c2, 'LineWidth', 2);

xline(0, '--k', 'LineWidth', 1);
xlim(xLimVals);
ylim([0, yMax])

text(-0.32, yMax*0.5, 'Rate', 'Color','k', 'FontWeight','bold', ...
    'FontSize', 13, 'Rotation', 90, 'HorizontalAlignment','center')

% scale bar: Hz
yScale = 10;
plot([sbX sbX], [0, yScale], 'k-', 'LineWidth', 2.5)
text(sbX+0.05, yScale/2, [num2str(yScale) ' Hz'], 'Color','k', 'FontSize', 10, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle')

axis off

% --- Make R1 and R4 traces independently draggable (with their fills) -----
makeLineDraggable(hRate8, hFill8, time);
makeLineDraggable(hRate9, hFill9, time);

% --- Row 3: Spike raster (R1 = blue, R4 = red) — IDENTICAL to Row 1 -------
ax3 = subplot('Position', [rasterL 0.10 rasterW rasterH]);
hold on
for i = 1:nDraw
    trial = sel1(i);
    spkTimes = time(obj.trialdat(:, neuron, trial)~=0);
    plot(spkTimes, i*ones(size(spkTimes)), '.', 'MarkerSize', ms, 'Color','b');
end
for i = 1:nDraw
    trial = sel4(i);
    spkTimes = time(obj.trialdat(:, neuron, trial)~=0);
    plot(spkTimes, (nDraw+i)*ones(size(spkTimes)), '.', 'MarkerSize', ms, 'Color','r');
end
xline(0, '--k', 'LineWidth', 1);
xlim(xLimVals);
ylim(rasterYLim);   % EXACTLY same as ax1

text(-0.32, nTrials/2, 'Spikes', 'Color','k', 'FontWeight','bold', ...
    'FontSize', 13, 'Rotation', 90, 'HorizontalAlignment','center')

plot([sbX sbX], [0, rasterScale], 'k-', 'LineWidth', 2.5)
text(sbX+0.05, rasterScale/2, [num2str(rasterScale) ' trials'], 'Color','k', ...
    'FontSize', 10, 'HorizontalAlignment','left', 'VerticalAlignment','middle')

axis off

% --- Row 4: thin label-only axes for "Go Cue" / "0.25 s" -------------------
ax4 = subplot('Position', [rasterL 0.01 rasterW 0.07]);
hold on
xlim(xLimVals);
ylim([0 1]);

plot([0.6 0.85], [0.8 0.8], 'k-', 'LineWidth', 2.5)
text(0.72, 0.5, '0.25 s', 'Color','k', 'FontSize', 11, ...
    'HorizontalAlignment','center', 'VerticalAlignment','top')

text(0, 0.5, 'Go Cue', 'Color','k', 'FontSize', 12, ...
    'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','top')

axis off

sgtitle(['Neuron ', num2str(nn), '  (R1 = blue, R4 = red)'], 'FontSize', 12)

end

% Helper: make a PSTH line + its CI fill draggable together (vertical offset)
function makeLineDraggable(hLine, hFill, time)
    origY      = hLine.YData;
    origFillY  = hFill.YData;
    hLine.ButtonDownFcn = @(src,evt) startDrag(src, hFill, evt, origY, origFillY);
end

function startDrag(hLine, hFill, ~, origY, origFillY)
    fig = ancestor(hLine, 'figure');
    startPoint = get(hLine.Parent, 'CurrentPoint');
    startY = startPoint(1,2);

    fig.WindowButtonMotionFcn = @(src,evt) dragLine(hLine, hFill, origY, origFillY, startY, hLine.Parent);
    fig.WindowButtonUpFcn = @(src,evt) stopDrag(fig);
end

function dragLine(hLine, hFill, origY, origFillY, startY, ax)
    currentPoint = get(ax, 'CurrentPoint');
    currentY = currentPoint(1,2);
    dy = currentY - startY;
    hLine.YData = origY + dy;
    hFill.YData = origFillY + dy;
end

function stopDrag(fig)
    fig.WindowButtonMotionFcn = '';
    fig.WindowButtonUpFcn = '';
end


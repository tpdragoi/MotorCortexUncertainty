% Finding "Kinematic Modes"
clear,clc

sz = 14;

%
% d = 'C:\Users\vdragoi\Desktop\uninstructedMovements_v2-main';
% addpath(genpath(fullfile(d,'utils')))
% addpath(genpath(fullfile(d,'DataLoadingScripts')))
% addpath(genpath(fullfile(d,'funcs')))
% rmpath(genpath(fullfile(d,'fig1')));
% addpath 'C:\Users\vdragoi\Desktop\uninstructedMovements_v2-main\ObjVis\warp'

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
params.tmax = 5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 35;

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




% % TD1d R1
% date = '2024-07-09';
% meta1 = loadTD10s_neur(meta1,datapth,date);
% date = '2024-07-10';
% meta2 = loadTD10s_neur(meta2,datapth,date);
date = '2024-07-11';
meta3 = loadTD10s_neur(meta3,datapth,date);
% % date = '2024-07-13';
% % meta4 = loadTD10s_neur713(meta4,datapth,date);
% % date = '2024-07-14';
% % meta5 = loadTD10s_neur(meta5,datapth,date);
% 
% % TD4d R1
% % date = '2024-07-05';
% % meta6 = loadTD9s_neur(meta6,datapth,date);
% % date = '2024-07-06';
% % meta7 = loadTD9s_neur(meta7,datapth,date);
% date = '2024-07-07';
% meta8 = loadTD9s_neur(meta8,datapth,date);
% % date = '2024-07-08';
% % meta9 = loadTD9s_neur(meta9,datapth,date);
date = '2024-07-09';
meta11 = loadTD9s_neur709(meta11,datapth,date);
% % date = '2024-07-10';
% % meta12 = loadTD9s_neur(meta12,datapth,date);% TD13d R1
% 
% 
% % date = '2025-02-18';
% % meta13 = loadTD3l_neur(meta13,datapth,date);
% % date = '2025-02-19';
% % meta14 = loadTD3l_neur(meta14,datapth,date);
% % date = '2025-02-20';
% % meta15 = loadTD3l_neur(meta15,datapth,date);
% % date = '2025-02-22';
% % meta16 = loadTD3l_neur222(meta16,datapth,date);
% % 
% % TD15d R1
% % date = '2025-02-18';
% % meta17 = loadTD2l_neur(meta17,datapth,date);
% % date = '2025-02-19';
% % meta18 = loadTD2l_neur(meta18,datapth,date);
% % date = '2025-02-20';
% % meta19 = loadTD2l_neur(meta19,datapth,date);
% % date = '2025-02-21';
% % meta20 = loadTD2l_neur(meta20,datapth,date);
% % date = '2025-02-24';
% % meta21 = loadTD2l_neur(meta21,datapth,date);
% 
% 
% 
% 
% % date = '2025-07-25';
% % meta22 = loadTD27_neur(meta22,datapth,date);
% % date = '2025-07-26';
% % meta23 = loadTD27_neur(meta23,datapth,date);
% % date = '2025-07-27';
% % meta24 = loadTD27_neur(meta24,datapth,date);
% date = '2025-07-28';
% meta25 = loadTD27_neur(meta25,datapth,date);
% date = '2025-07-29';
% meta26 = loadTD27_neur(meta26,datapth,date);
% % date = '2025-07-30';
% % meta27 = loadTD27_neur(meta27,datapth,date);
% % date = '2025-07-31';
% % meta28 = loadTD27_neur(meta28,datapth,date);


% date = '2025-07-30';
% meta29 = loadTD26_neur(meta29,datapth,date);
% date = '2025-07-31';
% meta30 = loadTD26_neur(meta30,datapth,date);
% date = '2025-08-01';
% meta31 = loadTD26_neur(meta31,datapth,date);
% date = '2025-08-02';
% meta32 = loadTD26_neur222(meta32,datapth,date);
% date = '2025-08-03';
% meta33 = loadTD26_neur222(meta33,datapth,date);
% date = '2025-08-05';
% meta34 = loadTD26_neur(meta34,datapth,date);
% date = '2025-08-06';
% meta35 = loadTD26_neur(meta35,datapth,date);
% date = '2025-08-07';
% meta36 = loadTD26_neur(meta36,datapth,date);

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



for sessnum = 1:length(all_meta)

clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin 

meta = all_meta(sessnum,1); 
params.probe = {meta.probe};

% LOAD DATA
clear obj 
params.cluid = {};
[obj,params] = loadSessionData(meta,params,params.behav_only);
trialSet = [1:obj.bp.Ntrials]';

for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end

% Get kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix),{' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
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
    Kin_r1  = Kin_all(:, P8_local);                  % [time x nR1]
    Kin_r4  = Kin_all(:, P9_local);                  % [time x nR4]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort R1 trials by last lick time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort R4 trials by last lick time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Colormaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    p_all   = prctile(allData, [1 99]);
    rang    = [p_all(1) p_all(2)];
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



end % kk loop

end % sessnum loop

a = 3;



% r1Trials = params.trialid{4};
% r4Trials = params.trialid{5};
% 
% condpsth_x = condpsth(:,numClu,r1Trials);
% mean_values = mean(condpsth_x, [1, 2]);
% std_values = std(condpsth_x, 0, [1, 2]);  % 0 indicates normalization by N-1 (sample standard deviation)
% tnsorp1 = squeeze(mean(condpsth_x, 2));
% 
% figure; plot(obj.time, movmean(mean(tnsorp1,2),20))
% 
% r1Trials = params.trialid{4};
% r4Trials = params.trialid{5};
% 
% condpsth_x = condpsth(:,numClu,r4Trials);
% mean_values = mean(condpsth_x, [1, 2]);
% std_values = std(condpsth_x, 0, [1, 2]);  % 0 indicates normalization by N-1 (sample standard deviation)
% tnsorp1 = squeeze(mean(condpsth_x, 2));
% 
% hold on; plot(obj.time, movmean(mean(tnsorp1,2),20))



%%

% r1Trials = params.trialid{4};
% r4Trials = params.trialid{5};
% 
% % condpsth: [time x trials x neurons] => 1300 x 200 x 86
% 
% % Extract data for trial 4 across all neurons
% trial4 = condpsth(:, 8, :);  % size = 1300 x 1 x 86
% 
% % Squeeze to 1300 x 86
% trial4 = squeeze(trial4);  % size = 1300 x 86
% 
% % Compute mean across neurons (dimension 2)
% mean_firing = mean(trial4, 2);  % size = 1300 x 1
% 
% % Optional: create time vector (e.g., assuming 5 ms per sample)
% % t = (0:1299) * 0.005 - 6.5;  % adjust if you know the actual sampling
% 
% % Plot
% figure;
% plot(t, mean_firing, 'k', 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Mean Firing Rate (a.u.)');
% title('Mean Firing Rate Across Neurons - Trial 4');
% grid on;

%%


%%
% 
% figure;
% subplot(2,1,1)
% 
% periodN = 'firstLick';
% imagesc(obj.time, [1:size(tnsorp_r1',1)], tnsorp_r1');  % use masked data
% colormap(cmp)
% caxis(rang);  % includes -1 for black
% colorbar; grid off; axis tight; 
% box off;
% xlabel(['time from ', num2str(params.alignEvent)]);
% ylabel('Trial #'); 
% xlim([-0.3 3])
% title(['MSR R1 Date : ', num2str(obj.pth.dt), ' Anm ', num2str(obj.pth.anm)]);
% xline(0, 'Color', [1 0 0 0.5], 'LineWidth', 2)
% set(gca, 'FontSize', sz)
% 
% orange = [1 0.5 0.1];
% 
% 
% for i = 1:size(tnsorp_r1',1)
%     yval(i) = i;
%     xval(i) = lastlickr1(i);
% hold on
%     line([xval(i), xval(i)], [yval(i)-0.5, yval(i)+0.5], ...
%         'LineWidth', 3, 'Color', orange);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % rang = [-7 7];
% sz = 20;
% 
% %%%%%%%%% Masking data after last lick time by setting to -1 %%%%%%%%%%%
% masked_r4 = tnsorp_r4';  % Transpose so trials x time
% 
% for i = 1:size(masked_r4,1)
%     if ~isnan(lastlickr4(i))  % Make sure the time is valid
%         mask_idx = obj.time > lastlickr4(i);  % Find timepoints after last lick
%         masked_r4(i, mask_idx) = -1;  % Set to -1 for black
%     end
% end
% 
% %%%%%%%%% Plotting %%%%%%%%%%%
% 
% % Extend colormap by adding black at the beginning
% 
% subplot(2,1,2)
% 
% periodN = 'firstLick';
% imagesc(obj.time, [1:size(tnsorp_r4',1)], tnsorp_r4');  % use masked data
% colormap(cmp)
% caxis(rang);  % includes -1 for black
% colorbar; grid off; axis tight; 
% box off;
% xlabel(['time from ', num2str(params.alignEvent)]);
% ylabel('Trial #'); 
% xlim([-0.3 3])
% title(['MSR R1 Date : ', num2str(obj.pth.dt), ' Anm ', num2str(obj.pth.anm)]);
% xline(0, 'Color', [1 0 0 0.5], 'LineWidth', 2)
% set(gca, 'FontSize', sz)
% 
% orange = [1 0.5 0.1];
% 
% for i = 1:size(masked_r4,1)
%     yval(i) = i;
%     xval(i) = lastlickr4(i);
%     line([xval(i), xval(i)], [yval(i)-0.5, yval(i)+0.5], ...
%         'LineWidth', 3, 'Color', orange);
% end



a = 2;

%%

% orange = [1 0.5 0.1 ];
% blue = [0 0.31 1 ];
% 
% 
% % me1 = MotionEnergy(:,params.trialid{1, 4});
% % me1 = normalize(me1,'range',[0,14.5]);
% % me2 = MotionEnergy(:,params.trialid{1, 5});
% % me2 = normalize(me2,'range',[0,14.5]);
% 
% y = tnsorp_r1;
% yhat = tnsorp_r4;
% 
% % y = normalize(y,'range',[0,7]);
% % yhat = normalize(yhat,'range',[0,7]);
% 
% figure;
% ax = nexttile;
% % ax = subplot(1,1,1);
% 
% 
% 
% alph = 0.2;
% 
% windowSize = 1;
% 
% % mu.me1 = nanmedian(me1,2);
% % sd.me1 = nanstd(me1,[],2) ./ sqrt(size(me1,2));
% % mu.me2 = nanmedian(me2,2);
% % sd.me2 = nanstd(me2,[],2) ./ sqrt(size(me2,2));
% 
% % mu.me1 = mu.me1 - 4;
% % mu.me2 = mu.me2 - 4;
% 
% % hold on
% % shadedErrorBar(obj.time,movmean(mu.me1,windowSize),sd.me1,{'--','Color',[blue 0.2],'LineWidth',2},alph,ax)
% % hold on
% % shadedErrorBar(obj.time,movmean(mu.me2,windowSize),sd.me2,{'--','Color',[orange 0.2],'LineWidth',2},alph,ax)
% 
% 
% 
% alph = 0.7;
% 
% hold on;
% 
% 
% windowSize = 1;
% mu.y = nanmedian(y,2);
% mu.yhat = nanmedian(yhat,2);
% sd.y = nanstd(y,[],2) ./ sqrt(size(y,2));
% sd.yhat = nanstd(yhat,[],2) ./ sqrt(size(yhat,2));
% hold on
% shadedErrorBar(obj.time,movmean(mu.y,windowSize),sd.y,{'Color',[blue 1],'LineWidth',2},alph,ax)
% shadedErrorBar(obj.time,movmean(mu.yhat,windowSize),sd.yhat,{'Color',[orange 1],'LineWidth',2},alph,ax)
% 
% 
% % ylabel(par.feats,'Interpreter','none')
% ylabel('Projection (a.u.)')
% xlabel(['time from ', num2str(params.alignEvent)]);
% set(gca,'FontSize',sz)
% axis tight;box off
% % ylim([0.3 7])
% xlim([-1. 3])
% 
% 

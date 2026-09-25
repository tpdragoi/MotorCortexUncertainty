% Finding "Kinematic Modes"
clear,clc

sz = 8;

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
meta37 = [];



% meta1 = loadTD_inactivation1(meta1,datapth);
% meta2 = loadTD_inactivation2(meta2,datapth);
% meta3 = loadTD_inactivation3(meta3,datapth);
% meta4 = loadTD_inactivation4(meta4,datapth);
% meta5 = loadTD_inactivation5(meta5,datapth);




% % TD1d R1
% date = '2023-02-21';
% meta1 = loadTD1_neural(meta1,datapth,date);
% date = '2023-02-22';
% meta2 = loadTD1_neural(meta2,datapth,date);
% date = '2023-02-23';
% meta3 = loadTD1_neural(meta3,datapth,date);
% date = '2023-02-24';
% meta4 = loadTD1_neural(meta4,datapth,date);
% % 
% % % TD4d R1
% date = '2023-02-21';
% meta5 = loadTD4_neural(meta5,datapth,date);
% date = '2023-02-24';
% meta6 = loadTD4_neural(meta6,datapth,date);
% date = '2023-02-25';
% meta7 = loadTD4_neural(meta7,datapth,date);
% date = '2023-03-19';
% meta8 = loadTD4_neural(meta8,datapth,date);
% 
% % TD13d R1
% date = '2024-11-11';
% meta9 = loadTD13_neural(meta9,datapth,date);
% date = '2024-11-12';
% meta11 = loadTD13_neural(meta11,datapth,date);
% date = '2024-11-13';
% meta12 = loadTD13_neural(meta12,datapth,date);
% date = '2024-11-21';
% meta13 = loadTD13_neural(meta13,datapth,date);
% date = '2024-11-22';
% meta14 = loadTD13_neural(meta14,datapth,date);
% date = '2024-11-24';
% meta15 = loadTD13_neural(meta15,datapth,date);
% date = '2024-11-25';
% meta16 = loadTD13_neural(meta16,datapth,date);
% % 
% % TD15d R1
% date = '2024-11-24';
% meta17 = loadTD15_neural(meta17,datapth,date);
% date = '2024-11-25';
% meta18 = loadTD15_neural(meta18,datapth,date);
% date = '2024-11-26';
% meta19 = loadTD15_neural(meta19,datapth,date);
% date = '2024-11-27';
% meta20 = loadTD15_neural(meta20,datapth,date);
% 
% % TD8d R1
% date = '2024-09-06';
% meta21 = loadTD8_neural(meta21,datapth,date);
% date = '2024-09-07';
% meta22 = loadTD8_neural(meta22,datapth,date);
% date = '2024-09-08';
% meta23 = loadTD8_neural(meta23,datapth,date);
% date = '2024-09-09';
% meta24 = loadTD8_neural(meta24,datapth,date);
% date = '2024-09-10';
% meta25 = loadTD8_neural(meta25,datapth,date);
% date = '2024-09-22';
% meta26 = loadTD8_neural(meta26,datapth,date);
% % 
% % 
% date = '2025-06-17';
% meta27 = loadTD22_neural(meta27,datapth,date);
% date = '2025-06-18';
% meta28 = loadTD22_neural(meta28,datapth,date);
% date = '2025-06-19';
% meta29 = loadTD22_neural(meta29,datapth,date);
% date = '2025-06-20';
% meta30 = loadTD22_neural(meta30,datapth,date);
% date = '2025-06-21';
% meta31 = loadTD22_neural(meta31,datapth,date);
% % 
% date = '2025-06-17';
% meta32 = loadTD23_neural(meta32,datapth,date);
% date = '2025-06-18';
% meta33 = loadTD23_neural(meta33,datapth,date);
% date = '2025-06-19';
% meta34 = loadTD23_neural(meta34,datapth,date);
% date = '2025-06-20';
% meta35 = loadTD23_neural(meta35,datapth,date);
% date = '2025-06-21';
% meta36 = loadTD23_neural(meta36,datapth,date);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TD1d R1
% date = '2023-02-21';
% meta1 = loadTD1_neural(meta1,datapth,date);
date = '2023-02-22';
meta2 = loadTD1_neural(meta2,datapth,date);
% date = '2023-02-23';
% meta3 = loadTD1_neural(meta3,datapth,date);
date = '2023-02-24';
meta4 = loadTD1_neural(meta4,datapth,date);
% 
% % TD4d R1
% date = '2023-02-21';
% meta5 = loadTD4_neural(meta5,datapth,date);
% date = '2023-02-24';
% meta6 = loadTD4_neural(meta6,datapth,date);
% date = '2023-02-25';
% meta7 = loadTD4_neural(meta7,datapth,date);
% date = '2023-03-19';
% meta8 = loadTD4_neural(meta8,datapth,date);

% TD13d R1
% date = '2024-11-11';
% meta9 = loadTD13_neural(meta9,datapth,date);
date = '2024-11-12';
meta11 = loadTD13_neural(meta11,datapth,date);
% date = '2024-11-13';
% meta12 = loadTD13_neural(meta12,datapth,date);
% date = '2024-11-21';
% meta13 = loadTD13_neural(meta13,datapth,date);
% date = '2024-11-22';
% meta14 = loadTD13_neural(meta14,datapth,date);
date = '2024-11-24';
meta15 = loadTD13_neural(meta15,datapth,date);
date = '2024-11-25';
meta16 = loadTD13_neural(meta16,datapth,date);

% TD15d R1
% date = '2024-11-24';
% meta17 = loadTD15_neural(meta17,datapth,date);
date = '2024-11-25';
meta18 = loadTD15_neural(meta18,datapth,date);
% date = '2024-11-26';
% meta19 = loadTD15_neural(meta19,datapth,date);
% % date = '2024-11-27';
% % meta20 = loadTD15_neural(meta20,datapth,date);

% TD8d R1
% date = '2024-09-06';
% meta21 = loadTD8_neural(meta21,datapth,date);
% date = '2024-09-07';
% meta22 = loadTD8_neural(meta22,datapth,date);
% date = '2024-09-08';
% meta23 = loadTD8_neural(meta23,datapth,date);
% date = '2024-09-09';
% meta24 = loadTD8_neural(meta24,datapth,date);
% date = '2024-09-10';
% meta25 = loadTD8_neural(meta25,datapth,date);
% date = '2024-09-22';
% meta26 = loadTD8_neural(meta26,datapth,date);
% 
% 
% date = '2025-06-17';
% meta27 = loadTD22_neural(meta27,datapth,date);
% date = '2025-06-18';
% meta28 = loadTD22_neural(meta28,datapth,date);
% date = '2025-06-19';
% meta29 = loadTD22_neural(meta29,datapth,date);
% date = '2025-06-20';
% meta30 = loadTD22_neural(meta30,datapth,date);
% date = '2025-06-21';
% meta31 = loadTD22_neural(meta31,datapth,date);
% 
% date = '2025-06-17';
% meta32 = loadTD23_neural(meta32,datapth,date);
% date = '2025-06-18';
% meta33 = loadTD23_neural(meta33,datapth,date);
% date = '2025-06-19';
% meta34 = loadTD23_neural(meta34,datapth,date);
date = '2025-06-20';
meta35 = loadTD23_neural(meta35,datapth,date);
% date = '2025-06-21';
% meta36 = loadTD23_neural(meta36,datapth,date);



y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

% all_meta = [meta1;meta2;meta3;meta4;meta5];
% all_meta = [meta2;meta1;meta3;meta1;meta1];
% all_meta = [meta1];

% all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12;meta13 ...
%     ;meta14;meta15;meta16];

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
    
% if sessnum ~= 1
%     params.trialid = [];
% end


meta = all_meta(sessnum,1); 

% meta = loadTD_inactivation(meta,datapth);
% meta = loadYH(meta,datapth);

params.probe = {meta.probe};

% LOAD DATA
clear obj 

% if sessnum > 1
%     clear params.trialid
% end
params.cluid = {};
[obj,params] = loadSessionData(meta,params,params.behav_only );
% [obj,params] = loadSessionData(meta,params);

trialSet = [1:obj.bp.Ntrials]';


for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end

% Get kinematic data
%---------------------------------------------
% kin (struct array) - one entry per session
%---------------------------------------------
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end

% clearvars -except kin meta obj params

%

conds2use = [1];                      % With reference to 'params.condition'
kinfeat = 'tongue_length';    % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
sessix = 1;

psthForProj = [];
for c = conds2use
    condtrix = trialSet;                                           % Get the trials from this condition
    %     condpsth = obj(sessix).trialdat(:,:,condtrix);                                  % Take the single trial PSTHs for these trials

    if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    condtrix(condtrix > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 298) = [];
     % elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
     %    condtrix(condtrix > 298) = [];
    end

end


kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
Length = kin.dat(:,condtrix,kinix);

% kinfeat = 'tongue_angle';    % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
% 
% psthForProj = [];
% for c = conds2use
%     condtrix = params(sessix).trialid{c};                                           % Get the trials from this condition
%     condpsth = obj(sessix).trialdat(:,:,condtrix);                                  % Take the single trial PSTHs for these trials
% end
% 
% kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
% angle = kin.dat(:,condtrix,kinix);
% 
% kinfeat = 'motion_energy';    % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
% 
% psthForProj = [];
% for c = conds2use
%     condtrix = params(sessix).trialid{c};                                           % Get the trials from this condition
%     condpsth = obj(sessix).trialdat(:,:,condtrix);                                  % Take the single trial PSTHs for these trials
% end
% 
% kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));
% % MotionEnergy = kin.dat(:,condtrix,kinix);

%%



% reg = clu_ALM;
% reg1 = clu_m1TJ;
% 
% 
% 
% sz = 10;
% NT = obj.bp.Ntrials;
% lw = 2;
% ct = 1;
% n = 20;
% 
% neurons = reg1;
% 
% col = {[1 0 0] [0 0 1] [0.5 0 0] [0 0 0.5]};
% % 2 is LR1, 4 RR1, 5 is LR16, 7 is RR16
% 
% for i = 1:n:length(neurons)
%     f = figure;
%     for j = 1:n
% 
%         try
%         ax = nexttile;
% 
% %         plot(obj.time,obj.psth(:,neurons(ct),2),'Color',col{1},'LineWidth',lw); hold on 
% %         plot(obj.time,obj.psth(:,neurons(ct),4),'Color',col{2},'LineWidth',lw); hold on 
% %         plot(obj.time,obj.psth(:,neurons(ct),5),'Color',col{3},'LineWidth',lw); hold on 
% %         plot(obj.time,obj.psth(:,neurons(ct),7),'Color',col{4},'LineWidth',lw); hold on 
% 
%         plot(obj.time,obj.psth(:,neurons(ct),4),'Color',col{1},'LineWidth',lw); hold on
% %         plot(obj.time,obj.psth(:,neurons(ct),9),'Color',col{2},'LineWidth',lw); hold on 
% 
% %         legend('LR1','LR16','RR1','RR16')
% %         plot(objs.time,objs.psth(:,neurons(ct),3),'Color','b','LineWidth',1.5); hold on 
% 
% 
% title(['Cell = ', num2str(neurons(ct))]);
% xlabel(['time from ', num2str(params.alignEvent)]);
% % xlabel('time')
% xline(-2)
% xline(0)
% xlim([-0.5 1.5]);
% set(gca,'FontSize',sz)
% % axis square
% 
% ylabel('firing rate')
% box off
%         ct = ct + 1;
%     end
% end
% end
% 
% a = 2;



%%

% figure; 
% subplot(1,2,1)
% cmp = jet;
% 
% imagesc(obj.time, [1:size(Length,2)], Length', 'Interpolation', 'bilinear');
% xline(0,'LineWidth', 1.5)
% colormap(jet)
% %     caxis([rang]);
% colorbar; grid off; axis tight; axis square;
% box off; xlabel('time(s)'); ylabel('Trial #'); title('length plot');
% set(gca,'FontSize',15)
% subplot(1,2,2)
% imagesc(obj.time, [1:size(angle,2)], angle', 'Interpolation', 'bilinear');
% xline(0,'LineWidth', 1.5)
% colormap(jet)
% %     caxis([rang]);
% colorbar; grid off; axis tight; axis square;
% box off; xlabel('time(s)'); ylabel('Trial #'); title('angle plot');
% set(gca,'FontSize',15)


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

%%

all11 = [1:obj.bp.Ntrials];
hit11 = all11((obj.bp.hit == 1));
r111 = all11((obj.bp.rewardedLick == 1));
r444 = all11((obj.bp.rewardedLick == 4));

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
end



%% TONGUE

task = 14;

conds2use = [1];                      % With reference to 'params.condition'
kinfeat = 'tongue_length';    % top_tongue_xvel_view2 | motion_energy | nose_xvel_view1 | jaw_yvel_view2 | trident_yvel_view1
% top_tongue_ydisp_view2 % topleft_tongue_ydisp_view2 % bottom_tongue_ydisp_view2 % bottomleft_tongue_ydisp_view2

sessix = 1;

Ncells = size(obj.psth, 2);
% clu_m1TJ = 1:size(params.cluid{1, 1},1);
% clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;

condpsth = obj.trialdat;

% clu_m1TJ = 1:size(params.cluid,1);
% 
% numClu = clu_m1TJ;

if strcmp(obj.pth.dt,'2024-07-13') && strcmp(obj.pth.anm,'TD10si') || strcmp(obj.pth.dt,'2024-07-09') && strcmp(obj.pth.anm,'TD9si') ...
        || strcmp(obj.pth.dt,'2025-02-22') && strcmp(obj.pth.anm,'TDl3')  || strcmp(obj.pth.dt,'2025-02-21') && strcmp(obj.pth.anm,'TDl2') ...
        || strcmp(obj.pth.dt,'2025-02-19') && strcmp(obj.pth.anm,'TDl2') || strcmp(obj.pth.dt,'2025-04-21') && strcmp(obj.pth.anm,'TD20d') ...
        || strcmp(obj.pth.dt,'2025-08-02') && strcmp(obj.pth.anm,'TD26d')  || strcmp(obj.pth.dt,'2025-08-03') && strcmp(obj.pth.anm,'TD26d')   


Ncells = size(obj.psth, 2);
clu_m1TJ = 1:size(params.cluid,1);
% clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;

reg = clu_m1TJ;
reg1 = clu_m1TJ;  

else

Ncells = size(obj.psth, 2);
clu_m1TJ = 1:size(params.cluid{1, 1},1);
clu_ALM = size(params.cluid{1, 1},1)+1:Ncells;

reg = clu_m1TJ;
reg1 = clu_ALM;  

end






allreg = {};

allreg{1} = [reg];
allreg{2} = [reg1];


%%
figure;

for kk = 1:2

% numClu = clu_ALM;
numClu = allreg{kk};

%%% Get neural data for projecting onto modes %%%
psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};            
    
        if strcmp(obj.pth.dt,'2024-11-11') && strcmp(obj.pth.anm,'TD13d')
    condtrix(condtrix > 278) = [];
    elseif strcmp(obj.pth.dt,'2024-09-07') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 313) = [];
    elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
        condtrix(condtrix > 298) = [];
     % elseif strcmp(obj.pth.dt,'2024-09-09') && strcmp(obj.pth.anm,'TD8d')
     %    condtrix(condtrix > 298) = [];
        end


    % Get the trials from this condition
    condpsth = obj(sessix).trialdat(:,:,condtrix);                                  % Take the single trial PSTHs for these trials
end

kinix =  find(strcmp(kin(sessix).featLeg,kinfeat));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Lick Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condpsth = condpsth(:,numClu,:);

tnsorp1 = squeeze(mean(condpsth, 2));

tnsorp_r1 = tnsorp1(:,P8);
tnsorp_r4 = tnsorp1(:,P9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Get Last Lick Time %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear mode

%%%%%%%% R1 %%%%%%%%%
Kinematics1 = kin.dat(:,condtrix,kinix);
Kinematics1 = Kinematics1(:, P8);

tol = 1e-9;
Kinematics1 = Kinematics1 - mode(Kinematics1);
last_nonzero_index = max((1:size(Kinematics1,1)).' .* (abs(Kinematics1) > tol), [], 1);
[~, idx1] = sort(last_nonzero_index);

last_nonzero_index = last_nonzero_index(idx1);
tnsorp_r1          = tnsorp_r1(:, idx1);
sorted_P8          = P8(idx1);

lastlickr1 = zeros(1, size(tnsorp_r1, 2));
rewardr1   = nan(1,   size(tnsorp_r1, 2));

for i = 1:size(tnsorp_r1, 2)
    temp = sorted_P8(i);
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

%%%%%%% R4 %%%%%%%%%%
Kinematics1 = kin.dat(:,condtrix,kinix);
Kinematics1 = Kinematics1(:, P9);

tol = 1e-9;
Kinematics1 = Kinematics1 - mode(Kinematics1);
last_nonzero_index = max((1:size(Kinematics1,1)).' .* (abs(Kinematics1) > tol), [], 1);
[~, idx1] = sort(last_nonzero_index);

last_nonzero_index = last_nonzero_index(idx1);
tnsorp_r4          = tnsorp_r4(:, idx1);
sorted_P9          = P9(idx1);

lastlickr4 = zeros(1, size(tnsorp_r4, 2));
rewardr4   = nan(1,   size(tnsorp_r4, 2));
rewardr41  = nan(1,   size(tnsorp_r4, 2));

for i = 1:size(tnsorp_r4, 2)
    temp = sorted_P9(i);
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
cmp    = blueToWhite;

allData = [tnsorp_r1(:); tnsorp_r4(:)];
p_all   = prctile(allData, [2 98]);
rang    = [p_all(1) p_all(2)];

% 1) R1 panel
subplot(2,2,kk)
imagesc(obj.time, 1:nR1, tnsorp_r1');
colormap(cmp); caxis(rang);
colorbar; axis tight;
xlabel('Time from Go Cue (s)')
ylabel('Trial # (R1)')
title(sprintf('Date %s  Animal %s  Probe %d (R1)', obj.pth.dt, obj.pth.anm, kk))
xlim([-1 3])
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
colorbar; axis tight;
xlabel('Time from Go Cue (s)')
ylabel('Trial # (R4)')
title(sprintf('Date %s  Animal %s  Probe %d (R4)', obj.pth.dt, obj.pth.anm, kk))
xlim([-1 3])
set(gca, 'FontSize', sz)

hold on
for i = 1:nR4
    line([lastlickr4(i) lastlickr4(i)], [i-0.5 i+0.5], 'Color', orange, 'LineWidth', 3);
    % if you want reward times uncomment:
    % r = rewardr4(i);
    % if r < 2
    %     line([r r],[i-0.5 i+0.5],'Color',pink,'LineWidth',2.5);
    % end
end
hold off

set(gcf, 'Position', [50 100 900 500]);

end

% %--- pull out your reward labels (I assume column 1 has the 1’s and 4’s) ---
% rew = obj.bp.rewardedLick(:,1);  
% trials1 = find(rew==1);   % indices of trials where label==1
% trials4 = find(rew==4);   % indices of trials where label==4
% 
% n1 = numel(trials1);
% n4 = numel(trials4);
% 
% gap = 1;  % blank row between the two chunks (optional)

% figure; 
% hold on;
% 
% %--- plot the “1” trials in blue ---
% for ii = 1:n1
%     tr = trials1(ii);
%     % align lick times to goCue of that trial
%     lt = obj.bp.ev.lickL{tr} - obj.bp.ev.goCue(tr);
%     y = ii;    % stack from 1 up to n1
%     scatter(lt, y*ones(size(lt)), 3, 'b', 'filled');
% end
% 
% %--- plot the “4” trials in red ---
% for jj = 1:n4
%     tr = trials4(jj);
%     lt = obj.bp.ev.lickL{tr} - obj.bp.ev.goCue(tr);
%     y = n1 + gap + jj;   % start just above the blue block
%     scatter(lt, y*ones(size(lt)), 3, 'r', 'filled');
% end
% 
% %--- styling ---
% % vertical line at t=0
% yl = [0, n1 + gap + n4 + 1];
% plot([0 0], yl, 'k--', 'LineWidth',1);
% 
% xlabel('Time from go cue (s)');
% ylabel('Trial (chunked by reward)');
% title('Lick Raster: blue = reward 1 trials; red = reward 4 trials');
% ylim(yl);
% xlim([-0.3 2])
% box on;
% hold off;

end

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

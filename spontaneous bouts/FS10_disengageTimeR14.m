%% FS10_disengageTimeR14.m
%  Motor cortical activity during spontaneous licking bouts.
%  Bouts that were not triggered by the go cue and were never rewarded are
%  aligned to their onset and compared with go cue-evoked bouts.
%
%  ENGAGEMENT MODE WEIGHTING. The mode is defined with d': each neuron's
%  pre-minus-post difference divided by the pooled across-trial SD of the two
%  windows, rather than the raw difference in means. See engagementModeWeights
%  at the bottom of this file.
%  READS
%    Data\<task>\<ANM>_<DATE>_obj.mat   spikes and behavior
%    Data\<task>\<ANM>_<DATE>_kin.mat   video kinematics
%    through shared\slimMeta and shared\slimToLegacy; nothing outside this folder
%  ANALYSIS SETTINGS
%    params.alignEvent  'goCue'
%    params.dt          1/300
%    params.smooth      10
%    params.quality     {'good','excellent'}
%    params.lowFR       0.01
%    params.window      -2.5 to 22 s
%  Run the whole file. Section headings below follow the order of the
%  analysis, from loading through fitting to the figures.

clear,clc

% Progress messages are silenced by default. To see them, set verbose = true
% in the logf helper at the bottom of this file.


sz = 26;

params.alignEvent          = 'goCue';   % 'fourthLick' 'goCue'  'moveOnset'  'firstLick' 'thirdLick' 'lastLick' 'reward'

% time warping only operates on neural data for now.
params.behav_only          = 0;

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
hmmRoot = fullfile(v2Root, 'Disengagement Times');   % HMM-GLM disengagement times, one folder per session and probe
assert(exist(hmmRoot, 'dir') == 7, 'No disengagement-time folder: %s', hmmRoot);

%% PARAMETERS
params.timeWarp            = 0;   % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20;   % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.01;   % minimum mean firing rate, Hz

% params.condition(1) = {'hit==1 | hit==0' };    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 1& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 2& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & trialTypes == 3& rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & rewardedLick == 1'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1 & rewardedLick == 4'};    % left to right         % right hits, no stim, aw off
% params.condition(end+1) = {'hit==1' };    % left to right         % right hits, no stim, aw off
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

params.tmin = -2.5;
params.tmax = 22;
params.dt = 1/300;

% smooth with causal gaussian kernel
params.smooth = 10;

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

%% CONSISTENCY CONFIG
% One place for the settings that were hard-coded in several spots, or that
% differ between the five scripts.

% ENGAGEMENT MODE. The Methods describe a 400 ms window either side of the
% transition; the code has always used 300 ms. Set this to 0.400 to follow the
% Methods, or change the Methods to 300 ms -- but the two have to agree.
cfg.engModeWin_s = 0.400;   % +/-400 ms around the transition

% A neuron whose rate barely varies across trials has a near-zero pooled SD and
% would otherwise dominate the weight vector. Its SD is floored here, in
% spikes/s, before the division.
cfg.sdFloor = 0.5;


% PROJECTION. Methods: p(t) = sum_i w_i r_i(t), weights normalized so
% sum|w_i| = 1. 'mean' divides that by the unit count, which is what the code
cfg.projMode = 'sum';   % 'sum' (Methods) | 'mean' (previous behavior)
cfg.plotPerSession = false;   % false = no per-session heatmap figure
% (one window per session and probe). The group
% figures at the end are made either way.

%% SPECIFY DATA TO LOAD

datapth = '';   % raw data folder not used (was: datapth = 'C:\Users\LabTech\Documents\Cortical Disengagement Code and Data\uninstructedMovements_v2-main\data';)

% one empty placeholder per session slot; the ones a task uses are
% filled in below and the rest drop out of the all_meta concatenation
[meta, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11, ...
    meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19, meta20, meta21, ...
    meta22, meta23, meta24, meta25, meta26, meta27, meta28, meta29, meta30, meta31, ...
    meta32, meta33, meta34, meta35, meta36, meta37, meta38, meta39] = deal([]);
date = '2023-02-21';
meta1 = slimMeta('loadTD1_neural', date);
date = '2023-02-22';
meta2 = slimMeta('loadTD1_neural', date);
date = '2023-02-23';
meta3 = slimMeta('loadTD1_neural', date);
date = '2023-02-24';
meta4 = slimMeta('loadTD1_neural', date);

date = '2023-02-21';
meta5 = slimMeta('loadTD4_neural', date);
% which paired each day's spikes with the NEXT day's HMM fit and dropped 03-19).
date = '2023-02-24';
meta6 = slimMeta('loadTD4_neural', date);
date = '2023-02-25';
meta7 = slimMeta('loadTD4_neural', date);
date = '2023-03-19';
meta8 = slimMeta('loadTD4_neural', date);

date = '2024-11-12';
meta9 = slimMeta('loadTD13_neural', date);
date = '2024-11-13';
meta10 = slimMeta('loadTD13_neural', date);
date = '2024-11-21';
meta11 = slimMeta('loadTD13_neural', date);

date = '2024-11-24';
meta12 = slimMeta('loadTD15_neural', date);
date = '2024-11-25';
meta13 = slimMeta('loadTD15_neural', date);
date = '2024-11-26';
meta14 = slimMeta('loadTD15_neural', date);
date = '2024-11-27';
meta15 = slimMeta('loadTD15_neural', date);

date = '2025-06-17';
meta16 = slimMeta('loadTD22_neural', date);
date = '2025-06-18';
meta17 = slimMeta('loadTD22_neural', date);
date = '2025-06-19';
meta18 = slimMeta('loadTD22_neural', date);
date = '2025-06-20';
meta19 = slimMeta('loadTD22_neural', date);
date = '2025-06-21';
meta20 = slimMeta('loadTD22_neural', date);

date = '2025-06-17';
meta21 = slimMeta('loadTD23_neural', date);
date = '2025-06-18';
meta22 = slimMeta('loadTD23_neural', date);
date = '2025-06-19';
meta23 = slimMeta('loadTD23_neural', date);
date = '2025-06-20';
meta24 = slimMeta('loadTD23_neural', date);
date = '2025-06-21';
meta25 = slimMeta('loadTD23_neural', date);

% SESSION × PROBE dataDir lookup
%  sessions  1-4  TD1d  (2023-02-21, 02-22, 02-23, 02-24)
%  sessions  5-8  TD4d  (2023-02-21, 02-24, 02-25, 03-19)
%  sessions  9-11 TD13d (2024-11-12, 11-13, 11-21)
%  sessions 12-15 TD15d (2024-11-24, 11-25, 11-26, 11-27)
%  sessions 16-20 TD22d (2025-06-17 .. 06-21)
%  sessions 21-25 TD23d (2025-06-17 .. 06-21)

% C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Disengagement_forTudor\R14_ALM

%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_21_P1',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_21_P2';  % sess 1
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_22_P1', 'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_22_P2';  % sess 2
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_23_P1',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_23_P2';  % sess 3
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_24_P1',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD1d_2023_02_24_P2';  % sess 4
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_21_P2',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_21_P2';  % sess 5
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_23_P2',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_23_P2';  % sess 5
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_24_P2',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_24_P2';  % sess 6
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_25_P2',  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD4d_2023_02_25_P2';  % sess 7
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_12_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_12_P1';% sess 9
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_13_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_13_P2';% sess 10
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_21_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD13d_2024_11_21_P2';% sess 10
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_24_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_24_P2';% sess 11
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_25_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_25_P1';% sess 12
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_26_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_26_P1';% sess 13
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_27_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD15d_2024_11_27_P2';% sess 13
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_17_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_17_P2';% sess 14
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_18_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_18_P2';% sess 15
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_19_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_19_P2';% sess 16
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_20_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_20_P1';% sess 17
%  'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_21_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD22d_2025_06_21_P1';% sess 17
% 'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_17_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_17_P2';% sess 18
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_18_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_18_P1';% sess 18
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_19_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_19_P1';% sess 19
%   'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_20_P2','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_20_P2';% sess 20
% 'C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_21_P1','C:\Users\LabTech\Documents\Cortical Disengagement HMM Results\Results_Final\TD23d_2025_06_20_P1';% sess 20
% };

dataDirs = {
  fullfile(hmmRoot, 'TD1d_2023_02_21_P1'),  fullfile(hmmRoot, 'TD1d_2023_02_21_P2');
  fullfile(hmmRoot, 'TD1d_2023_02_22_P1'),  fullfile(hmmRoot, 'TD1d_2023_02_22_P2');
  fullfile(hmmRoot, 'TD1d_2023_02_23_P1'),  fullfile(hmmRoot, 'TD1d_2023_02_23_P2');
  fullfile(hmmRoot, 'TD1d_2023_02_24_P1'),  fullfile(hmmRoot, 'TD1d_2023_02_24_P2');
  fullfile(hmmRoot, 'TD4d_2023_02_21_P2'),  fullfile(hmmRoot, 'TD4d_2023_02_21_P2');
  fullfile(hmmRoot, 'TD4d_2023_02_24_P2'),  fullfile(hmmRoot, 'TD4d_2023_02_24_P2');
  fullfile(hmmRoot, 'TD4d_2023_02_25_P2'),  fullfile(hmmRoot, 'TD4d_2023_02_25_P2');
  fullfile(hmmRoot, 'TD4d_2023_03_19_P1'),  fullfile(hmmRoot, 'TD4d_2023_03_19_P2');
  fullfile(hmmRoot, 'TD13d_2024_11_12_P1'), fullfile(hmmRoot, 'TD13d_2024_11_12_P1');
  fullfile(hmmRoot, 'TD13d_2024_11_13_P2'), fullfile(hmmRoot, 'TD13d_2024_11_13_P2');
  fullfile(hmmRoot, 'TD13d_2024_11_21_P2'), fullfile(hmmRoot, 'TD13d_2024_11_21_P2');
  fullfile(hmmRoot, 'TD15d_2024_11_24_P2'), fullfile(hmmRoot, 'TD15d_2024_11_24_P2');
  fullfile(hmmRoot, 'TD15d_2024_11_25_P1'), fullfile(hmmRoot, 'TD15d_2024_11_25_P1');
  fullfile(hmmRoot, 'TD15d_2024_11_26_P1'), fullfile(hmmRoot, 'TD15d_2024_11_26_P1');
  fullfile(hmmRoot, 'TD15d_2024_11_27_P2'), fullfile(hmmRoot, 'TD15d_2024_11_27_P2');
  fullfile(hmmRoot, 'TD22d_2025_06_17_P2'), fullfile(hmmRoot, 'TD22d_2025_06_17_P2');
  fullfile(hmmRoot, 'TD22d_2025_06_18_P1'), fullfile(hmmRoot, 'TD22d_2025_06_18_P2');
  fullfile(hmmRoot, 'TD22d_2025_06_19_P1'), fullfile(hmmRoot, 'TD22d_2025_06_19_P2');
  fullfile(hmmRoot, 'TD22d_2025_06_20_P1'), fullfile(hmmRoot, 'TD22d_2025_06_20_P2');
  fullfile(hmmRoot, 'TD22d_2025_06_21_P1'), fullfile(hmmRoot, 'TD22d_2025_06_21_P1');
  fullfile(hmmRoot, 'TD23d_2025_06_17_P1'), fullfile(hmmRoot, 'TD23d_2025_06_17_P2');
  fullfile(hmmRoot, 'TD23d_2025_06_18_P1'), fullfile(hmmRoot, 'TD23d_2025_06_18_P2');
  fullfile(hmmRoot, 'TD23d_2025_06_19_P1'), fullfile(hmmRoot, 'TD23d_2025_06_19_P2');
  fullfile(hmmRoot, 'TD23d_2025_06_20_P1'), fullfile(hmmRoot, 'TD23d_2025_06_20_P2');
  fullfile(hmmRoot, 'TD23d_2025_06_21_P1'), fullfile(hmmRoot, 'TD23d_2025_06_21_P2');   % ALM now uses the _P2 fit
};

y1_A = [];
y2_A = [];
y3_A = [];
y4_A = [];

%     ;meta14;meta15;meta16];

all_meta = [meta1;meta2;meta3;meta4;meta5;meta6;meta7;meta8;meta9;meta10;meta11;meta12 ...
    ;meta13;meta14;meta15;meta16;meta17;meta18;meta19;meta20;meta21;meta22;meta23;meta24;meta25];

y1_all = [];
y2_all = [];
y3_all = [];
y4_all = [];
allP1 = {};
allP4 = {};
allmoveP1 = [];
allmoveP4 = [];

sess_y = {};
sess_y4 = {};
sess_yhat = {};
sess_yhat4 = {};

params.behav_only = 0;

% Probe maps cross-checked against tongue_r14.m (spec.groupMaps, matched by
% animal and date): all 25 sessions and every M1/ALM probe agree. One change:
% Session 25 (TD23d 2025-06-21) ALM is back on probe 2: the HMM Alignments
% folder now holds a TD23d_2025_06_21_P2 fit, so ALM is aligned to its own
% probe's transitions. HMM panels: M1 21 sessions / 6 animals, ALM 15 sessions /
% 4 animals.
% now uses the maps to skip probes that carry neither region.
m1  = [1 1 0 1  2 0 2 1  1 2 2  2 1 1 2  0 1 2 2 1  2 2 2 1 0];
alm = [2 2 2 2  0 2 0 2  0 0 0  0 0 0 0  2 2 1 1 0  1 1 1 2 2];   % session 25 ALM restored to probe 2 (a _P2 fit now exists)
assert(numel(m1) == numel(all_meta) && numel(alm) == numel(all_meta), ...
    'probe maps have %d / %d entries but all_meta has %d sessions.', numel(m1), numel(alm), numel(all_meta));

% trials whose disengagement shift is NaN are dropped; they are tallied here
% and reported once after the loop instead of one warning per session.
nanDropR1 = 0;  nanDropR4 = 0;
for sessnum = 1:length(all_meta)

clear allTrials L_ctrl R_ctrl L_stim R_stim S21c S21 Length angle obj aa aaa idxHit kin

meta = all_meta(sessnum,1);

params.probe = {meta.probe};
params.cluid = {};

[obj, params, kin] = slimToLegacy(meta, params);

fprintf('Session %d\n', sessnum);

conds2use = [1];
kinfeat   = 'tongue_length';
sessix    = 1;

psthForProj = [];
for c = conds2use
    condtrix = params(sessix).trialid{c};
    condpsth = obj(sessix).trialdat(:,:,condtrix);
    condtrix(end) = [];
end

kinix  = find(strcmp(kin(sessix).featLeg, kinfeat));
Length = kin.dat(:, condtrix, kinix);

[reg1, reg, isSingleProbe] = regionSplit(obj(1), params(1));   % was: params.cluid{1,1} with no guard

allRegions = {reg1 reg};

for aa = 1:2

% Only probes the maps use (same rule as tongue_r14's probesOf). aa = PROBE.
if m1(sessnum) ~= aa && alm(sessnum) ~= aa, continue; end

brainRegion = allRegions{aa};

dataDir = dataDirs{sessnum, aa};

fileBases = {'R1_Trial_Track','R4_Trial_Track','R4_dt','R1_dt'};

HMM = readHMM(dataDir, fileBases);   % was: assignin into the base workspace

r1Trials = table2array(HMM.R1_Trial_Track);
r4Trials = table2array(HMM.R4_Trial_Track);

% HMM dt files are ms from the go cue (confirmed via the bigPlot session log).
% * 0.01/params.dt here and /10 below, which composes to the same thing
% (up to double rounding) but left a 10x intermediate under this name.
r1DtsIdx = round(table2array(HMM.R1_dt) * 0.001 / params.dt);
r4DtsIdx = round(table2array(HMM.R4_dt) * 0.001 / params.dt);
logf('  [hmm] sess %2d probe %d | median transition re: GC  R1 %.3f s  R4 %.3f s\n', ...
    sessnum, aa, median(r1DtsIdx,'omitnan')*params.dt, median(r4DtsIdx,'omitnan')*params.dt);

% SETUP
neurons = brainRegion;
[~, idx_zero] = min(abs(obj.time - 0));
allDat  = obj.trialdat(:,neurons,:);
T       = size(allDat,1);
N       = size(allDat,2);

nColors = 256;
cmap_bw = [linspace(0,   1, nColors)', ...
           linspace(0,   1, nColors)', ...
           linspace(0.741, 1, nColors)'];

trials1  = r1Trials(:);
align1   = r1DtsIdx(:);
dtBins1  = align1;   % already bins (was round(align1/10), see above)

trials4  = r4Trials(:);
align4   = r4DtsIdx(:);
dtBins4  = align4;

bad1 = isnan(dtBins1);
if any(bad1)
    nanDropR1 = nanDropR1 + sum(bad1);
    trials1(bad1) = [];  dtBins1(bad1) = [];
end

bad4 = isnan(dtBins4);
if any(bad4)
    nanDropR4 = nanDropR4 + sum(bad4);
    trials4(bad4) = [];  dtBins4(bad4) = [];
end

aligned1 = doAlign(allDat, trials1, dtBins1);
aligned4 = doAlign(allDat, trials4, dtBins4);

psth1 = nan(T,N);
psth4 = nan(T,N);
for nn = 1:N
    psth1(:,nn) = nanmean(aligned1(:,nn,trials1),3);
    psth4(:,nn) = nanmean(aligned4(:,nn,trials4),3);
end

[T,N] = size(psth4);

win1 = obj.time >= -cfg.engModeWin_s & obj.time < 0;
win2 = obj.time >= 0    & obj.time <= cfg.engModeWin_s;
w = engagementModeWeights(aligned4, trials4, win1, win2, cfg.sdFloor);

mod4 = mean(psth4 .* reshape(w,1,N), 2);
mod1 = mean(psth1 .* reshape(w,1,N), 2);

% SHARED SETUP
binSz    = obj.time(2) - obj.time(1);
nL       = size(Length, 1);
time_L   = obj.time(1) + (0:nL-1) * binSz;
nT1      = numel(trials1);
nT4      = numel(trials4);
T_neural = size(allDat,1);
T_safe   = min(T_neural, nL);

proj1_goCue = nan(T_neural, nT1);
proj4_goCue = nan(T_neural, nT4);
for k = 1:nT1
    tr = trials1(k);
    X  = allDat(:,:,tr);
    proj1_goCue(:,k) = projectMode(X, w, cfg.projMode);
end
for k = 1:nT4
    tr = trials4(k);
    X  = allDat(:,:,tr);
    proj4_goCue(:,k) = projectMode(X, w, cfg.projMode);
end
heatmap_mat  = [proj1_goCue'; proj4_goCue'];
nTrialsTotal = nT1 + nT4;

trialOrder = [trials1(:); trials4(:)];
L_sub  = Length(:, min(trialOrder, size(Length,2)));
L_min  = nanmin(L_sub(:));
L_max  = nanmax(L_sub(:));
if L_max == L_min, L_max = L_min + 1; end
row_scale = 0.85;

% MEAN FIRING RATE PER TRIAL
meanFR_R1 = nan(T_neural, nT1);
meanFR_R4 = nan(T_neural, nT4);
for k = 1:nT1
    tr = trials1(k);
    meanFR_R1(:,k) = nanmean(allDat(:,:,tr), 2);
end
for k = 1:nT4
    tr = trials4(k);
    meanFR_R4(:,k) = nanmean(allDat(:,:,tr), 2);
end
msr_mat = [meanFR_R1'; meanFR_R4'];   % nTrials x T

% LICK EVENT DETECTION (for circle overlays)
withinBout_sec_loc    = 0.25;
betweenBout_sec_loc   = 0.5;
minLickLen_sec_loc    = 0.04;
minBoutLicks_loc      = 1;
preGC_search_sec_loc  = 2.5;
preGC_min_sec_loc     = 0.7;
preGC_silence_sec_loc = 0.7;

binSz_L_loc          = time_L(2) - time_L(1);
withinBout_bins_loc  = round(withinBout_sec_loc    / binSz_L_loc);
betweenBout_bins_loc = round(betweenBout_sec_loc   / binSz_L_loc);
minLickLen_bins_loc  = round(minLickLen_sec_loc     / binSz_L_loc);
silence_bins_loc     = round(preGC_silence_sec_loc  / binSz_L_loc);
[~, i_gc_L_loc]      = min(abs(time_L - 0));
i_pre_start_loc      = find(time_L >= -preGC_search_sec_loc, 1, 'first');
i_pre_end_loc        = i_gc_L_loc - 1;

hm_init_times  = [];  hm_init_rows  = [];
hm_reeng_times = [];  hm_reeng_rows = [];
hm_preGC_times = [];  hm_preGC_rows = [];

all_trials_hm = [trials1(:); trials4(:)];

for k = 1:nTrialsTotal
    tr = all_trials_hm(k);
    if tr > size(Length,2), continue; end

% POST-GC bouts
    sig = Length(:, tr)';
    sig(1 : i_gc_L_loc-1) = NaN;

    contact_idx = find(~isnan(sig) & sig > 0);
    if ~isempty(contact_idx)
        d     = diff(contact_idx);
        sp    = find(d > 1);
        s_all = contact_idx([1,  sp+1]);
        e_all = contact_idx([sp, end]);

        ok    = (e_all - s_all + 1) >= minLickLen_bins_loc;
        s_all = s_all(ok);  e_all = e_all(ok);

        if ~isempty(s_all)
            nRuns  = numel(s_all);
            boutID = ones(1,nRuns); curBout = 1;
            for li = 2:nRuns
                if (s_all(li) - e_all(li-1) - 1) > withinBout_bins_loc
                    curBout = curBout + 1;
                end
                boutID(li) = curBout;
            end
            nRaw = curBout;

            rbs_v = nan(1,nRaw);
            rbS_v = cell(1,nRaw);
            rbe_v_last = nan(1,nRaw);
            for b = 1:nRaw
                mem           = boutID == b;
                rbs_v(b)      = s_all(find(mem,1,'first'));
                rbe_v_last(b) = e_all(find(mem,1,'last'));
                rbS_v{b}      = s_all(mem);
            end

            mS_v = {rbS_v{1}}; mEndAny = rbe_v_last(1); nM = 1;
            for b = 2:nRaw
                if (rbs_v(b) - mEndAny - 1) < betweenBout_bins_loc
                    mS_v{nM} = [mS_v{nM}, rbS_v{b}];
                    mEndAny  = rbe_v_last(b);
                else
                    nM = nM+1;
                    mS_v{nM} = rbS_v{b};
                    mEndAny  = rbe_v_last(b);
                end
            end

            finalBouts_s = {};
            for b = 1:nM
                if numel(mS_v{b}) >= minBoutLicks_loc
                    finalBouts_s{end+1} = mS_v{b};
                end
            end

            for bi = 1:numel(finalBouts_s)
                s_bin   = finalBouts_s{bi}(1);
                t_event = time_L(s_bin);
                if bi == 1
                    hm_init_times(end+1) = t_event;
                    hm_init_rows(end+1)  = k;
                else
                    hm_reeng_times(end+1) = t_event;
                    hm_reeng_rows(end+1)  = k;
                end
            end
        end
    end

% PRE-GC bouts
    sig_pre = Length(:, tr)';
    sig_pre(1 : i_pre_start_loc-1) = NaN;
    sig_pre(i_pre_end_loc+1 : end) = NaN;

    contact_idx_pre = find(~isnan(sig_pre) & sig_pre > 0);
    if isempty(contact_idx_pre), continue; end

    d     = diff(contact_idx_pre);
    sp    = find(d > 1);
    s_all = contact_idx_pre([1,  sp+1]);
    e_all = contact_idx_pre([sp, end]);

    ok    = (e_all - s_all + 1) >= minLickLen_bins_loc;
    s_all = s_all(ok);  e_all = e_all(ok);
    if isempty(s_all), continue; end

    nRuns  = numel(s_all);
    boutID = ones(1,nRuns); curBout = 1;
    for li = 2:nRuns
        if (s_all(li) - e_all(li-1) - 1) > withinBout_bins_loc
            curBout = curBout + 1;
        end
        boutID(li) = curBout;
    end
    nRaw = curBout;

    rbS_v = cell(1,nRaw); rbE_v = cell(1,nRaw);
    rbs_v = nan(1,nRaw);  rbe_v = nan(1,nRaw);
    for b = 1:nRaw
        mem      = boutID == b;
        rbs_v(b) = s_all(find(mem,1,'first'));
        rbe_v(b) = e_all(find(mem,1,'last'));
        rbS_v{b} = s_all(mem);
        rbE_v{b} = e_all(mem);
    end

    mS_v = {rbS_v{1}}; mEndAny = rbe_v(1); nM = 1;
    for b = 2:nRaw
        if (rbs_v(b) - mEndAny - 1) < betweenBout_bins_loc
            mS_v{nM} = [mS_v{nM}, rbS_v{b}];
            mEndAny  = rbe_v(b);
        else
            nM = nM+1;
            mS_v{nM} = rbS_v{b};
            mEndAny  = rbe_v(b);
        end
    end

    finalBouts_pre_s = {};
    for b = 1:nM
        if numel(mS_v{b}) >= minBoutLicks_loc
            finalBouts_pre_s{end+1} = mS_v{b};
        end
    end

    for bi = 1:numel(finalBouts_pre_s)
        s_bin        = finalBouts_pre_s{bi}(1);
        t_lick_re_gc = (s_bin - i_gc_L_loc) * binSz_L_loc;
        if t_lick_re_gc > -preGC_min_sec_loc, continue; end

        guard_start = max(1, s_bin - silence_bins_loc);
        guard_end   = s_bin - 1;
        if guard_end >= guard_start
            guard_sig = Length(:, tr)';
            guard_sig = guard_sig(guard_start:guard_end);
            guard_sig(isnan(guard_sig)) = 0;
            if any(guard_sig > 0), continue; end
        end

        hm_preGC_times(end+1) = time_L(s_bin);
        hm_preGC_rows(end+1)  = k;
    end
end

% FIGURE — Side-by-side: Mode Projection (left) & Mean Spike Rate (right)
clim_proj = prctile(abs(heatmap_mat(:)), 95);
msr_vals  = msr_mat(~isnan(msr_mat(:)));
clim_msr  = prctile(msr_vals, 95);

if cfg.plotPerSession
figure('Color','w','Units','normalized','Position',[.02 .05 .90 .55]);

% helper: draw tongue + circles on an axes
mk_sz = 30;

% ---- LEFT: Mode projection ----
ax_proj = axes('Position',[0.05 0.12 0.38 0.78]);
imagesc(ax_proj, obj.time, 1:nTrialsTotal, heatmap_mat);
hold(ax_proj,'on');
yline(nT1+0.5, 'w-', 'LineWidth', 2);
xline(0,       'w--','LineWidth', 1.5);
caxis(ax_proj, [-clim_proj  clim_proj]);
colormap(ax_proj, cmap_bw);
cb1 = colorbar(ax_proj,'eastoutside');
cb1.Label.String = 'Mode projection';

for k = 1:nT1
    tr = trials1(k);
    if tr > size(Length,2), continue; end
    ln = (Length(:,tr) - L_min) / (L_max - L_min);
    plot(ax_proj, time_L, k + ln*row_scale, 'r-', 'LineWidth',0.5);
end
for k = 1:nT4
    tr = trials4(k);
    if tr > size(Length,2), continue; end
    ln = (Length(:,tr) - L_min) / (L_max - L_min);
    plot(ax_proj, time_L, (nT1+k) + ln*row_scale, 'r-', 'LineWidth',0.5);
end
if ~isempty(hm_init_times)
    scatter(ax_proj, hm_init_times, hm_init_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[1.00 0.50 0.00],'MarkerFaceColor','none','LineWidth',2);
end
if ~isempty(hm_reeng_times)
    scatter(ax_proj, hm_reeng_times, hm_reeng_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[1.00 0.41 0.71],'MarkerFaceColor','none','LineWidth',2);
end
if ~isempty(hm_preGC_times)
    scatter(ax_proj, hm_preGC_times, hm_preGC_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[0.20 0.63 0.17],'MarkerFaceColor','none','LineWidth',2);
end
ylabel(ax_proj,'Trial #'); xlabel(ax_proj,'Time (s)');
ylim(ax_proj,[0.5 nTrialsTotal+0.5]); xlim(ax_proj,[-1 2]);
set(ax_proj,'TickDir','out','FontSize',10);
text(ax_proj,-0.9, nT1/2,     'R1','Color','w','FontSize',12,'FontWeight','bold');
text(ax_proj,-0.9, nT1+nT4/2, 'R4','Color','w','FontSize',12,'FontWeight','bold');
title(ax_proj, sprintf('Mode proj.  |  %s  %s  probe%d', ...
    num2str(obj.pth.dt), obj.pth.anm, aa), 'FontSize',10);
box off;

% ---- RIGHT: Mean Spike Rate ----
ax_msr = axes('Position',[0.52 0.12 0.38 0.78]);
imagesc(ax_msr, obj.time, 1:nTrialsTotal, msr_mat);
hold(ax_msr,'on');
yline(nT1+0.5, 'w-', 'LineWidth', 2);
xline(0,       'w--','LineWidth', 1.5);
caxis(ax_msr, [0  clim_msr]);
colormap(ax_msr, parula);
cb2 = colorbar(ax_msr,'eastoutside');
cb2.Label.String = 'Mean spike rate (spk/s)';

for k = 1:nT1
    tr = trials1(k);
    if tr > size(Length,2), continue; end
    ln = (Length(:,tr) - L_min) / (L_max - L_min);
    plot(ax_msr, time_L, k + ln*row_scale, 'r-', 'LineWidth',0.5);
end
for k = 1:nT4
    tr = trials4(k);
    if tr > size(Length,2), continue; end
    ln = (Length(:,tr) - L_min) / (L_max - L_min);
    plot(ax_msr, time_L, (nT1+k) + ln*row_scale, 'r-', 'LineWidth',0.5);
end
if ~isempty(hm_init_times)
    scatter(ax_msr, hm_init_times, hm_init_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[1.00 0.50 0.00],'MarkerFaceColor','none','LineWidth',2);
end
if ~isempty(hm_reeng_times)
    scatter(ax_msr, hm_reeng_times, hm_reeng_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[1.00 0.41 0.71],'MarkerFaceColor','none','LineWidth',2);
end
if ~isempty(hm_preGC_times)
    scatter(ax_msr, hm_preGC_times, hm_preGC_rows, mk_sz, 'o', ...
        'MarkerEdgeColor',[0.20 0.63 0.17],'MarkerFaceColor','none','LineWidth',2);
end
ylabel(ax_msr,'Trial #'); xlabel(ax_msr,'Time (s)');
ylim(ax_msr,[0.5 nTrialsTotal+0.5]); xlim(ax_msr,[-1 2]);
set(ax_msr,'TickDir','out','FontSize',10);
text(ax_msr,-0.9, nT1/2,     'R1','Color','w','FontSize',12,'FontWeight','bold');
text(ax_msr,-0.9, nT1+nT4/2, 'R4','Color','w','FontSize',12,'FontWeight','bold');
title(ax_msr, sprintf('Mean spike rate  |  %s  %s  probe%d', ...
    num2str(obj.pth.dt), obj.pth.anm, aa), 'FontSize',10);
box off;
end   % cfg.plotPerSession

% SAVE
All_R1Modes{sessnum,aa} = mod1;
All_R4Modes{sessnum,aa} = mod4;

AllProj_R1{sessnum, aa} = proj1_goCue;
AllProj_R4{sessnum, aa} = proj4_goCue;

len_R1 = nan(nL, nT1);
len_R4 = nan(nL, nT4);
for k = 1:nT1
    tr = trials1(k);
    if tr <= size(Length,2), len_R1(:,k) = Length(:,tr); end
end
for k = 1:nT4
    tr = trials4(k);
    if tr <= size(Length,2), len_R4(:,k) = Length(:,tr); end
end

AllLen_R1{sessnum, aa} = len_R1;
AllLen_R4{sessnum, aa} = len_R4;

AllSessID_R1{sessnum, aa} = repmat(sessnum, 1, nT1);
AllSessID_R4{sessnum, aa} = repmat(sessnum, 1, nT4);

AllTime{sessnum, aa}  = obj.time;
AllTimeL{sessnum, aa} = time_L;

AllMeanFR_R1{sessnum, aa} = meanFR_R1;
AllMeanFR_R4{sessnum, aa} = meanFR_R4;

end   % aa loop

end   % sessnum loop


% GLOBAL PARAMETERS
withinBout_sec     = 0.25;
betweenBout_sec    = 0.5;
minLickLen_sec     = 0.025;
minBoutLicks       = 1;
preBins_sec        = 1.0;
postBins_sec       = 1.0;
preGC_search_sec   = 2.5;
preGC_min_sec      = 0.7;
preGC_silence_sec  = 0.7;
bl_start_sec       = -0.6;
bl_end_sec         = -0.2;
xlim_hmap          = [-0.25  0.15];
xlim_line          = [-0.25  0.25];
pct_clim           = 95;
nColors            = 256;

% SETUP
m1  = m1(:);
alm = alm(:);
nSess = numel(m1);

t_ax = []; t_L = [];
for s = 1:nSess
    for aa = 1:2
        if size(AllProj_R1,1)>=s && ~isempty(AllProj_R1{s,aa})
            t_ax = AllTime{s,aa};
            t_L  = AllTimeL{s,aa};
            break;
        end
    end
    if ~isempty(t_ax), break; end
end

binSz_L            = t_L(2) - t_L(1);
withinBout_bins_L  = round(withinBout_sec    / binSz_L);
betweenBout_bins_L = round(betweenBout_sec   / binSz_L);
minLickLen_L       = round(minLickLen_sec     / binSz_L);
preBins            = round(preBins_sec        / binSz_L);
postBins           = round(postBins_sec       / binSz_L);
i_gc_L             = find(t_L >= 0,           1, 'first');
i_pre_start        = find(t_L >= -preGC_search_sec, 1, 'first');
i_pre_end          = i_gc_L - 1;
silence_bins       = round(preGC_silence_sec  / binSz_L);
snip_len           = preBins + postBins + 1;
t_common           = (-preBins : postBins) * binSz_L;
bl_mask            = t_common >= bl_start_sec & t_common <= bl_end_sec;

% CONCATENATE — M1 and ALM, with per-session 0-1 FR normalization
regionDefs = {m1, 'M1'; alm, 'ALM'};

init_cells     = cell(1,2);
reeng_cells    = cell(1,2);
pre_init_cells = cell(1,2);

for rr = 1:2
    map   = regionDefs{rr,1};
    rname = regionDefs{rr,2};

    idx_sess = find(map ~= 0);

    all_init     = {};
    all_reeng    = {};
    all_pre_init = {};

    for ki = 1:numel(idx_sess)
        sess  = idx_sess(ki);
        probe = map(sess);

        if size(AllProj_R1,1) < sess || isempty(AllProj_R1{sess,probe}), continue; end

% Pull per-session data
        t_L_s  = AllTimeL{sess, probe};
        t_ax_s = AllTime{sess,  probe};
        if isempty(t_L_s) || isempty(t_ax_s), continue; end

        binSz_L_s          = t_L_s(2) - t_L_s(1);
        withinBout_bins_s  = round(withinBout_sec    / binSz_L_s);
        betweenBout_bins_s = round(betweenBout_sec   / binSz_L_s);
        minLickLen_bins_s  = round(minLickLen_sec     / binSz_L_s);
        silence_bins_s     = round(preGC_silence_sec  / binSz_L_s);
        [~, i_gc_L_s]      = min(abs(t_L_s - 0));
        i_pre_start_s      = find(t_L_s >= -preGC_search_sec, 1, 'first');
        i_pre_end_s        = i_gc_L_s - 1;

        binSz_ax_s   = t_ax_s(2) - t_ax_s(1);
        T_neural_s   = numel(t_ax_s);
        nL_s         = numel(t_L_s);
        preBins_s    = round(preBins_sec  / binSz_ax_s);
        postBins_s   = round(postBins_sec / binSz_ax_s);
        preBins_L_s  = round(preBins_sec  / binSz_L_s);
        postBins_L_s = round(postBins_sec / binSz_L_s);

% proj: T x nTrials -> transpose to nTrials x T
        p_R1 = AllProj_R1{sess, probe}';
        p_R4 = AllProj_R4{sess, probe}';
        l_R1 = AllLen_R1{sess, probe}';
        l_R4 = AllLen_R4{sess, probe}';
        f_R1 = AllMeanFR_R1{sess, probe}';
        f_R4 = AllMeanFR_R4{sess, probe}';

% Per-session 0-1 normalization of FR
        fr_all = [f_R1(:); f_R4(:)];
        fr_min_s = nanmin(fr_all);
        fr_max_s = nanmax(fr_all);
        if isempty(fr_min_s) || fr_max_s == fr_min_s, fr_max_s = fr_min_s + 1; end
        f_R1 = (f_R1 - fr_min_s) / (fr_max_s - fr_min_s);
        f_R4 = (f_R4 - fr_min_s) / (fr_max_s - fr_min_s);

        for trType = 1:2
            if trType == 1
                proj = p_R1; len = l_R1; fr = f_R1;
            else
                proj = p_R4; len = l_R4; fr = f_R4;
            end
            if isempty(proj) || isempty(len), continue; end

            nTr = min(size(proj,1), size(len,1));

            for tr = 1:nTr

% POST-GC bouts
                sig = len(tr, :);
                sig(1 : i_gc_L_s-1) = NaN;

                contact_idx = find(~isnan(sig) & sig > 0);
                if ~isempty(contact_idx)
                    d_diff = diff(contact_idx);
                    sp     = find(d_diff > 1);
                    s_all  = contact_idx([1,  sp+1]);
                    e_all  = contact_idx([sp, end]);

                    ok    = (e_all - s_all + 1) >= minLickLen_bins_s;
                    s_all = s_all(ok);  e_all = e_all(ok);

                    if ~isempty(s_all)
                        nRuns  = numel(s_all);
                        boutID = ones(1,nRuns); curBout = 1;
                        for li = 2:nRuns
                            if (s_all(li) - e_all(li-1) - 1) > withinBout_bins_s
                                curBout = curBout + 1;
                            end
                            boutID(li) = curBout;
                        end
                        nRaw = curBout;

                        rbs_v = nan(1,nRaw); rbe_v = nan(1,nRaw);
                        rbS_v = cell(1,nRaw); rbE_v = cell(1,nRaw);
                        for b = 1:nRaw
                            mem      = boutID == b;
                            rbs_v(b) = s_all(find(mem,1,'first'));
                            rbe_v(b) = e_all(find(mem,1,'last'));
                            rbS_v{b} = s_all(mem);
                            rbE_v{b} = e_all(mem);
                        end

                        mS_v = {rbS_v{1}}; mE_v = {rbE_v{1}}; mEndAny = rbe_v(1); nM = 1;
                        for b = 2:nRaw
                            if (rbs_v(b) - mEndAny - 1) < betweenBout_bins_s
                                mS_v{nM} = [mS_v{nM}, rbS_v{b}];
                                mE_v{nM} = [mE_v{nM}, rbE_v{b}];
                                mEndAny  = rbe_v(b);
                            else
                                nM = nM+1;
                                mS_v{nM} = rbS_v{b};
                                mE_v{nM} = rbE_v{b};
                                mEndAny  = rbe_v(b);
                            end
                        end

                        finalBouts_s = {}; finalBouts_e = {};
                        for b = 1:nM
                            if numel(mS_v{b}) >= minBoutLicks
                                finalBouts_s{end+1} = mS_v{b};
                                finalBouts_e{end+1} = mE_v{b};
                            end
                        end

                        for bi = 1:numel(finalBouts_s)
                            s_kin   = finalBouts_s{bi}(1);
                            e_kin   = finalBouts_e{bi}(1);
                            lickDur = (e_kin - s_kin + 1) * binSz_L_s;

                            t_lick = t_L_s(s_kin);
                            s_ax   = round((t_lick - t_ax_s(1)) / binSz_ax_s) + 1;

% NaN-padded proj snippet (on kinematic axis)
                            snip_proj = nan(1, snip_len);
                            i1r = s_kin - preBins_L_s;  i2r = s_kin + postBins_L_s;
                            i1  = max(1, i1r);           i2  = min(nL_s, i2r);
                            if i2 >= i1
                                o1 = i1 - i1r + 1;  o2 = o1 + (i2 - i1);
                                snip_proj(o1:o2) = proj(tr, i1:i2);
                            end

% NaN-padded FR snippet (on neural axis)
                            snip_fr = nan(1, snip_len);
                            i1r_ax = s_ax - preBins_s;  i2r_ax = s_ax + postBins_s;
                            i1_ax  = max(1, i1r_ax);     i2_ax  = min(T_neural_s, i2r_ax);
                            if i2_ax >= i1_ax
                                o1 = i1_ax - i1r_ax + 1;  o2 = o1 + (i2_ax - i1_ax);
                                snip_fr(o1:o2) = fr(tr, i1_ax:i2_ax);
                            end

% NaN-padded len snippet
                            snip_len_vec = nan(1, snip_len);
                            if i2 >= i1
                                snip_len_vec(o1:o2) = len(tr, i1:i2);
                            end
                            snip_len_vec(isnan(snip_len_vec)) = 0;

                            t_gc_rel = (i_gc_L_s - s_kin) * binSz_L_s;
                            entry = {snip_proj, t_common, lickDur, snip_fr, snip_len_vec, s_kin, t_gc_rel};

                            if bi == 1
                                all_init{end+1} = entry;
                            else
                                all_reeng{end+1} = entry;
                            end
                        end
                    end
                end

% PRE-GC bouts
                sig_pre = len(tr, :);
                sig_pre(1 : i_pre_start_s-1) = NaN;
                sig_pre(i_pre_end_s+1 : end)  = NaN;

                contact_idx_pre = find(~isnan(sig_pre) & sig_pre > 0);
                if isempty(contact_idx_pre), continue; end

                d_diff = diff(contact_idx_pre);
                sp     = find(d_diff > 1);
                s_all  = contact_idx_pre([1,  sp+1]);
                e_all  = contact_idx_pre([sp, end]);

                ok    = (e_all - s_all + 1) >= minLickLen_bins_s;
                s_all = s_all(ok);  e_all = e_all(ok);
                if isempty(s_all), continue; end

                nRuns  = numel(s_all);
                boutID = ones(1,nRuns); curBout = 1;
                for li = 2:nRuns
                    if (s_all(li) - e_all(li-1) - 1) > withinBout_bins_s
                        curBout = curBout + 1;
                    end
                    boutID(li) = curBout;
                end
                nRaw = curBout;

                rbs_v = nan(1,nRaw); rbe_v = nan(1,nRaw);
                rbS_v = cell(1,nRaw); rbE_v = cell(1,nRaw);
                for b = 1:nRaw
                    mem      = boutID == b;
                    rbs_v(b) = s_all(find(mem,1,'first'));
                    rbe_v(b) = e_all(find(mem,1,'last'));
                    rbS_v{b} = s_all(mem);
                    rbE_v{b} = e_all(mem);
                end

                mS_v = {rbS_v{1}}; mE_v = {rbE_v{1}}; mEndAny = rbe_v(1); nM = 1;
                for b = 2:nRaw
                    if (rbs_v(b) - mEndAny - 1) < betweenBout_bins_s
                        mS_v{nM} = [mS_v{nM}, rbS_v{b}];
                        mE_v{nM} = [mE_v{nM}, rbE_v{b}];
                        mEndAny  = rbe_v(b);
                    else
                        nM = nM+1;
                        mS_v{nM} = rbS_v{b};
                        mE_v{nM} = rbE_v{b};
                        mEndAny  = rbe_v(b);
                    end
                end

                finalBouts_pre_s = {}; finalBouts_pre_e = {};
                for b = 1:nM
                    if numel(mS_v{b}) >= minBoutLicks
                        finalBouts_pre_s{end+1} = mS_v{b};
                        finalBouts_pre_e{end+1} = mE_v{b};
                    end
                end
                if isempty(finalBouts_pre_s), continue; end

                for bi = 1:numel(finalBouts_pre_s)
                    s_kin = finalBouts_pre_s{bi}(1);
                    e_kin = finalBouts_pre_e{bi}(1);

% Guard 1
                    t_lick_re_gc = (s_kin - i_gc_L_s) * binSz_L_s;
                    if t_lick_re_gc > -preGC_min_sec, continue; end

% Guard 2 — silence before lick
                    guard_start = max(1, s_kin - silence_bins_s);
                    guard_end   = s_kin - 1;
                    if guard_end >= guard_start
                        guard_sig = len(tr, guard_start:guard_end);
                        guard_sig(isnan(guard_sig)) = 0;
                        if any(guard_sig > 0), continue; end
                    end

                    lickDur = (e_kin - s_kin + 1) * binSz_L_s;

                    t_lick = t_L_s(s_kin);
                    s_ax   = round((t_lick - t_ax_s(1)) / binSz_ax_s) + 1;

% Guard 3 — at least 0.1s of recording before lick
                    min_pre_bins = round(0.1 / binSz_ax_s);
                    if s_ax - min_pre_bins < 1, continue; end

% NaN-padded proj snippet
                    snip_proj = nan(1, snip_len);
                    i1r = s_kin - preBins_L_s;  i2r = s_kin + postBins_L_s;
                    i1  = max(1, i1r);           i2  = min(nL_s, i2r);
                    if i2 >= i1
                        o1 = i1 - i1r + 1;  o2 = o1 + (i2 - i1);
                        snip_proj(o1:o2) = proj(tr, i1:i2);
                    end

% NaN-padded FR snippet
                    snip_fr = nan(1, snip_len);
                    i1r_ax = s_ax - preBins_s;  i2r_ax = s_ax + postBins_s;
                    i1_ax  = max(1, i1r_ax);     i2_ax  = min(T_neural_s, i2r_ax);
                    if i2_ax >= i1_ax
                        o1_ax = i1_ax - i1r_ax + 1;  o2_ax = o1_ax + (i2_ax - i1_ax);
                        snip_fr(o1_ax:o2_ax) = fr(tr, i1_ax:i2_ax);
                    end

% NaN-padded len snippet
                    snip_len_vec = nan(1, snip_len);
                    if i2 >= i1
                        snip_len_vec(o1:o2) = len(tr, i1:i2);
                    end
                    snip_len_vec(isnan(snip_len_vec)) = 0;

                    t_gc_rel = (i_gc_L_s - s_kin) * binSz_L_s;
                    entry = {snip_proj, t_common, lickDur, snip_fr, snip_len_vec, s_kin, t_gc_rel};

                    all_pre_init{end+1} = entry;
                end

            end   % tr
        end   % trType
    end   % ki (sessions)

    init_cells{rr}     = all_init;
    reeng_cells{rr}    = all_reeng;
    pre_init_cells{rr} = all_pre_init;
end

% Pool M1 + ALM
init_all     = [init_cells{1},     init_cells{2}];
reeng_all    = [reeng_cells{1},    reeng_cells{2}];
pre_init_all = [pre_init_cells{1}, pre_init_cells{2}];

if nanDropR1 > 0 || nanDropR4 > 0
    fprintf('Dropped %d R1 and %d R4 trials with a NaN disengagement shift.\n', ...
        nanDropR1, nanDropR4);
end

fprintf('Initiation      — M1: %d,  ALM: %d\n', numel(init_cells{1}),     numel(init_cells{2}));
fprintf('Re-engagement   — M1: %d,  ALM: %d\n', numel(reeng_cells{1}),    numel(reeng_cells{2}));
fprintf('Pre-go-cue init — M1: %d,  ALM: %d\n', numel(pre_init_cells{1}), numel(pre_init_cells{2}));

evData = { init_all,     'Initiation  (1st lick after go cue)'; ...
           reeng_all,    'Re-engagement  (1st lick, bouts 2-N)'; ...
           pre_init_all, 'Pre-go-cue initiation'};

% BUILD MATRICES — MODE PROJECTION
rng(42);

hmaps           = cell(1,3);
mats_line       = cell(1,3);
mats_len        = cell(1,3);
durs_all_sorted = cell(1,3);
ns              = zeros(1,3);

for ev = 1:3
    cells  = evData{ev,1};
    n      = numel(cells);
    ns(ev) = n;
    if n == 0, continue; end

    mat_h   = nan(n, snip_len);
    mat_l   = nan(n, snip_len);
    mat_len = nan(n, snip_len);
    durs    = nan(1,n);

    for k = 1:n
        seg     = cells{k}{1};
        lickDur = cells{k}{3};
        seg_len = cells{k}{5};

        mat_l(k,:)   = seg;
        mat_len(k,:) = seg_len;

        seg_h = seg;
        seg_h(t_common > lickDur) = NaN;
        mat_h(k,:) = seg_h;

        durs(k) = lickDur;
    end

    mats_line{ev} = mat_l;

% Randomize then sort by duration
    rand_order   = randperm(n);
    mat_h_rand   = mat_h(rand_order, :);
    mat_len_rand = mat_len(rand_order, :);
    durs_rand    = durs(rand_order);

    [durs_sorted, si]   = sort(durs_rand, 'ascend');
    hmaps{ev}           = mat_h_rand(si, :);
    durs_all_sorted{ev} = durs_sorted;
    mats_len{ev}        = mat_len_rand(si, :);
end

% Clim: non-NaN only
all_proj_vals = [];
for ev = 1:3
    if ~isempty(hmaps{ev})
        v = hmaps{ev}(:);
        all_proj_vals = [all_proj_vals; v(~isnan(v))];
    end
end
clim_lo_proj = prctile(all_proj_vals, 5);
clim_hi_proj = prctile(all_proj_vals, 95);
if clim_lo_proj == clim_hi_proj, clim_hi_proj = clim_lo_proj + 1; end

% The projection heatmaps below are drawn on a 0-1 scale instead of in
% raw projection units. The pooled 5th-95th percentile, which is the range the
% color axis already spanned, maps to 0 and 1, so the picture is unchanged and
% only the numbers on the colorbar are. Values outside that range are clipped
% to 0 or 1 so nothing falls off the end of the color map. NaNs stay NaN and
% are still drawn with a sentinel below the axis.
hmaps_norm = cell(1,3);
for ev = 1:3
    if isempty(hmaps{ev}), continue; end
    z = (hmaps{ev} - clim_lo_proj) ./ (clim_hi_proj - clim_lo_proj);
    z(z < 0) = 0;
    z(z > 1) = 1;
    hmaps_norm{ev} = z;
end

% Tongue length normalization
len_all_vals = [];
for ev = 1:3
    if ~isempty(mats_len{ev}), len_all_vals = [len_all_vals; mats_len{ev}(:)]; end
end
len_gmin = nanmin(len_all_vals);  len_gmax = nanmax(len_all_vals);
if len_gmax == len_gmin, len_gmax = len_gmin + 1; end
mats_len_norm = cell(1,3);
for ev = 1:3
    if ~isempty(mats_len{ev})
        mats_len_norm{ev} = (mats_len{ev} - len_gmin) / (len_gmax - len_gmin);
    end
end

% BUILD MATRICES — MEAN SPIKE RATE
rng(42);

hmaps_fr           = cell(1,3);
mats_line_fr       = cell(1,3);
mats_len_fr        = cell(1,3);
durs_all_sorted_fr = cell(1,3);
ns_fr              = zeros(1,3);

for ev = 1:3
    cells     = evData{ev,1};
    n         = numel(cells);
    ns_fr(ev) = n;
    if n == 0, continue; end

    mat_h_fr  = nan(n, snip_len);
    mat_l_fr  = nan(n, snip_len);
    mat_len_fr = nan(n, snip_len);
    durs_fr   = nan(1,n);

    for k = 1:n
        seg_fr  = cells{k}{4};
        lickDur = cells{k}{3};
        seg_len = cells{k}{5};

        mat_l_fr(k,:)   = seg_fr;
        mat_len_fr(k,:) = seg_len;

        seg_h_fr = seg_fr;
        seg_h_fr(t_common > lickDur) = NaN;
        mat_h_fr(k,:) = seg_h_fr;

        durs_fr(k) = lickDur;
    end

    mats_line_fr{ev} = mat_l_fr;

% Randomize then sort by duration
    rand_order_fr   = randperm(n);
    mat_h_fr_rand   = mat_h_fr(rand_order_fr, :);
    mat_len_fr_rand = mat_len_fr(rand_order_fr, :);
    durs_fr_rand    = durs_fr(rand_order_fr);

    [durs_sorted_fr, si_fr]   = sort(durs_fr_rand, 'ascend');
    hmaps_fr{ev}              = mat_h_fr_rand(si_fr, :);
    durs_all_sorted_fr{ev}    = durs_sorted_fr;
    mats_len_fr{ev}           = mat_len_fr_rand(si_fr, :);
end

% Clim: non-NaN only, 5th-95th percentile
all_fr_vals = [];
for ev = 1:3
    if ~isempty(hmaps_fr{ev})
        v = hmaps_fr{ev}(:);
        all_fr_vals = [all_fr_vals; v(~isnan(v))];
    end
end
clim_lo_fr = prctile(all_fr_vals, 5);
clim_hi_fr = prctile(all_fr_vals, 95);
if clim_lo_fr == clim_hi_fr, clim_hi_fr = clim_lo_fr + 1; end

% Tongue normalization for FR plots
len_fr_vals = [];
for ev = 1:3
    if ~isempty(mats_len_fr{ev}), len_fr_vals = [len_fr_vals; mats_len_fr{ev}(:)]; end
end
len_gmin_fr = nanmin(len_fr_vals);  len_gmax_fr = nanmax(len_fr_vals);
if len_gmax_fr == len_gmin_fr, len_gmax_fr = len_gmin_fr + 1; end
mats_len_norm_fr = cell(1,3);
for ev = 1:3
    if ~isempty(mats_len_fr{ev})
        mats_len_norm_fr{ev} = (mats_len_fr{ev} - len_gmin_fr) / (len_gmax_fr - len_gmin_fr);
    end
end

% SHARED STYLE
col_init     = [1.00 0.50 0.00];
col_reeng    = [1.00 0.41 0.71];
col_pre_init = [0.47 0.67 0.19];
col_tick     = [0.85 0.10 0.10];
tickH        = 0.4;
tickLW       = 3;
lw_main      = 3;
cols         = {col_init, col_reeng, col_pre_init};
leg_names    = {'Go cue init', 'Re-engagement', 'Pre-GC init'};
ev_titles    = {'Initiation  (1st lick after go cue)', ...
                'Re-engagement  (1st lick, bouts 2-N)', ...
                sprintf('Pre-go-cue init  (>=%.1fs before GC)', preGC_min_sec)};

% FIGURES 1-3 — MODE PROJECTION HEATMAPS (subplot 1x3)
figure('Color','w','Units','normalized','Position',[0.02 0.10 0.95 0.80]);
sgtitle('Mode projection  —  sorted by lick duration', 'FontSize',13,'FontWeight','bold');

for ev = 1:3
    ax = subplot(1, 3, ev);
    n  = ns(ev);

    if n == 0
        title(ax, sprintf('%s\nn=0', ev_titles{ev}), 'FontSize',10);
        axis(ax,'off'); continue;
    end

    hmap        = hmaps_norm{ev};   % 0-1 scale, not raw projection units
    durs_sorted = durs_all_sorted{ev};

% Replace NaN with sentinel below clim
    hmap_plot = hmap;
    hmap_plot(isnan(hmap_plot)) = -1;   % below the 0-1 color axis

    imagesc(ax, t_common, 1:n, hmap_plot);
    hold(ax,'on');
    for k = 1:n
        plot(ax, [0 0],   [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
        tE = durs_sorted(k);
        plot(ax, [tE tE], [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
    end
    caxis(ax, [0  1]);   % 0-1 scale
    colormap(ax, parula);
    cb = colorbar(ax,'eastoutside');
    cb.Ticks      = [0 0.5 1];
    cb.TickLabels = {'0','0.5','1'};
    cb.TickLength = 0;
    cb.Label.String = 'Mode proj. (norm.)';
    xlabel(ax,'Time re. lick start (s)');
    if ev == 1, ylabel(ax,'Event # (sorted by lick duration)'); end
    xlim(ax, xlim_hmap); ylim(ax, [0.5 n+0.5]);
    set(ax,'YDir','normal','TickDir','out','FontSize',10);
    title(ax, sprintf('%s\nn=%d (M1+ALM, R1+R4)', ev_titles{ev}, n), 'FontSize',10);
    box(ax,'off');
end

% FIGURE 4 — Mean mode projection, all 3 conditions
figure('Color','w','Units','normalized','Position',[.25 .10 .45 .55]);
ax_proj = axes();
hold(ax_proj,'on');

for ev = 1:3
    if isempty(mats_line{ev}), continue; end
    mat = mats_line{ev};
    n   = ns(ev);
    mu  = nanmean(mat, 1);
    n_valid = sum(~isnan(mat), 1);
    sem = nanstd(mat, 0, 1) ./ sqrt(max(n_valid, 1));
    ci  = 1.96 * sem;
    valid = n_valid >= 2;
    t_p = t_common(valid);  mu_p = mu(valid);  ci_p = ci(valid);
    in_xl = t_p >= xlim_line(1) & t_p <= xlim_line(2);
    t_p = t_p(in_xl);  mu_p = mu_p(in_xl);  ci_p = ci_p(in_xl);

    if numel(t_p) > 1
        fill(ax_proj, [t_p fliplr(t_p)], [mu_p+ci_p fliplr(mu_p-ci_p)], cols{ev}, ...
             'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax_proj, t_p, mu_p, '-', 'Color',cols{ev}, 'LineWidth',lw_main, ...
         'DisplayName', sprintf('%s (n=%d)', leg_names{ev}, n));
end

yyaxis(ax_proj,'right');
for ev = 1:3
    if isempty(mats_len_norm{ev}), continue; end
    ln_mu = nanmean(mats_len_norm{ev}, 1);
    ln_mu(isnan(ln_mu)) = 0;
    in_xl = t_common >= xlim_line(1) & t_common <= xlim_line(2);
    plot(ax_proj, t_common(in_xl), ln_mu(in_xl), ':', 'Color',cols{ev}, 'LineWidth',2, ...
         'HandleVisibility','off');
end
ax_proj.YAxis(2).Color = [0 0 0];
ylim(ax_proj, [0 1]);
ylabel(ax_proj,'Tongue length (norm)');
yyaxis(ax_proj,'left');

xline(ax_proj, 0, 'k--', 'LineWidth',1,'HandleVisibility','off');
xlim(ax_proj, xlim_line);
xlabel(ax_proj,'Time re. lick start (s)');
ylabel(ax_proj,'Mode proj.');
title(ax_proj,'Mode projection','FontSize',11);
legend(ax_proj, 'Location','northwest','FontSize',9);
box(ax_proj,'off');
set(ax_proj,'TickDir','out','FontSize',11);

% FIGURES 5-7 — MEAN SPIKE RATE HEATMAPS (subplot 1x3)
figure('Color','w','Units','normalized','Position',[0.02 0.10 0.95 0.80]);
sgtitle('Mean spike rate (0-1 norm)  —  sorted by lick duration', 'FontSize',13,'FontWeight','bold');

for ev = 1:3
    ax = subplot(1, 3, ev);
    n  = ns_fr(ev);

    if n == 0
        title(ax, sprintf('%s\nn=0', ev_titles{ev}), 'FontSize',10);
        axis(ax,'off'); continue;
    end

    hmap        = hmaps_fr{ev};
    durs_sorted = durs_all_sorted_fr{ev};

% Replace NaN with sentinel below clim
    hmap_plot = hmap;
    hmap_plot(isnan(hmap_plot)) = clim_lo_fr - 1;

    imagesc(ax, t_common, 1:n, hmap_plot);
    hold(ax,'on');
    for k = 1:n
        plot(ax, [0 0],   [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
        tE = durs_sorted(k);
        plot(ax, [tE tE], [k-tickH k+tickH], '-', 'Color',col_tick, 'LineWidth',tickLW);
    end
    caxis(ax, [clim_lo_fr  clim_hi_fr]);
    colormap(ax, parula);
    cb = colorbar(ax,'eastoutside');
    cb.Ticks      = [clim_lo_fr  clim_hi_fr];
    cb.TickLabels = {sprintf('%.2f',clim_lo_fr), sprintf('%.2f',clim_hi_fr)};
    cb.TickLength = 0;
    cb.Label.String = 'Mean spike rate (0-1 norm)';
    xlabel(ax,'Time re. lick start (s)');
    if ev == 1, ylabel(ax,'Event # (sorted by lick duration)'); end
    xlim(ax, xlim_hmap); ylim(ax, [0.5 n+0.5]);
    set(ax,'YDir','normal','TickDir','out','FontSize',10);
    title(ax, sprintf('%s\nn=%d (M1+ALM, R1+R4)', ev_titles{ev}, n), 'FontSize',10);
    box(ax,'off');
end

% FIGURE 8 — Mean spike rate + tongue length, all 3 conditions
figure('Color','w','Units','normalized','Position',[.25 .10 .45 .55]);
ax_fr = axes();
hold(ax_fr,'on');

for ev = 1:3
    if isempty(mats_line_fr{ev}), continue; end
    mat = mats_line_fr{ev};
    n   = ns_fr(ev);
    mu  = nanmean(mat, 1);
    n_valid = sum(~isnan(mat), 1);
    sem = nanstd(mat, 0, 1) ./ sqrt(max(n_valid, 1));
    ci  = 1.96 * sem;
    valid = n_valid >= 2;
    t_p = t_common(valid);  mu_p = mu(valid);  ci_p = ci(valid);
    in_xl = t_p >= xlim_line(1) & t_p <= xlim_line(2);
    t_p = t_p(in_xl);  mu_p = mu_p(in_xl);  ci_p = ci_p(in_xl);

    if numel(t_p) > 1
        fill(ax_fr, [t_p fliplr(t_p)], [mu_p+ci_p fliplr(mu_p-ci_p)], cols{ev}, ...
             'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax_fr, t_p, mu_p, '-', 'Color',cols{ev}, 'LineWidth',lw_main, ...
         'DisplayName', sprintf('%s (n=%d)', leg_names{ev}, n));
end

yyaxis(ax_fr,'right');
for ev = 1:3
    if isempty(mats_len_norm_fr{ev}), continue; end
    ln_mu = nanmean(mats_len_norm_fr{ev}, 1);
    ln_mu(isnan(ln_mu)) = 0;
    in_xl = t_common >= xlim_line(1) & t_common <= xlim_line(2);
    plot(ax_fr, t_common(in_xl), ln_mu(in_xl), ':', 'Color',cols{ev}, 'LineWidth',2, ...
         'HandleVisibility','off');
end
ax_fr.YAxis(2).Color = [0 0 0];
ylim(ax_fr, [0 1]);
ylabel(ax_fr,'Tongue length (norm)');
yyaxis(ax_fr,'left');

xline(ax_fr, 0, 'k--', 'LineWidth',1,'HandleVisibility','off');
xlim(ax_fr, xlim_line);
xlabel(ax_fr,'Time re. lick start (s)');
ylabel(ax_fr,'Mean spike rate (0-1 norm)');
title(ax_fr,'Mean spike rate','FontSize',11);
legend(ax_fr, 'Location','northwest','FontSize',9);
box(ax_fr,'off');
set(ax_fr,'TickDir','out','FontSize',11);


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
% 'mean' is the previous behavior: the same sum divided by the number of
% units. Because the weights are already normalized, that extra division just
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


function w = engagementModeWeights(aligned4, trials4, win1, win2, sdFloor)
% Weights defining the engagement mode, from the C4-reward trials.
%
% For each neuron: the difference between its mean firing rate in the window
% before the state transition and in the window after it, divided by the pooled
% across-trial standard deviation of those two windows. That is d', so a neuron
% counts for more only when it separates the two states reliably, rather than
% merely by a large number of spikes. Weights are normalized so sum|w| = 1.
%
% Each trial contributes one value per window (the bins inside the window are
% averaged first), so the standard deviations are taken across trials.

    N   = size(aligned4, 2);
    nTr = numel(trials4);

    preTr  = reshape(nanmean(aligned4(win1,:,trials4), 1), N, nTr);
    postTr = reshape(nanmean(aligned4(win2,:,trials4), 1), N, nTr);

% pooled SD: square each window's across-trial SD, average the two, square root
    sdPooled = sqrt((nanstd(preTr,0,2).^2 + nanstd(postTr,0,2).^2) / 2)';

    w = (nanmean(preTr,2)' - nanmean(postTr,2)') ./ max(sdPooled, sdFloor);

    w(~isfinite(w)) = 0;
    w = w / sum(abs(w));
end

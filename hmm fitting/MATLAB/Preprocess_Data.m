%% Load in the data
clear; clc;
datapth = '.\DataObjects'; % raw sesion data. Change to correct path as needed
outputBase = '.\Processed_Sessions'; % where per-session outputs are saved
addpath('.\zutils')

% Define loading params per session
anm = 'TD1d';
date = '2023-02-21';

% Load the datastructures
kinraw = load(fullfile(datapth, sprintf('%s_%s_kin.mat', anm, strrep(date, '-', '_'))));
objraw = load(fullfile(datapth, sprintf('%s_%s_obj.mat', anm, strrep(date, '-', '_'))));

% Select GC aligned datastructure
kin = kinraw.kinGC;
obj = objraw.objGC;

clear kinraw objraw;

%% Set parameters for processing
% the datafile was generated with -2.5 to 32s
tmin=-2.0;  % s
tmax=5.0;  % s
dt = 1/100;
pre_gc_points = -tmin / dt;
GaussianKernelSize = 10;  % samples

%% Extract the keypoint data
traj_features = [{'tongue_xdisp_view1'}, {'left_tongue_xdisp_view1'}, {'right_tongue_xdisp_view1'}, ...
    {'jaw_xdisp_view1'}, {'trident_xdisp_view1'}, {'nose_xdisp_view1'}, {'tongue_ydisp_view1'}, ...
    {'left_tongue_ydisp_view1'}, {'right_tongue_ydisp_view1'}, {'jaw_ydisp_view1'}, {'trident_ydisp_view1'}, ...
    {'nose_ydisp_view1'}, {'top_tongue_xdisp_view2'}, {'topleft_tongue_xdisp_view2'}, ...
    {'bottom_tongue_xdisp_view2'}, {'bottomleft_tongue_xdisp_view2'}, {'jaw_xdisp_view2'}, ...
    {'top_nostril_xdisp_view2'}, {'bottom_nostril_xdisp_view2'}, {'top_tongue_ydisp_view2'}, ...
    {'topleft_tongue_ydisp_view2'}, {'bottom_tongue_ydisp_view2'}, {'bottomleft_tongue_ydisp_view2'}, ...
    {'jaw_ydisp_view2'}, {'top_nostril_ydisp_view2'}, {'bottom_nostril_ydisp_view2'}];

[found, idx] = ismember(traj_features, kin.featLeg);

% pos = kin.dat100(:, :, idx);   % time x trials x numel(traj_features)

pos_full = kin.dat100(:, :, idx);   % time x trials x numel(traj_features), full -2.5:32s window

t = kin.time100(:);                 % actual time vector for dat100
keep = t >= tmin & t <= tmax;

pos  = pos_full(keep, :, :);
pos = permute(pos, [1 3 2]);   % time x features x trials
time = t(keep);

%% Get tongue lengths and contacts
% Get the tongue lengths for the desired trials
% Condition 8 is R1 Hits, Condition 9 is R4 Hits

condix = 1;  % all trials
sessix = 1;  % 1 session at a time
kinix = find(strcmp(kin(sessix).featLeg,'tongue_length'));  
all_length = kin.dat100(:, :, kinix);
all_length =  all_length(keep,:,:);


all_contacts = obj.bp.ev.lickL;
for i = 1:obj.bp.Ntrials
    contacts = obj.bp.ev.lickL{i, 1};  % lickL contains all contacts
    gc = obj.bp.ev.goCue(i);
    all_contacts{i} = contacts - gc;
end

R1_Trials = obj.trialid{8};  % Rmv'd a check here <= NTRIALS 9/17/26
R4_Trials = obj.trialid{9};

R1_Trial_Track = R1_Trials;
R4_Trial_Track = R4_Trials;

R1_Contacts = all_contacts(R1_Trials);
R4_Contacts = all_contacts(R4_Trials);

R1_Tongue = all_length((pre_gc_points-100+1):end, R1_Trials);
R4_Tongue = all_length((pre_gc_points-100+1):end, R4_Trials);

%% Filter trials based on lick port contacts
[trials2removeR1, FCs_R1_clean, ~, ~, ~, LRCs_R1_clean] = filter_trials_by_licking(R1_Contacts, 100, min_licks=3);
[trials2removeR4, FCs_R4_clean, ~, ~, ~, LRCs_R4_clean] = filter_trials_by_licking(R4_Contacts, 100, min_licks=5);

FCs_R1_clean(trials2removeR1) = [];
LRCs_R1_clean(trials2removeR1) = [];

FCs_R4_clean(trials2removeR4) = [];
LRCs_R4_clean(trials2removeR4) = [];

R1_Trial_Track(trials2removeR1) = [];
R4_Trial_Track(trials2removeR4) = [];

%% Keypoints
R1_Keypoints = pos(101:end, :, R1_Trials);  % trim first 100 points to get -1 to 5
R4_Keypoints = pos(101:end, :, R4_Trials);
R1_Keypoints(:,:,trials2removeR1) = [];
R4_Keypoints(:,:,trials2removeR4) = [];

R1_Keypoints_Uncut = reshape(permute(R1_Keypoints, [1, 3, 2]), [], size(R1_Keypoints, 2));
R4_Keypoints_Uncut = reshape(permute(R4_Keypoints, [1, 3, 2]), [], size(R4_Keypoints, 2));
R1K = zscore(R1_Keypoints_Uncut);
R4K = zscore(R4_Keypoints_Uncut);

n_time = size(R1_Keypoints, 1);
n_keypoints = size(R1_Keypoints, 2);
n_trials = size(R1_Keypoints, 3);
R1K_final = permute(reshape(R1K, n_time, n_trials, n_keypoints), [1, 3, 2]);

n_time = size(R4_Keypoints, 1);
n_keypoints = size(R4_Keypoints, 2);
n_trials = size(R4_Keypoints, 3);
R4K_final = permute(reshape(R4K, n_time, n_trials, n_keypoints), [1, 3, 2]);

R1_Keypoints_Uncut = R1K;
R4_Keypoints_Uncut = R4K;
R1_Keypoints_Cut = chop_and_stack_neural_data(R1K_final, LRCs_R1_clean, 100);
R4_Keypoints_Cut = chop_and_stack_neural_data(R4K_final, LRCs_R4_clean, 100);

%% Neural Data Processing
Ncells = size(obj.psth100, 2);

if iscell(obj.cluid) && numel(obj.cluid) > 1
    % Two probes
    probe1 = 1:numel(obj.cluid{1, 1});
    probe2 = numel(obj.cluid{1, 1}) + 1 : Ncells;
else
    % One probe
    probe1 = 1:numel(obj.cluid{1,1});
    probe2 = []; % No second probe
end

% Build trialdat slices based on detected probes
probe1_trialdat = obj.trialdat100(:, probe1, :);

% Cut to tmin tmax
probe1_trialdat = probe1_trialdat(keep,:,:);


% Gaussian filter the neural data
sz = size(probe1_trialdat);
x2 = reshape(probe1_trialdat, sz(1), []);
out2 = mySmooth(x2, GaussianKernelSize, 'reflect');
probe1_trialdat = reshape(out2, sz);


if ~isempty(probe2)
    probe2_trialdat = obj.trialdat100(:, probe2, :);
    
    probe2_trialdat = probe2_trialdat(keep,:,:);

    sz = size(probe2_trialdat);
    x2 = reshape(probe2_trialdat, sz(1), []);
    out2 = mySmooth(x2, GaussianKernelSize, 'reflect');
    probe2_trialdat = reshape(out2, sz);
end


% Display the number of probes detected
if ~isempty(probe2)
    nProbes = 2;
else
    nProbes = 1;
end
disp(['Number of probes: ', num2str(nProbes)]);

%% PCA for probe1
probe1_segment = probe1_trialdat(101:end, :, :);
[num_timepoints1, ~, num_trials1] = size(probe1_segment);
probe1_PCA = reshape(permute(probe1_segment, [1, 3, 2]), [], size(probe1_segment, 2));
[coeff1, score1] = pca(zscore(probe1_PCA));
score1 = score1(:, 1:10);
score1_reshaped = reshape(score1, num_timepoints1, num_trials1, 10);

Probe1_PCs_R1 = score1_reshaped(:, R1_Trials, :);
Probe1_PCs_R4 = score1_reshaped(:, R4_Trials, :);
Probe1_PCs_R1(:, trials2removeR1, :) = [];
Probe1_PCs_R4(:, trials2removeR4, :) = [];
Probe1_PCs_R1_Uncut = reshape(Probe1_PCs_R1, [], size(Probe1_PCs_R1, 3));
Probe1_PCs_R4_Uncut = reshape(Probe1_PCs_R4, [], size(Probe1_PCs_R4, 3));
Probe1_PCs_R4_Cut = chop_and_stack_neural_data(permute(Probe1_PCs_R4, [1, 3, 2]), LRCs_R4_clean, 100);
Probe1_PCs_R1_Cut = chop_and_stack_neural_data(permute(Probe1_PCs_R1, [1, 3, 2]), LRCs_R1_clean, 100);

if ~isempty(probe2)
    % PCA for probe2
    probe2_segment = probe2_trialdat(101:end, :, :);
    [num_timepoints2, ~, num_trials2] = size(probe2_segment);
    probe2_PCA = reshape(permute(probe2_segment, [1, 3, 2]), [], size(probe2_segment, 2));
    [coeff2, score2] = pca(zscore(probe2_PCA));
    score2 = score2(:, 1:10);
    score2_reshaped = reshape(score2, num_timepoints2, num_trials2, 10);

    Probe2_PCs_R1 = score2_reshaped(:, R1_Trials, :);
    Probe2_PCs_R4 = score2_reshaped(:, R4_Trials, :);
    Probe2_PCs_R1(:, trials2removeR1, :) = [];
    Probe2_PCs_R4(:, trials2removeR4, :) = [];
    Probe2_PCs_R1_Uncut = reshape(Probe2_PCs_R1, [], size(Probe2_PCs_R1, 3));
    Probe2_PCs_R4_Uncut = reshape(Probe2_PCs_R4, [], size(Probe2_PCs_R4, 3));
    Probe2_PCs_R4_Cut = chop_and_stack_neural_data(permute(Probe2_PCs_R4, [1, 3, 2]), LRCs_R4_clean, 100);
    Probe2_PCs_R1_Cut = chop_and_stack_neural_data(permute(Probe2_PCs_R1, [1, 3, 2]), LRCs_R1_clean, 100);
end

%% Get the tongue length, FCs, and LRCs times to save
SR = 100;
% Filter Tongue Length
R1_Tongue_Uncut = all_length(pre_gc_points-100+1:end, R1_Trials);  % -1s through 4s (500 points)
R4_Tongue_Uncut = all_length(pre_gc_points-100+1:end, R4_Trials);

R1_Tongue_Uncut(:,trials2removeR1) = [];
R4_Tongue_Uncut(:,trials2removeR4) = [];

FCs_Adj_R1 = ceil(FCs_R1_clean*SR + SR);
FCs_Adj_R4 = ceil(FCs_R4_clean*SR + SR);

LRCs_Adj_R1 = ceil(LRCs_R1_clean*SR + SR);
LRCs_Adj_R4 = ceil(LRCs_R4_clean*SR + SR);

%% Save Data
sessionName = anm;
sessionDate = date;

outputFolder = fullfile(outputBase, [sessionName '_' sessionDate]);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Save Trial Tracks
csvwrite(fullfile(outputFolder, "R1_Trial_Track.csv"), R1_Trial_Track);
csvwrite(fullfile(outputFolder, "R4_Trial_Track.csv"), R4_Trial_Track);

%% Save Key Point Features
csvwrite(fullfile(outputFolder, "Keypoint_Feats_R1_Uncut.csv"), R1_Keypoints_Uncut);
csvwrite(fullfile(outputFolder, "Keypoint_Feats_R4_Uncut.csv"), R4_Keypoints_Uncut);
csvwrite(fullfile(outputFolder, "Keypoint_Feats_R1_Cut.csv"), R1_Keypoints_Cut);
csvwrite(fullfile(outputFolder, "Keypoint_Feats_R4_Cut.csv"), R4_Keypoints_Cut);

%% Save Neural PCs for Probe 1
csvwrite(fullfile(outputFolder, "PCA_Probe1_R1_Uncut.csv"), Probe1_PCs_R1_Uncut);
csvwrite(fullfile(outputFolder, "PCA_Probe1_R4_Uncut.csv"), Probe1_PCs_R4_Uncut);
csvwrite(fullfile(outputFolder, "PCA_Probe1_R1_Cut.csv"), Probe1_PCs_R1_Cut);
csvwrite(fullfile(outputFolder, "PCA_Probe1_R4_Cut.csv"), Probe1_PCs_R4_Cut);

%% Save Neural PCs for Probe 2 (only if nProbes > 1)
if nProbes > 1
    csvwrite(fullfile(outputFolder, "PCA_Probe2_R1_Uncut.csv"), Probe2_PCs_R1_Uncut);
    csvwrite(fullfile(outputFolder, "PCA_Probe2_R4_Uncut.csv"), Probe2_PCs_R4_Uncut);
    csvwrite(fullfile(outputFolder, "PCA_Probe2_R1_Cut.csv"), Probe2_PCs_R1_Cut);
    csvwrite(fullfile(outputFolder, "PCA_Probe2_R4_Cut.csv"), Probe2_PCs_R4_Cut);
end

%% Save Tongue Length for visualizations
csvwrite(fullfile(outputFolder, "Tongue_R1.csv"), R1_Tongue_Uncut);
csvwrite(fullfile(outputFolder, "Tongue_R4.csv"), R4_Tongue_Uncut);

%% Save FCs and LRCs
csvwrite(fullfile(outputFolder, "FCs_R1.csv"), FCs_Adj_R1);
csvwrite(fullfile(outputFolder, "FCs_R4.csv"), FCs_Adj_R4);
csvwrite(fullfile(outputFolder, "LRCs_R1.csv"), LRCs_Adj_R1);
csvwrite(fullfile(outputFolder, "LRCs_R4.csv"), LRCs_Adj_R4);

%% Save metadata as a .txt file for record-keeping
metadataFile = fullfile(outputFolder, 'metadata.txt');
fid = fopen(metadataFile, 'w');
fprintf(fid, 'Processing Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Script Name: %s\n', mfilename('fullpath'));
fprintf(fid, 'Session ID: %s\n', sessionName);
fprintf(fid, 'Session Date: %s\n', sessionDate);
fclose(fid);

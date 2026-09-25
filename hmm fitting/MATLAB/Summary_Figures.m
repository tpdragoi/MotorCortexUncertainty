% Cross-session, cross-animal figure creation for disengagement analysis:
%   - Per session, each trial's inferred state is aligned to its first
%     contact (FC) time and averaged within that session
%   - Sessions are then averaged per animal, and animals are averaged
%     together for a grand mean across animals
%   - Line plot: faded per-animal averages overlaid with the bold grand
%     mean across animals
%   - Cloud plot: grand mean +/- std cloud (across animals) in place of
%     the per-animal lines

clear; clc; close all;

% Configuration
base_dir     = ".\Results";
alt_base_dir = ".\Processed_Sessions";
subfolder      = '';

fig_dir = fullfile(base_dir, 'Summary_Figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Load session metadata
metadata   = readtable('MATLAB\Session_Metadata.xlsx');
sessionList = metadata.Filename;
nSessions   = numel(sessionList);

% Alignment params
pre_points     = 10;
post_points    = 200;
window_size    = pre_points + 1 + post_points;
GCidx          = pre_points + 1;

% Plotting params
R1_color       = [0 0 1];
R4_color       = [0 0.76 1];
alpha_sess     = 0.5;
alpha_mean     = 1.0;
maxPostVisible = 150;
timeSec        = (-pre_points : post_points) / 100;
vis_idx        = timeSec <= maxPostVisible/100;

% Loop over sessions and create summary inference structures
avgR1_sessions = nan(nSessions, window_size);
avgR4_sessions = nan(nSessions, window_size);

for s = 1:nSessions
    sessName = sessionList{s};
    fprintf('Processing session %s ...\n', sessName);

    save_dir = fullfile(base_dir, sessName, subfolder);
    alt_dir  = fullfile(alt_base_dir, sessName);

    if ~isfolder(save_dir) || ~isfolder(alt_dir)
        warning('Missing dirs for %s, skipping.', sessName);
        continue;
    end

    R1s = readmatrix(fullfile(save_dir, 'R1_States_Reg.csv'));
    R4s = readmatrix(fullfile(save_dir, 'R14_States_Reg.csv'));
    C1  = readmatrix(fullfile(alt_dir,  'FCs_R1.csv')) - 100;
    C4  = readmatrix(fullfile(alt_dir,  'FCs_R4.csv')) - 100;

    % Align R1 trials 
    nT1 = size(R1s, 1);
    aligned1 = nan(nT1, window_size);
    for t = 1:nT1
        ap  = GCidx + round(C1(t));
        st  = ap - pre_points;  en = ap + post_points;
        vs  = max(st, 1);       ve = min(en, size(R1s, 2));
        is_ = vs - st + 1;      ie_ = is_ + (ve - vs);
        aligned1(t, is_:ie_) = exp(R1s(t, vs:ve));
    end
    avgR1_sessions(s, :) = nanmean(aligned1, 1);

    % Align R4 trials 
    nT4 = size(R4s, 1);
    aligned4 = nan(nT4, window_size);
    for t = 1:nT4
        ap  = GCidx + round(C4(t));
        st  = ap - pre_points;  en = ap + post_points;
        vs  = max(st, 1);       ve = min(en, size(R4s, 2));
        is_ = vs - st + 1;      ie_ = is_ + (ve - vs);
        aligned4(t, is_:ie_) = exp(R4s(t, vs:ve));
    end
    avgR4_sessions(s, :) = nanmean(aligned4, 1);

    fprintf('  -> done\n');
end

%% Combine across sessions: average + std, then plot
% Drop skipped sessions
validRows = ~all(isnan(avgR1_sessions), 2);
matR1 = avgR1_sessions(validRows, :);
matR4 = avgR4_sessions(validRows, :);
validNames = sessionList(validRows);

% Calculate average inference per animal across sessions
animalNames = cellfun(@(s) extractBefore(s, '_'), validNames, 'UniformOutput', false);
uniqueAnimals = unique(animalNames, 'stable');

nAnimals = numel(uniqueAnimals);
animalR1 = nan(nAnimals, window_size);
animalR4 = nan(nAnimals, window_size);

for a = 1:nAnimals
    idx = strcmp(animalNames, uniqueAnimals{a});
    animalR1(a, :) = mean(matR1(idx, :), 1);
    animalR4(a, :) = mean(matR4(idx, :), 1);
end

% Overall mean and std across animals
meanR1 = mean(animalR1, 1);   stdR1 = std(animalR1, 0, 1);
meanR4 = mean(animalR4, 1);   stdR4 = std(animalR4, 0, 1);

% Line plot
f1 = figure; hold on;
for a = 1:nAnimals
    h1 = plot(timeSec(vis_idx), animalR1(a, vis_idx), 'LineWidth', 1.5);
    h1.Color = [R1_color, alpha_sess];
    h4 = plot(timeSec(vis_idx), animalR4(a, vis_idx), 'LineWidth', 1.5);
    h4.Color = [R4_color, alpha_sess];
end
hGM1 = plot(timeSec(vis_idx), meanR1(vis_idx), 'LineWidth', 4, 'Color', [R1_color, alpha_mean]);
hGM4 = plot(timeSec(vis_idx), meanR4(vis_idx), 'LineWidth', 4, 'Color', [R4_color, alpha_mean]);
xline(0, '--k', 'LineWidth', 1);
xlabel('Time from first port contact (s)');
ylabel('Engaged State Probability');
title('Average Inference across animals', 'Interpreter', 'none');
xlim([timeSec(find(vis_idx,1,'first')), timeSec(find(vis_idx,1,'last'))]);
ylim([0, 1]);
xticks([-0.1, 0, 0.5, 1.0, 1.5]);
yticks(0:0.2:0.8);
grid on;
legend([hGM1, hGM4], {'R1 mean across animals', 'R4 mean across animals'}, ...
       'Location', 'Best', 'Interpreter', 'none');
saveas(f1,  fullfile(fig_dir, 'AverageInference_Line.png'));
savefig(f1, fullfile(fig_dir, 'AverageInference_Line.fig'));

% Std Cloud Plot
f2 = figure; hold on;
alpha_cloud = 0.3;
x   = timeSec(vis_idx);
y1u = meanR1(vis_idx) + stdR1(vis_idx);
y1l = meanR1(vis_idx) - stdR1(vis_idx);
hC1 = fill([x, fliplr(x)], [y1u, fliplr(y1l)], R1_color, 'FaceAlpha', alpha_cloud, 'EdgeColor', 'none');
y4u = meanR4(vis_idx) + stdR4(vis_idx);
y4l = meanR4(vis_idx) - stdR4(vis_idx);
hC4 = fill([x, fliplr(x)], [y4u, fliplr(y4l)], R4_color, 'FaceAlpha', alpha_cloud, 'EdgeColor', 'none');
hGM1 = plot(x, meanR1(vis_idx), 'LineWidth', 4, 'Color', [R1_color, alpha_mean]);
hGM4 = plot(x, meanR4(vis_idx), 'LineWidth', 4, 'Color', [R4_color, alpha_mean]);
xline(0, '--k', 'LineWidth', 1);
xlabel('Time from first port contact (s)');
ylabel('Engaged State Probability');
title('Average Inference across animals', 'Interpreter', 'none');
xlim([x(1), x(end)]);
ylim([0, 1]);
xticks([-0.1, 0, 0.5, 1.0, 1.5]);
yticks(0:0.2:0.8);
grid on;
legend([hGM1, hGM4], {'R1 mean ±1 SD', 'R4 mean ±1 SD'}, ...
       'Location', 'Best', 'Interpreter', 'none');
saveas(f2,  fullfile(fig_dir, 'AverageInference_Cloud.png'));
savefig(f2, fullfile(fig_dir, 'AverageInference_Cloud.fig'));savefig(f2, fullfile(fig_dir, 'AverageInference_Cloud.fig'));

%% Summary: animals and session counts
fprintf('\n--- Sessions used per animal ---\n');
for a = 1:nAnimals
    idx = strcmp(animalNames, uniqueAnimals{a});
    fprintf('  %s: %d sessions\n', uniqueAnimals{a}, sum(idx));
end
fprintf('  Total sessions: %d\n', numel(validNames));
fprintf('  Total animals:  %d\n', nAnimals);
% Per-session figure creation for disengagement analysis:
%   - Trial-averaged inference (mean +/- std cloud), GC-aligned
%   - Trial-averaged inference (mean +/- std cloud), aligned to first
%     contact (FC)
%   - Single-trial state heatmap (GC-aligned) with tongue traces and
%     estimated disengagement timepoints overlaid

clear; clc; close all

%% Setup
base_dir     = ".\\Results";
alt_base_dir = ".\\Processed_Sessions";
subfolder    = '';

% Get list of all subfolders in base_dir
session_dirs = dir(base_dir);
session_dirs = session_dirs([session_dirs.isdir]);
session_dirs = session_dirs(~ismember({session_dirs.name}, {'.', '..'}));

% Plotting params
nColors      = 256;
pure_blue    = [0, 0, 1];
pure_red     = [1, 0, 0];
white        = [1, 1, 1];
blend_factor = 0.5;
blue_side    = (1 - blend_factor) * pure_blue + blend_factor * white;
red_side     = (1 - blend_factor) * pure_red  + blend_factor * white;
custom_cmap  = [linspace(blue_side(1), red_side(1), nColors)', ...
                linspace(blue_side(2), red_side(2), nColors)', ...
                linspace(blue_side(3), red_side(3), nColors)'];

pregc = 100;
start_time = 90;  % matches Julia's START_TIME

for ij = 1:length(session_dirs)
    session_name = session_dirs(ij).name;
    save_dir = fullfile(base_dir, session_name, subfolder);
    alt_session_dir = fullfile(alt_base_dir, session_name);

    if isfolder(alt_session_dir)
        fprintf('Found matching folder for %s in alt_base_dir.\n', session_name);
    else
        warning('No matching folder for %s in alt_base_dir.', session_name);
        continue;
    end

    %% Load state/tongue traces and first-contact times
    R1_States = readmatrix(fullfile(save_dir, 'R1_States_Reg.csv'));
    R4_States = readmatrix(fullfile(save_dir, 'R14_States_Reg.csv'));
    R1_Tongue = readmatrix(fullfile(save_dir, 'R1_Tongue_Reg.csv'));
    R4_Tongue = readmatrix(fullfile(save_dir, 'R14_Tongue_Reg.csv'));

    FCs_R1 = readmatrix(fullfile(alt_session_dir, 'FCs_R1.csv')) - pregc;
    FCs_R4 = readmatrix(fullfile(alt_session_dir, 'FCs_R4.csv')) - pregc;

    %% Combine R4/R1 for the heatmap, normalize tongue traces per trial
    All_States = exp([R4_States; R1_States]);

    R4_Tongue_norm = (R4_Tongue - nanmin(R4_Tongue, [], 2)) ./ (nanmax(R4_Tongue, [], 2) - nanmin(R4_Tongue, [], 2));
    R4_Tongue_norm(R4_Tongue_norm == 0) = NaN;

    R1_Tongue_norm = (R1_Tongue - nanmin(R1_Tongue, [], 2)) ./ (nanmax(R1_Tongue, [], 2) - nanmin(R1_Tongue, [], 2));
    R1_Tongue_norm(R1_Tongue_norm == 0) = NaN;

    R1_Tongue = R1_Tongue_norm';
    R4_Tongue = R4_Tongue_norm';
    R1_States = R1_States';
    R4_States = R4_States';

    %% Estimate disengagement time per trial (needed for the heatmap markers below)
    % Time_Reg.csv is seconds relative to true GC, index-aligned column-for-column
    time_reg = readmatrix(fullfile(save_dir, 'Time_Reg.csv'));
    GCtime = find(time_reg >= 0, 1, 'first');
    tick_times = ceil(time_reg(1)*4)/4 : 0.25 : floor(time_reg(end)*4)/4;
    tick_positions = interp1(time_reg, 1:numel(time_reg), tick_times, 'nearest');

    engage_thresh = 0.5;        % Probability threshold for engagement
    min_block_length = 2;       % Minimum consecutive timepoints for detecting largest engagement block
    min_disengaged_gap = 20;    % Number of consecutive timepoints below threshold to define a final disengagement

    R1_idx = estimate_disengage_times(exp(R1_States'), engage_thresh, min_block_length, min_disengaged_gap);
    R4_idx = estimate_disengage_times(exp(R4_States'), engage_thresh, min_block_length, min_disengaged_gap);

    % ms relative to the GC bin (column GCtime)
    R1_disengage = nan(size(R1_idx));
    R1_disengage(~isnan(R1_idx)) = (time_reg(R1_idx(~isnan(R1_idx))) - time_reg(GCtime)) * 1000;

    R4_disengage = nan(size(R4_idx));
    R4_disengage(~isnan(R4_idx)) = (time_reg(R4_idx(~isnan(R4_idx))) - time_reg(GCtime)) * 1000;

    csvwrite(fullfile(save_dir, "R1_dt.csv"), R1_disengage);
    csvwrite(fullfile(save_dir, "R4_dt.csv"), R4_disengage);

    %% Trial-averaged inference (GC-aligned, mean +/- std cloud)
    R1_Inf_Mean = mean(exp(R1_States'), 1);
    R1_Inf_Std  = std(exp(R1_States'), 0, 1);

    R4_Inf_Mean = mean(exp(R4_States'), 1);
    R4_Inf_Std  = std(exp(R4_States'), 0, 1);

    x = 1:length(R1_Inf_Mean);

    figure
    hold on

    fill([x, fliplr(x)], ...
         [R1_Inf_Mean + R1_Inf_Std, fliplr(R1_Inf_Mean - R1_Inf_Std)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    fill([x, fliplr(x)], ...
         [R4_Inf_Mean + R4_Inf_Std, fliplr(R4_Inf_Mean - R4_Inf_Std)], ...
         [0.7 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    plot(x, R1_Inf_Mean, 'b', 'LineWidth', 2)
    plot(x, R4_Inf_Mean, 'Color', [0 0.76 1], 'LineWidth', 2)

    ylabel("State 1 Probability")
    xlabel("Time (s)")
    xticks(tick_positions)
    xticklabels(string(tick_times))
    ylim([0, 1])
    title("Trial Averaged Inference (GC-aligned)")

    xline(GCtime, "--k", "LineWidth", 1)
    text(GCtime, min([R1_Inf_Mean - R1_Inf_Std, R4_Inf_Mean - R4_Inf_Std], [], 'all') - 0.02, 'GC', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10);

    legend(["R1 Std", "R4 Std", "R1 Mean", "R4 Mean"])
    saveas(gcf, fullfile(save_dir, 'Ave_Inference.png'))
    saveas(gcf, fullfile(save_dir, 'Ave_Inference.fig'))

    %% Trial-averaged inference aligned to first contact (FC), mean +/- std cloud
    window_size = 211;
    pre_points  = 10;
    post_points = 200;

    n_trials = size(R1_States, 2);
    R1_FC_aligned = nan(window_size, n_trials);
    for i = 1:n_trials
        % FCs_R1(i) = raw_row - pregc; R1_States is indexed by local column
        % (raw_row - start_time + 1), so convert: + (pregc - start_time + 1).
        align_point = round(FCs_R1(i)) + (pregc - start_time + 1);
        start_idx = align_point - pre_points;
        end_idx   = align_point + post_points;

        valid_start = max(start_idx, 1);
        valid_end   = min(end_idx, size(R1_States, 1));

        insert_start = valid_start - start_idx + 1;
        insert_end   = insert_start + (valid_end - valid_start);

        R1_FC_aligned(insert_start:insert_end, i) = R1_States(valid_start:valid_end, i);
    end

    n_trials = size(R4_States, 2);
    R4_FC_aligned = nan(window_size, n_trials);
    for i = 1:n_trials
        align_point = round(FCs_R4(i)) + (pregc - start_time + 1);
        start_idx = align_point - pre_points;
        end_idx   = align_point + post_points;

        valid_start = max(start_idx, 1);
        valid_end   = min(end_idx, size(R4_States, 1));

        insert_start = valid_start - start_idx + 1;
        insert_end   = insert_start + (valid_end - valid_start);

        R4_FC_aligned(insert_start:insert_end, i) = R4_States(valid_start:valid_end, i);
    end

    R1_Inf_Mean = nanmean(exp(R1_FC_aligned'), 1);
    R1_Inf_Std  = nanstd(exp(R1_FC_aligned'), 0, 1);

    R4_Inf_Mean = nanmean(exp(R4_FC_aligned'), 1);
    R4_Inf_Std  = nanstd(exp(R4_FC_aligned'), 0, 1);

    % Replace NaNs (timepoints with no valid trials) with 0.0
    R1_Inf_Mean(isnan(R1_Inf_Mean)) = 0.0;
    R1_Inf_Std(isnan(R1_Inf_Std))   = 0.0;
    R4_Inf_Mean(isnan(R4_Inf_Mean)) = 0.0;
    R4_Inf_Std(isnan(R4_Inf_Std))   = 0.0;

    x = 1:length(R1_Inf_Mean);

    figure
    hold on

    % Clamp shaded areas between 0 and 1
    R1_Shade_Upper = min(R1_Inf_Mean + R1_Inf_Std, 1);
    R1_Shade_Lower = max(R1_Inf_Mean - R1_Inf_Std, 0);
    R4_Shade_Upper = min(R4_Inf_Mean + R4_Inf_Std, 1);
    R4_Shade_Lower = max(R4_Inf_Mean - R4_Inf_Std, 0);

    r1_s = fill([x, fliplr(x)], ...
         [R1_Shade_Upper, fliplr(R1_Shade_Lower)], ...
         [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    r4_s = fill([x, fliplr(x)], ...
         [R4_Shade_Upper, fliplr(R4_Shade_Lower)], ...
         [0.7 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    r1_p = plot(x, R1_Inf_Mean, 'b', 'LineWidth', 2);
    r4_p = plot(x, R4_Inf_Mean, 'Color', [0 0.76 1], 'LineWidth', 2);

    ylabel("State 1 Probability")
    xlabel("Time (s)")
    xticks([0 60 110 160 210]);
    xticklabels({'-0.1', '0.5', '1.0', '1.5', '2.0'});
    xlim([0, 200])
    ylim([0, 1])
    title("Trial Averaged Inference (FC-aligned)")

    xline(11, "--k", "LineWidth", 1)
    text(11, min([R1_Inf_Mean - R1_Inf_Std, R4_Inf_Mean - R4_Inf_Std], [], 'all') - 0.02, 'FC', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10);

    legend([r1_p, r4_p, r1_s, r4_s], {'R1 Mean', 'R4 Mean', 'R1 Std', 'R4 Std'})
    saveas(gcf, fullfile(save_dir, 'Ave_Inference_FCAligned.png'))
    saveas(gcf, fullfile(save_dir, 'Ave_Inference_FCAligned.fig'))

    %% Single-trial heatmap (GC-aligned) with tongue traces + disengagement markers
    figure;
    hold on;

    lw = 0.75;    % Line width
    px = 75; py = 75;
    width = 700; height = 800;
    set(gcf, 'Position', [px, py, width, height]);

    R4_color = 'k';
    R1_color = 'k';

    [trs, ~] = size(All_States);

    % Overlay heatmap
    imagesc(1:size(All_States, 2), 1:size(All_States, 1), All_States);

    % Tongue traces
    for j = 1:size(R4_Tongue, 2)
        plot(1:length(R4_Tongue(:, j)), j-1 + R4_Tongue(:, j), R4_color, 'LineWidth', lw);
    end
    for j = 1:size(R1_Tongue, 2)
        plot(1:length(R1_Tongue(:, j)), j-1 + size(R4_Tongue, 2) + R1_Tongue(:, j), R1_color, 'LineWidth', lw);
    end

    % Disengagement time markers
    for i = 1:length(R4_idx)
        if ~isnan(R4_idx(i))
            plot(R4_idx(i), i, 'wo', 'MarkerSize', 4, 'LineWidth', 1.2);
        end
    end
    for i = 1:length(R1_idx)
        if ~isnan(R1_idx(i))
            y = i + size(R4_Tongue, 2);    % Offset for R1 trials
            plot(R1_idx(i), y, 'wo', 'MarkerSize', 4, 'LineWidth', 1.2);
        end
    end

    colormap(custom_cmap);
    colorbar;

    % Horizontal separator between R4 and R1 trials
    yline(size(R4_Tongue, 2), 'k-', 'LineWidth', 3);

    text(-30, size(R4_Tongue, 2) + 2, 'R1', 'FontSize', 12, 'Color', "k");
    text(-30, size(R4_Tongue, 2) - 2, 'R4', 'FontSize', 12, 'Color', "k");

    set(gca, 'YTick', 0:10:trs);
    xlabel('Time (s)');
    ylabel('Trial Number');

    xticks(tick_positions);
    xticklabels(string(tick_times));
    nTrials = size(All_States, 1);
    ylim([0 nTrials]);

    box off;
    axis tight;
    set(gca, 'TickLength', [0 0]);
    hold off;

    title("Single Trial State Estimates: R1 and R4 with Disengagement Times");
    xline(GCtime, '--k', 'LineWidth', 1);
    text(GCtime, -5, 'GC', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12, ...
        'Color', 'k');

    saveas(gcf, fullfile(save_dir, 'Inference_Heatmap_Disengagement.png'));
    saveas(gcf, fullfile(save_dir, 'Inference_Heatmap_Disengagement.fig'));

    close all;
end

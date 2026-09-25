function disengage_times = estimate_disengage_times(state_probs, engaged_thresh, min_block_length, min_disengaged_gap)
% estimate_disengage_times
% -------------------------
% Estimate disengagement time for each trial as the end of the primary engagement period.
%
% Arguments:
% - state_probs: trials x timepoints matrix of probabilities (e.g., exp(R1_States))
% - engaged_thresh: probability threshold above which a timepoint is considered "engaged"
% - min_block_length: minimum consecutive timepoints above threshold to count as "main" engagement
% - min_disengaged_gap: number of consecutive timepoints below threshold to define a final disengagement
%
% Returns:
% - disengage_times: [num_trials x 1] vector of time indices when the main engagement ends

[num_trials, num_timepoints] = size(state_probs);
disengage_times = nan(num_trials, 1);

for trial = 1:num_trials
    p = state_probs(trial, :);
    
    % Binary mask of engagement
    engaged = p > engaged_thresh;
    
    % Find start and end indices of engagement blocks
    d = diff([0 engaged 0]);
    start_inds = find(d == 1);
    end_inds = find(d == -1) - 1;
    run_lengths = end_inds - start_inds + 1;
    
    % Keep only main engagement blocks
    is_main = run_lengths >= min_block_length;
    main_blocks = [start_inds(is_main)', end_inds(is_main)'];
    
    if isempty(main_blocks)
        continue
    end
    
    % Take first main engagement block
    main_block = main_blocks(1, :);
    last_idx = main_block(2);
    
    % Look for sustained disengagement after main block
    for i = last_idx+1:num_timepoints - min_disengaged_gap
        if all(p(i:i+min_disengaged_gap-1) < engaged_thresh)
            last_idx = i;
            break
        end
    end
    
    disengage_times(trial) = last_idx;
end
end

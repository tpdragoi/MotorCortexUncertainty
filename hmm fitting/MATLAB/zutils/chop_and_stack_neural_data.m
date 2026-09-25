function final_neural_data = chop_and_stack_neural_data(neural_data, trial_end_times, SR)
    % neural_data: time x neurons x trials
    % trial_end_times: vector of chop points per trial
    [~, num_neurons, num_trials] = size(neural_data);
    total_rows = sum(ceil(trial_end_times*SR));
    
    final_neural_data = NaN(total_rows, num_neurons);
    current_row = 1;
    
    for trial_idx = 1:num_trials
        % Plus SR at the end because the end times are relative to the GC
        % and we have 1 second of data before the GC we want to include
        trial_length = ceil(SR*trial_end_times(trial_idx)) + SR;
        for neuron_idx = 1:num_neurons
            
            neural_trial = neural_data(1:trial_length, neuron_idx, trial_idx);
            
            end_row = current_row + trial_length - 1;
            final_neural_data(current_row:end_row, neuron_idx) = neural_trial;
        end
        
        current_row = current_row + trial_length;
    end
end

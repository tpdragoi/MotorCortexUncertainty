function [trials2remove, first_lick_times, second_lick_times, fourth_lick_times, sixth_lick_times, trial_end_times] = filter_trials_by_licking(contacts, SR, varargin)
%FILTER_TRIALS_BY_LICKING Filters trials based on lick behavior criteria.
%
%   [trials2remove, first_lick_times, trial_end_times] = filter_trials_by_licking(contacts)
%   [trials2remove, first_lick_times, trial_end_times] = filter_trials_by_licking(contacts, 'min_licks', 8, 'time_max', 480, 'gap_multiplier', 2)
%
% Inputs:
%   contacts - Cell array, each cell contains lick times for a trial
%
% Optional name-value pairs:
%   'min_licks'      - Minimum number of licks required (default: 8)
%   'time_max'       - Max time threshold in ms (default: 480)
%   'gap_multiplier' - Multiplier for inter-lick interval threshold (default: 2)
%
% Outputs:
%   trials2remove   - Vector of trial indices to remove
%   first_lick_times - Vector of first lick times (NaN if trial removed)
%   trial_end_times  - Vector of trial end times (NaN if trial removed)

    % Parse inputs
    p = inputParser;
    addParameter(p, 'min_licks', 5, @isnumeric);
    addParameter(p, 'time_max', 480, @isnumeric);
    addParameter(p, 'gap_multiplier', 2, @isnumeric);
    parse(p, varargin{:});

    min_licks = p.Results.min_licks;
    time_max = p.Results.time_max;
    gap_multiplier = p.Results.gap_multiplier;

    % Initialize outputs
    num_trials = length(contacts);
    trials2remove = [];
    first_lick_times = NaN(1, num_trials);
    second_lick_times = NaN(1, num_trials);
    fourth_lick_times = NaN(1, num_trials);
    sixth_lick_times = NaN(1, num_trials);
    trial_end_times = NaN(1, num_trials);

    for tr = 1:num_trials
        trial_contacts = contacts{tr};
        % Filter to time window
        trial_contacts = trial_contacts(trial_contacts > 0 & trial_contacts < time_max / SR);

        % Skip if not enough licks
        if length(trial_contacts) < min_licks
            trials2remove = [trials2remove, tr];
            continue
        end

        % Record first lick
        first_lick_times(tr) = trial_contacts(1);

        % Record the second lick
        second_lick_times(tr) = trial_contacts(2);

        if length(trial_contacts) >= 4
            fourth_lick_times(tr) = trial_contacts(4);
        else
            fourth_lick_times(tr) = NaN;
        end

        if length(trial_contacts) >= 6
            sixth_lick_times(tr) = trial_contacts(6);
        else
            sixth_lick_times(tr) = NaN;
        end

        % Check gaps between licks
        ILIs = diff(trial_contacts);
        median_ILI = median(ILIs);
        threshold = gap_multiplier * median_ILI;

        offense = find(ILIs > threshold, 1);

        if ~isempty(offense)
            if offense < min_licks
                trials2remove = [trials2remove, tr];
                continue
            else
                % Cut trial at last valid lick
                trial_contacts = trial_contacts(1:offense);
            end
        end

        % Record trial end time
        trial_end_times(tr) = trial_contacts(end);
    end
end

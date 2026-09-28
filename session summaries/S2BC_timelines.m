%% S2BC_timelines.m
%  Fig. S2B -- days of training before the first recording session, one point
%              per animal. Learning animals were recorded from first exposure,
%              so their training time is zero by definition, not missing.
%  Fig. S2C -- sessions each animal performed in each epoch of the optogenetic
%              experiments. The same four VGAT-ChR2-EYFP animals run through
%              every epoch in order.
%
%  READS  session summaries/sessionCounts.m

clear; clc; close all

rng(0);   % jitter is cosmetic; fixed so the panels redraw identically
jitterWidth = 0.4;

S = sessionCounts();

cols = [1.0 0.0 0.0;    % Single Reward
        0.0 0.7 1.0;    % Delayed Reward
        0.6 0.0 0.8;    % Double Reward
        1.0 0.4 0.0;    % Learning
        0.0 0.6 0.0];   % VTA

%% Fig. S2B

train = { [19 14 15 18], [16 12 17 20 33 23 22], [31 16 26 17 20 19], ...
          [0 0 0 0], [47 33 23 68] };

for i = 1:numel(train)
    assert(numel(train{i}) == S.nAnimals(i), ...
        'S2B: %s has %d animals, expected %d.', S.tasks{i}, numel(train{i}), S.nAnimals(i));
end

stripPlot(train, S.tasks, cols, jitterWidth, 50, 'Training time (days)');

%% Fig. S2C

epochs = {'Pretraining','Simple Reward L4','Simple Reward GC','Retrain','Double Reward L4'};
opto   = { [14 15 6 6], [5 7 3 3], [1 1 2 1], [6 4 4 4], [5 5 2 6] };

optoCols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.93 0.69 0.13];

for i = 1:numel(opto)
    assert(numel(opto{i}) == 4, 'S2C: %s has %d animals, expected 4.', epochs{i}, numel(opto{i}));
end

stripPlot(opto, epochs, optoCols, jitterWidth, 80, 'Number of sessions');

%% CHECK THE EPOCH TOTALS AGAINST THE PHOTOINHIBITION PANELS

fprintf('\nFig. S2B  training time before the first recording\n');
for i = 1:numel(train)
    fprintf('  %-16s n = %d   median %.0f\n', S.tasks{i}, numel(train{i}), median(train{i}));
end

fprintf('\nFig. S2C  sessions per epoch\n');
want = [NaN 18 4 NaN 18];   % Fig. 1J/1K, Fig. 1L/1M, Fig. 3J/3K/3L
for i = 1:numel(opto)
    got = sum(opto{i});
    flag = '';
    if ~isnan(want(i)) && got ~= want(i), flag = sprintf('   <-- panels use %d', want(i)); end
    fprintf('  %-18s %2d sessions, %d animals%s\n', epochs{i}, got, numel(opto{i}), flag);
end
fprintf('\n');


%% LOCAL FUNCTIONS

function stripPlot(groups, labels, cols, jitterWidth, markerSz, yLab)
    figure('Color','w'); hold on
    for i = 1:numel(groups)
        y = groups{i}(:)';
        scatter(i + (rand(1,numel(y)) - 0.5)*jitterWidth, y, markerSz, cols(i,:), ...
                'filled', 'MarkerFaceAlpha', 0.6);
    end
    xlim([0.5 numel(groups)+0.5]);
    set(gca, 'XTick', 1:numel(groups), 'XTickLabel', labels, 'FontSize', 12, ...
             'Box','off', 'TickDir','out');
    ylabel(yLab);
    hold off
    set(gcf, 'Position', [50 50 600 300]);
end

%% S2A_percentCompleted.m
%  Fig. S2A -- proportion of trials completed per session, one point per session.
%  Value plotted is sum(obj.bp.hit) / obj.bp.Ntrials, the same quantity
%  FS10_spontaneous.m stores in percCompleted. Values are pasted from those runs
%  so this file needs no data on disk; regenerate from the loaders if a session
%  list changes. Counts are checked against sessionCounts.m.
%
%  READS  session summaries/sessionCounts.m

clear; clc; close all

pc.simple   = [0.9963 0.9964 1.0000 0.9966 1.0000 0.9752 0.9905 0.9944 0.9922 ...
               0.9800 0.9769 0.9375 0.8525 0.9517 0.9050 0.9015 0.9696 0.8814 ...
               0.9771 0.8210 0.9727 0.9521];

pc.delayed  = [0.7289 0.8153 0.8421 0.6804 0.8049 0.8631 0.7835 0.8596 0.9893 ...
               1.0000 1.0000 0.9922 0.9717 0.9746 0.9719 0.9885 1.0000 0.9500 ...
               0.9874 0.9788 0.9841 0.9667 0.9800 0.9961 1.0000 0.9636 0.9653 ...
               0.9699 0.9840 0.9679 0.9173 0.9145 0.9598 0.9867 0.9808];

pc.double   = [0.9802 0.8113 0.8385 0.8138 0.7073 0.8576 0.8950 0.9766 0.9558 ...
               0.9771 0.9531 0.9305 0.9542 0.9145 0.9730 0.9922 0.9915 0.9917 ...
               0.9018 0.9689 0.9113 0.9451 0.8713 0.8955 0.9602 1.0000];

pc.learning = [0.6923 1.0000 1.0000 1.0000 1.0000 0.8034 0.9875 1.0000 1.0000 ...
               0.9955 0.8976 1.0000 1.0000 1.0000 1.0000 0.6797 0.9905 0.9909 ...
               0.9917 1.0000];

pc.vta      = [0.3756 0.5046 0.7273 0.6403 0.6183 0.7027 0.7154 0.5185 0.8716 ...
               0.7152 0.4843 0.5541 0.7450 0.7205 0.8500 0.7786 0.7727 0.5280];

%% SETTINGS

S = sessionCounts();
G = {'simple','delayed','double','learning','vta'};

cols = [1.0 0.0 0.0;    % Single Reward
        0.0 0.7 1.0;    % Delayed Reward
        0.6 0.0 0.8;    % Double Reward
        1.0 0.4 0.0;    % Learning
        0.0 0.6 0.0];   % VTA

for i = 1:numel(G)
    assert(numel(pc.(G{i})) == S.nSessions(i), ...
        '%s has %d values, expected %d sessions.', S.tasks{i}, numel(pc.(G{i})), S.nSessions(i));
    assert(all(pc.(G{i}) >= 0 & pc.(G{i}) <= 1), '%s: value outside [0 1].', S.tasks{i});
end

%% PLOT

rng(0);   % jitter is cosmetic; fixed so the panel redraws identically

figure('Color','w'); hold on
for i = 1:numel(G)
    y = pc.(G{i})(:);
    scatter(i + (rand(size(y)) - 0.5)*0.14, y, 60, 'MarkerFaceColor', cols(i,:), ...
            'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.5);
end
xlim([0.5 numel(G)+0.5]);  ylim([0 1]);
set(gca, 'XTick', 1:numel(G), 'XTickLabel', S.tasks, 'TickDir','out', 'Box','off', 'FontSize',12);
ylabel('Proportion of trials completed');
hold off
set(gcf, 'Position', [50 50 600 300]);

fprintf('\nFig. S2A\n');
for i = 1:numel(G)
    y = pc.(G{i});
    fprintf('  %-16s n = %2d   median %.3f   range %.3f-%.3f\n', ...
            S.tasks{i}, numel(y), median(y), min(y), max(y));
end
fprintf('\n');

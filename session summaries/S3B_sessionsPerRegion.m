%% S3B_sessionsPerRegion.m
%  Fig. S3B -- recorded sessions per region and per animal, across all tasks.
%  One point per animal, regions separated within each task. Counts come from
%  sessionCounts.m so this panel cannot drift from the sessions the analyses
%  ran on. A region an animal was not recorded in contributes no point.
%
%  READS  session summaries/sessionCounts.m

clear; clc; close all

S = sessionCounts();

tasks   = S.tasks;
regions = S.regions;

%% NUMBERS FOR THE LEGEND

fprintf('\nFig. S3B\n');
for r = 1:numel(regions)
    for k = 1:numel(tasks)
        v = S.(S.fields{k}).(regions{r});
        if isempty(v), continue; end
        fprintf('  %-5s %-16s %2d sessions, %d animals   %s\n', ...
                regions{r}, tasks{k}, sum(v), numel(v), mat2str(v));
    end
end
fprintf('\n');

%% PLOT

colors.tjM1 = [0.0000 0.4470 0.7410];
colors.ALM  = [0.4660 0.6740 0.1880];
colors.tjS1 = [0.9290 0.6940 0.1250];

spacing     = 15;        % task to task
offset.tjM1 = -5.6;      % region within a task
offset.ALM  = -0.8;
offset.tjS1 = +4.8;
offsetStep  = 1;         % sideways nudge for animals with the same count
ms          = 150;
faceAlpha   = 0.6;

yL = [-0.5 7];
xCenters = 1 + (0:numel(tasks)-1) * spacing;

figure('Color','w','Position',[200 200 900 450]); hold on

for k = 1:numel(tasks)-1                              % task separators
    xs = mean(xCenters(k:k+1));
    line([xs xs], yL, 'Color',[0.8 0.8 0.8], 'LineStyle','--', 'LineWidth',1);
end
for k = 1:numel(tasks)                                % region separators
    sx = xCenters(k) + [offset.tjM1, offset.ALM, offset.tjS1];
    line(repmat(mean(sx(1:2)),1,2), yL, 'Color',[0.6 0.6 0.6], 'LineStyle','--');
    line(repmat(mean(sx(2:3)),1,2), yL, 'Color',[0.6 0.6 0.6], 'LineStyle','--');
end

for k = 1:numel(tasks)
    for r = 1:numel(regions)
        reg = regions{r};
        y = S.(S.fields{k}).(reg);
        if isempty(y), continue; end      % region not recorded in this task

        % animals with the same count would overlap, so step repeats sideways
        nudge = zeros(size(y));
        for u = unique(y)
            idx = find(y == u);
            nudge(idx) = (0:numel(idx)-1) * offsetStep;
        end

        for i = 1:numel(y)
            scatter(xCenters(k) + offset.(reg) + nudge(i), y(i), ms, 'o', ...
                'MarkerEdgeColor', colors.(reg), 'MarkerFaceColor', colors.(reg), ...
                'MarkerFaceAlpha', faceAlpha, 'LineWidth', 1.2);
        end
    end
end

h = gobjects(1, numel(regions));
for r = 1:numel(regions)
    h(r) = scatter(NaN, NaN, ms, 'o', 'MarkerEdgeColor', colors.(regions{r}), ...
                   'MarkerFaceColor', colors.(regions{r}), 'MarkerFaceAlpha', faceAlpha);
end
legend(h, regions, 'Location','northeast', 'Box','off');

xlim([-7 xCenters(end) + spacing*0.5]);
ylim(yL);
set(gca, 'XTick', xCenters, 'XTickLabel', tasks, 'XTickLabelRotation', 30, ...
         'FontSize', 14, 'Box','off', 'TickDir','out');
ylabel('Number of recorded sessions','FontSize',14);
set(gca, 'YGrid','on', 'XGrid','off', 'GridLineStyle','--', 'GridAlpha',0.3);
hold off

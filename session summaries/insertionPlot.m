%% insertionPlot.m
%  Probe insertion coordinates, one figure per task.
%
%  Each row is ONE SESSION and gives the region and [AP, ML] coordinate (mm from
%  bregma, both positive) of each probe, plus a hemisphere string whose first
%  character applies to probe 1 and second to probe 2 ('L' flips ML negative).
%  A session with one probe leaves the other entry empty. Every probe with both
%  a region and a coordinate is one dot, so a task's dots are its (session,
%  region) pairs. Dot counts are checked against the expected counts below.
%
%  Region names appear as recorded: the Delayed Reward table uses tjM1/alm/tjS1
%  and marks an absent probe 'n/a', the others use m1/m2 and an empty string.
%  canonRegion maps both onto tjM1 / ALM / tjS1.

clear; clc; close all

%% EXPECTED DOT COUNTS

expected.SimpleReward  = struct('tjM1',18, 'ALM',12);
expected.DelayedReward = struct('tjM1',21, 'ALM',15, 'tjS1',14);
expected.DoubleReward  = struct('tjM1',18, 'ALM',15);
expected.VTA           = struct('tjM1',18);
expected.Learning      = struct('tjM1',20);

nSessExp.SimpleReward = 22;  nSessExp.DelayedReward = 35;
nSessExp.DoubleReward = 26;  nSessExp.VTA = 18;  nSessExp.Learning = 20;

colors = containers.Map({'tjM1','ALM','tjS1'}, ...
    {[0 0.4470 0.7410], [0 0.6 0], [0.9290 0.6940 0.1250]});

%% ------------------------------------------------ SIMPLE REWARD TASK (22)

data = {
    '',     'm2', [],         [2.6, 1.6],  'R';
    'm1',   'm2', [2.0, 2.5], [2.6, 1.4],  'LR';
    '',     'm2', [],         [2.6, 1.4],  'R';
    '',     'm1', [],         [2.1, 2.6],  'L';
    'm2',   'm1', [2.7, 1.4], [2.2, 2.5],  'RL';

    'm1',   'm2', [2.1, 2.5], [2.5, 1.4],  'LR';
    '',     'm1', [],         [2.0, 2.7],  'L';
    '',     'm1', [],         [1.9, 2.5],  'L';
    '',     'm1', [],         [1.9, 2.6],  'L';
    '',     'm2', [],         [2.6, 1.5],  'R';

    '',     'm1', [],         [2.0, 2.6],  'L';
    '',     'm1', [],         [2.1, 2.3],  'L';
    'm2',   'm1', [2.6, 1.5], [2.0, 2.5],  'LL';
    '',     'm2', [],         [2.7, 1.4],  'LL';
    'm2',   'm1', [2.5, 1.5], [2.3, 2.7],  'LL';
    '',     'm1', [],         [2.2, 2.7],  'L';

    'm1',   'm2', [2.0, 2.8], [2.5, 1.6],  'RL';
    '',     'm1', [],         [2.1, 2.7],  'R';
    '',     'm1', [],         [2.1, 2.4],  'R';
    'm2',   'm1', [2.6, 1.4], [2.3, 2.5],  'LR';
    'm2',   'm1', [2.6, 1.7], [2.2, 2.6],  'LR';
    '',     'm1', [],         [2.1, 2.7],  'R';
};
plotInsertions(data, 'Simple Reward Task', colors, expected.SimpleReward, nSessExp.SimpleReward);

%% ----------------------------------------------- DELAYED REWARD TASK (35)

data = {
    'tjM1', 'alm', [2.2, 2.8], [2.5, 1.5], 'LL';
    'tjM1', 'alm', [2.3, 2.7], [2.6, 1.4], 'LL';
    'n/a',  'alm', [],         [2.6, 1.5], 'L';
    'tjM1', 'alm', [2.3, 2.8], [2.5, 1.5], 'LL';

    'n/a',  'tjM1', [],        [2.3, 2.8], 'R';
    'n/a',  'alm',  [],        [2.6, 1.6], 'R';
    'n/a',  'tjM1', [],        [2.4, 2.8], 'R';
    'tjM1', 'alm',  [2.5, 1.4],[2.5, 1.6], 'RR';

    'n/a',  'tjS1', [],        [0.3, 3.7], 'R';
    'tjM1', 'n/a',  [2.3, 2.4],[],         'L';
    'tjS1', 'tjM1', [0.6, 3.5],[2.2, 2.6], 'LR';
    'n/a',  'tjM1', [],        [2.3, 2.6], 'R';
    'tjS1', 'n/a',  [0.6, 3.5],[],         'L';
    'tjS1', 'n/a',  [0.6, 3.8],[],         'L';
    'tjS1', 'n/a',  [0.4, 3.6],[],         'L';

    'n/a',  'tjM1', [],        [2.4, 2.4], 'R';
    'tjM1', 'tjS1', [2.3, 2.5],[0.5, 3.6], 'RR';
    'tjM1', 'tjS1', [2.3, 2.4],[0.3, 3.4], 'RR';
    'tjS1', 'tjM1', [0.3, 3.5],[2.2, 2.6], 'RR';

    'tjS1', 'n/a',  [0.5, 3.6],[],         'L';
    'n/a',  'tjS1', [],        [0.4, 3.5], 'R';
    'n/a',  'tjS1', [],        [0.5, 3.6], 'R';
    'n/a',  'tjS1', [],        [0.4, 3.7], 'R';
    'tjS1', 'n/a',  [0.5, 3.5],[],         'L';
    'tjS1', 'n/a',  [0.4, 3.8],[],         'L';

    'n/a',  'alm',  [],        [2.6, 1.5], 'R';
    'tjM1', 'alm',  [2.2, 2.6],[2.7, 1.4], 'RR';
    'alm',  'tjM1', [2.5, 1.4],[2.3, 2.5], 'RR';
    'alm',  'tjM1', [2.7, 1.5],[2.2, 2.5], 'RR';
    'tjM1', 'n/a',  [2.4, 2.6],[],         'R';

    'alm',  'tjM1', [2.7, 1.5],[2.1, 2.5], 'LR';
    'alm',  'tjM1', [2.6, 1.5],[2.2, 2.5], 'LR';
    'alm',  'tjM1', [2.6, 1.4],[2.3, 2.5], 'LR';
    'tjM1', 'alm',  [2.3, 2.6],[2.5, 1.5], 'RL';
    'n/a',  'alm',  [],        [2.7, 1.4], 'L';
};
plotInsertions(data, 'Delayed Reward Task', colors, expected.DelayedReward, nSessExp.DelayedReward);

%% ------------------------------------------------ DOUBLE REWARD TASK (26)

data = {
    'm1',  'alm', [1.8, 2.6], [2.7, 1.6], 'LL';
    'm1',  '',    [1.9, 2.5], [],         'L';
    'm1',  'alm', [2.0, 2.6], [2.4, 1.4], 'LL';
    'm1',  'alm', [2.1, 2.6], [2.5, 1.5], 'LL';

    'm1',  '',    [1.8, 2.7], [],         'L';
    'm1',  '',    [1.9, 2.7], [],         'L';
    'm1',  '',    [2.2, 2.4], [],         'L';
    'm1',  'alm', [2.1, 2.6], [2.7, 1.4], 'RR';
    'm1',  'alm', [2.2, 2.5], [2.6, 1.5], 'RR';

    'm1',  '',    [1.8, 2.7], [],         'L';
    'm1',  'alm', [2.1, 2.6], [2.5, 1.4], 'LL';
    'm1',  'alm', [2.0, 2.6], [2.7, 1.5], 'LL';
    'm1',  '',    [2.0, 2.5], [],         'R';
    'm2',  '',    [2.6, 1.3], [],         'R';

    'm1',  '',    [1.9, 2.7], [],         'R';
    'm1',  '',    [2.1, 2.6], [],         'L';
    'm1',  '',    [2.1, 2.6], [],         'L';
    '',    'm1',  [],         [2.2, 2.6], 'R';
    '',    'm1',  [],         [2.3, 2.5], 'R';

    '',    'alm', [],         [2.5, 1.5], 'L';
    '',    'alm', [],         [2.6, 1.4], 'L';
    '',    'alm', [],         [2.6, 1.3], 'L';
    '',    'alm', [],         [2.5, 1.5], 'R';
    '',    'alm', [],         [2.6, 1.5], 'R';
    '',    'alm', [],         [2.6, 1.4], 'R';
    '',    'alm', [],         [2.5, 1.3], 'R';
};
plotInsertions(data, 'Double Reward Task', colors, expected.DoubleReward, nSessExp.DoubleReward);

%% ---------------------------------------------------- VTA REWARD TASK (18)

data = {
    'm1', '',    [2.3, 2.6], [],         'L';
    'm1', '',    [2.4, 2.6], [],         'L';
    'm1', '',    [2.1, 2.7], [],         'L';
    'm1', '',    [2.5, 2.5], [],         'L';

    'm1', '',    [1.9, 2.6], [],         'L';
    'm1', '',    [2.3, 2.4], [],         'L';
    'm1', '',    [2.3, 2.3], [],         'L';
    'm1', '',    [2.4, 2.5], [],         'L';

    '',   'm1',  [],         [1.9, 2.6], 'R';
    '',   'm1',  [],         [2.1, 2.6], 'R';
    '',   'm1',  [],         [2.1, 2.4], 'R';
    'm1', '',    [2.2, 2.6], [],         'R';
    '',   'm1',  [],         [1.9, 2.4], 'R';
    '',   'm1',  [],         [2.0, 2.5], 'R';

    '',   'm1',  [],         [2.1, 2.7], 'R';
    '',   'm1',  [],         [2.2, 2.5], 'R';

    'm1', '',    [2.0, 2.7], [],         'R';
    '',   'm1',  [],         [2.0, 2.6], 'R';
};
plotInsertions(data, 'VTA Reward Task', colors, expected.VTA, nSessExp.VTA);

%% ------------------------------------------------------ LEARNING TASK (20)

data = {
    'm1',  '',    [2.1, 2.6], [],         'L';
    '',    'm1',  [],         [1.9, 2.7], 'L';
    'm1',  '',    [2.2, 2.5], [],         'L';
    'm1',  '',    [2.1, 2.6], [],         'L';
    '',    'm1',  [],         [2.0, 2.5], 'L';

    'm1',  '',    [2.0, 2.6], [],         'R';
    '',    'm1',  [],         [2.1, 2.6], 'R';
    'm1',  '',    [2.0, 2.4], [],         'R';
    '',    'm1',  [],         [1.9, 2.6], 'R';
    '',    'm1',  [],         [2.0, 2.6], 'R';

    '',    'm1',  [],         [2.1, 2.5], 'L';
    'm1',  '',    [1.9, 2.6], [],         'L';
    'm1',  '',    [2.1, 2.4], [],         'L';
    'm1',  '',    [2.2, 2.6], [],         'L';
    'm1',  '',    [2.1, 2.4], [],         'L';

    '',    'm1',  [],         [1.9, 2.7], 'R';
    'm1',  '',    [2.0, 2.5], [],         'R';
    '',    'm1',  [],         [2.2, 2.5], 'R';
    'm1',  '',    [2.1, 2.4], [],         'R';
    'm1',  '',    [2.0, 2.6], [],         'R';
};
plotInsertions(data, 'Learning', colors, expected.Learning, nSessExp.Learning);


%% LOCAL FUNCTIONS

function plotInsertions(data, taskName, colors, expected, nSessExp)
% One figure for this task, and a check of the dot counts against expected.

    nRow = size(data,1);
    assert(nRow == nSessExp, ...
        ['%s: the table has %d rows but %d sessions are expected ' ...
         '(one row per session).'], taskName, nRow, nSessExp);

    X = []; Y = []; C = []; regionOf = {};

    for i = 1:nRow
        names  = {strtrim(data{i,1}), strtrim(data{i,2})};
        coords = {data{i,3},          data{i,4}};
        side   = strtrim(data{i,5});

        for pr = 1:2
            reg = canonRegion(names{pr});
            if isempty(reg) || isempty(coords{pr}), continue; end

            ap = coords{pr}(1);
            ml = coords{pr}(2);

            % first character of the side string is probe 1, second is probe 2;
            % a one-character string applies to whichever probe is present
            if numel(side) >= pr
                thisSide = upper(side(pr));
            elseif ~isempty(side)
                thisSide = upper(side(1));
            else
                thisSide = 'R';
            end
            if thisSide == 'L', ml = -abs(ml); else, ml = abs(ml); end

            X(end+1,1) = ml;             %#ok<AGROW>
            Y(end+1,1) = ap;             %#ok<AGROW>
            C(end+1,:) = colors(reg);    %#ok<AGROW>
            regionOf{end+1,1} = reg;     %#ok<AGROW>
        end
    end

    figure('Color','w'); hold on
    scatter(X, Y, 100, C, 'LineWidth', 3, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.4);
    plot([-1 0], [0.5 0.5], 'k', 'LineWidth', 2);
    text(-0.5, 0.35, '1 mm', 'HorizontalAlignment', 'center');
    xline(-5:1:5);  yline(-5:1:5);
    xlim([-5 5]);   ylim([-1 4]);
    xlabel('ML coordinate (mm)');
    ylabel('AP coordinate (mm)');
    title(taskName);
    hold off
    set(gcf, 'Position', [1100 400 300 300]);

    % ---- check the dots against the expected counts ----
    fprintf('\n%s | %d sessions, %d probe insertions\n', taskName, nRow, numel(X));
    regs = fieldnames(expected);
    bad = false;
    for r = 1:numel(regs)
        got  = sum(strcmp(regionOf, regs{r}));
        want = expected.(regs{r});
        if got ~= want, bad = true; end
        fprintf('  %-5s %3d dots   expected %3d   %s\n', regs{r}, got, want, ...
                ternary(got == want, 'ok', '<-- MISMATCH'));
    end
    seen = unique(regionOf);
    for k = 1:numel(seen)
        if ~isfield(expected, seen{k})
            fprintf('  %-5s %3d dots   NOT EXPECTED in this task\n', ...
                    seen{k}, sum(strcmp(regionOf, seen{k})));
            bad = true;
        end
    end
    if bad
        warning('%s: probe counts do not match the expected counts.', taskName);
    end
end


function reg = canonRegion(name)
% Map every spelling used in the tables onto tjM1 / ALM / tjS1.
% Returns '' for an absent probe.
    switch lower(strtrim(name))
        case {'', 'n/a'},      reg = '';
        case {'m1', 'tjm1'},   reg = 'tjM1';
        case {'m2', 'alm'},    reg = 'ALM';
        case {'s1', 'tjs1'},   reg = 'tjS1';
        otherwise
            error('canonRegion: unrecognized region name ''%s''.', name);
    end
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

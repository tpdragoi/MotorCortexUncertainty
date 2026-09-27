%% F3H_FS4F_crossTask.m
%  Simple Reward Task vs Double Reward Task, decoding index, for the tongue
%  (Fig. 3H) and the jaw (Fig. S4F).
%
%  Each of the four decoding scripts writes one value per session per region
%  when it finishes: that session's mean decoding index on C1 trials. This file
%  runs the four in order and then compares them, so the two halves of each
%  comparison are guaranteed to come from the same run rather than from
%  whatever happened to be left in tempdir.
%
%  The comparison is an unpaired one-tailed rank-sum, one test per region,
%  uncorrected. The sessions are different recordings in the two tasks, so
%  there is nothing to pair.
%
%  READS  the four summary files the decoding scripts write (see PATHS below)
%  Run the whole file.

clear; clc

%% RUN THE FOUR DECODING SCRIPTS
% Set this to false to compare summaries already on disk, which is what you
% want while iterating on the statistics -- the four scripts take a long time.
%
% Each decoding script begins with `clear`, so nothing defined above this
% block survives it. That is why every variable this file needs is built
% after the last script has finished, and why the four calls are written out
% rather than looped over a list.

RUN_SCRIPTS = true;

if RUN_SCRIPTS
    F1H_decodeTongue       % Simple Reward, tongue   -> Fig. 1H
    F3H_decodeTongue       % Double Reward, tongue   -> Fig. 3H
    FS4B_decodeJaw         % Simple Reward, jaw      -> Fig. S4B
    FS4F_decodeJaw         % Double Reward, jaw      -> Fig. S4F
end

%% PATHS
% These must match spec.summaryFile in the four scripts above. If you change
% one there, change it here.

simpleTongue = fullfile(tempdir, 'decodeSummary_tongue_r1_clean.mat');
doubleTongue = fullfile(tempdir, 'decodeSummary_tongue_r16_clean.mat');
simpleJaw    = fullfile(tempdir, 'decodeSummary_jaw_r1_clean.mat');
doubleJaw    = fullfile(tempdir, 'decodeSummary_jaw_r16_clean.mat');

%% SETTINGS

alpha = 0.05;

% Tail describes the Simple Reward Task relative to the Double Reward Task.
% 'left' asks whether decoding is LOWER in the Simple Reward Task, which is the
% direction Fig. 3H reports.
tail = 'left';

%% THE TWO COMPARISONS

fprintf('\nCROSS-TASK DECODING INDEX, Simple Reward vs Double Reward\n');

tongue = crossTaskCompare(simpleTongue, doubleTongue, ...
                          'Simple Reward, tongue', 'Double Reward, tongue', ...
                          alpha, tail);

jaw    = crossTaskCompare(simpleJaw, doubleJaw, ...
                          'Simple Reward, jaw', 'Double Reward, jaw', ...
                          alpha, tail);

%% ONE LINE PER PANEL, READY FOR THE LEGEND

fprintf('%s\n', repmat('=', 1, 78));
fprintf('FOR THE FIGURE LEGENDS\n');
fprintf('%s\n', repmat('-', 1, 78));
legendLine('Fig. 3H  (tongue)', tongue, alpha);
legendLine('Fig. S4F (jaw)',    jaw,    alpha);
fprintf('%s\n\n', repmat('=', 1, 78));


%% LOCAL FUNCTIONS

function legendLine(panel, res, alpha)
% One sentence per region, in the form a legend needs it.
    if isempty(res)
        fprintf('%-18s no region could be tested\n', panel);
        return
    end
    for k = 1:numel(res)
        r = res(k);
        if r.p < alpha, verdict = 'p < 0.05'; else, verdict = 'n.s.'; end
        fprintf('%-18s %-5s : Simple n = %d vs Double n = %d, ', ...
            panel, r.region, r.nA, r.nB);
        fprintf('one-tailed rank-sum p = %.4f (%s)\n', r.p, verdict);
        panel = '';   % the panel name only on its first line
    end
end

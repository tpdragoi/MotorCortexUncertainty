function out = crossTaskCompare(fileA, fileB, labelA, labelB, alpha, tail)
% crossTaskCompare  Compare the per-session decoding index of one task against
% another, region by region.
%
%   out = crossTaskCompare(fileA, fileB, labelA, labelB, alpha, tail)
%
% fileA and fileB are the summary .mat files the decoding scripts write to
% spec.summaryFile. Each holds one value per session per region, and that value
% is a DELTA, not a decoding index: the session's mean change in decoding index
% from contact 1, averaged over the contacts in spec.groupSummaryLicks (3 to 8
% in the tongue and jaw scripts). It is normally negative, because decoding
% declines over a bout. The sessions in the two files are different recordings,
% so the comparison is an UNPAIRED rank-sum, one test per region, uncorrected.
%
% tail is passed to ranksum and describes A relative to B. 'left' asks whether
% A is lower than B, which is the Fig. 3H / S4F question: is decoding weaker in
% the Simple Reward Task than in the Double Reward Task.
%
% Returns a struct array, one entry per region compared, with the region label,
% both n, both means, the one-tailed p, and the smallest p the test
% could have returned at those sample sizes.
%
% A rank-sum on n_A vs n_B sessions cannot return a one-tailed p below
% 1/nchoosek(n_A+n_B, n_A), so that floor is reported next to every p. With
% small session counts the floor is often close to alpha, and a p sitting on it
% means only that the two groups separated completely.

    if nargin < 5 || isempty(alpha), alpha = 0.05;    end
    if nargin < 6 || isempty(tail),  tail  = 'left';  end

    A = loadSummary(fileA, labelA);
    B = loadSummary(fileB, labelB);

    nG = min(numel(A.summaryByGroup), numel(B.summaryByGroup));
    out = struct('region', {}, 'nA', {}, 'nB', {}, 'meanA', {}, 'meanB', {}, ...
                 'p', {}, 'floor', {});

    fprintf('\n%s\n', repmat('=', 1, 78));
    fprintf('%s  vs  %s\n', A.header, B.header);
    fprintf('one value per session: mean change in decoding index from contact 1%s\n', ...
        contactsPhrase(A, B));
    fprintf('unpaired rank-sum, tail = %s\n', tail);
    fprintf('%s\n', repmat('-', 1, 78));
    fprintf('%-8s %6s %10s %6s %10s %12s %10s\n', ...
        'region', 'n A', 'mean A', 'n B', 'mean B', 'p one-tail', 'floor');

    for g = 1:nG
        region = regionName(A, B, g);

        a = A.summaryByGroup{g};  a = a(isfinite(a));
        b = B.summaryByGroup{g};  b = b(isfinite(b));

        if numel(a) < 2 || numel(b) < 2
            fprintf('%-8s %6d %10s %6d %10s %12s %10s\n', ...
                region, numel(a), '--', numel(b), '--', 'not run', '--');
            continue
        end

        p1 = ranksum(a, b, 'tail', tail);
        pFloor = 1 / nchoosek(numel(a) + numel(b), numel(a));

        fprintf('%-8s %6d %10.4f %6d %10.4f %12.4f%s %10.4f\n', ...
            region, numel(a), mean(a), numel(b), mean(b), ...
            p1, pStars(p1, [alpha 0.01 0.001]), pFloor);

        out(end+1) = struct('region', region, ...                       %#ok<AGROW>
            'nA', numel(a), 'nB', numel(b), ...
            'meanA', mean(a), 'meanB', mean(b), ...
            'p', p1, 'floor', pFloor);
    end

    onFloor = arrayfun(@(r) abs(r.p - r.floor) < 1e-12, out);
    if any(onFloor)
        fprintf('%s\n', repmat('-', 1, 78));
        fprintf('NOTE: %s sits on the exact floor, so the two groups separated\n', ...
            strjoin({out(onFloor).region}, ', '));
        fprintf('      completely and no smaller p was available at this n.\n');
    end
    fprintf('%s\n\n', repmat('=', 1, 78));
end


function S = loadSummary(f, label)
% Read one summary file and fail with a message that names the script to run.
    assert(exist(f, 'file') == 2, ...
        ['No summary for %s at\n  %s\nRun that task''s decoding script first; ' ...
         'it writes this file when it finishes.'], label, f);

    S = load(f);
    assert(isfield(S, 'summaryByGroup'), ...
        '%s is not a decoding summary (no summaryByGroup).', f);

    if isfield(S, 'figRef') && isfield(S, 'specName')
        S.header = sprintf('%s [%s, %s]', label, S.figRef, S.specName);
    else
        S.header = label;
    end
end


function name = regionName(A, B, g)
% Region g's label, checked against the other file when both recorded one.
    hasA = isfield(A, 'groupLabels') && numel(A.groupLabels) >= g;
    hasB = isfield(B, 'groupLabels') && numel(B.groupLabels) >= g;

    if hasA && hasB
        assert(strcmp(A.groupLabels{g}, B.groupLabels{g}), ...
            ['Region %d is %s in one summary and %s in the other, so the two ' ...
             'files do not list regions in the same order and this comparison ' ...
             'would pair the wrong ones.'], g, A.groupLabels{g}, B.groupLabels{g});
        name = A.groupLabels{g};
    elseif hasA
        name = A.groupLabels{g};
    elseif hasB
        name = B.groupLabels{g};
    else
        name = sprintf('group %d', g);   % summary without region labels
    end
end


function phrase = contactsPhrase(A, B)
% Name the contacts the delta was averaged over, when both summaries record
% them and they agree. Returns '' if either summary lacks summaryLicks.
    phrase = '';
    if ~isfield(A, 'summaryLicks') || ~isfield(B, 'summaryLicks'), return; end
    if ~isequal(A.summaryLicks, B.summaryLicks)
        phrase = sprintf(' (WARNING: averaged over %s in one file and %s in the other)', ...
            mat2str(A.summaryLicks), mat2str(B.summaryLicks));
        return
    end
    phrase = sprintf(', averaged over contacts %s', mat2str(A.summaryLicks));
end

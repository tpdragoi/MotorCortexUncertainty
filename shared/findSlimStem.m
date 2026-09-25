function stem = findSlimStem(dataDir, anmHint, dateStr)
% FINDSLIMSTEM  File stem '<ANM>_<YYYY_MM_DD>' of an exported session in dataDir.
%
%   stem = findSlimStem(fullfile(root,'Data','Learning'), 'TD3l', '2025-02-03')
%
% The stem uses obj.pth.anm / obj.pth.dt exactly as the pipeline stored them, and
% those do not always match the name in a script's session list:
%   - loader names: loadTD3l_many loads animal 'TDl3', loadTD27_many 'TD27d', ...
%   - one date was stored with a typo: TDl3 '2025-02-044'
% So: every *_obj.mat whose date part STARTS WITH the requested date is a
% candidate; if more than one animal recorded that day, the one whose name matches
% anmHint (case-insensitive, or with the same letters in another order, e.g.
% TD3l = TDl3) is used, then one whose name starts with anmHint (TD27 -> TD27d).
% Anything ambiguous or missing is an error, never a guess.
    d   = strrep(dateStr, '-', '_');
    L   = dir(fullfile(dataDir, '*_obj.mat'));
    stems = regexprep({L.name}, '_obj\.mat$', '');
    tok = regexp(stems, '^(.*?)_(\d{4}_\d{2}_\d+)$', 'tokens', 'once');
    ok  = ~cellfun(@isempty, tok);
    stems = stems(ok);  tok = tok(ok);
    anms  = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
    dates = cellfun(@(t) t{2}, tok, 'UniformOutput', false);

    cand = find(startsWith(dates, d));
    if isempty(cand)
        error('findSlimStem: no exported session on %s in %s.', dateStr, dataDir);
    end
    if numel(cand) > 1 && ~isempty(anmHint)
        a  = anms(cand);
        m1 = strcmpi(a, anmHint) | strcmp(cellfun(@(s) sort(lower(s)), a, 'UniformOutput', false), sort(lower(anmHint)));
        if any(m1), cand = cand(m1);
        else
            m2 = startsWith(lower(a), lower(anmHint));
            if any(m2), cand = cand(m2); end
        end
    end
    if numel(cand) > 1
        % exact date wins over a longer (typo) date string
        ex = strcmp(dates(cand), d);
        if nnz(ex) == 1, cand = cand(ex); end
    end
    if numel(cand) ~= 1
        error('findSlimStem: %d sessions match %s %s in %s: %s', numel(cand), anmHint, dateStr, ...
            dataDir, strjoin(stems(cand), ', '));
    end
    stem = stems{cand};
    if ~isempty(anmHint) && ~strcmpi(anms{cand}, anmHint)
        fprintf('  [stem] %s %s -> %s\n', anmHint, dateStr, stem);
    end
end

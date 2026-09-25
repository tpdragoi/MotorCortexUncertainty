function [obj, params, kin] = slimToLegacy(meta, params)
% SLIMTOLEGACY  Rebuild what loadSessionData + loadMotionEnergy + getKinematics
% returned for THESE params, from the exported slim files.
%
%   meta = slimMeta('loadTD26_neur', '2025-08-06');
%   [obj, params, kin] = slimToLegacy(meta, params);
%   [obj, params]      = slimToLegacy(meta, params);   % kin file not read
%
% meta may be an array (all_meta): obj / params / kin are then struct arrays,
% one element per session, as loadSessionData returned.
%
% replaces
%   [obj, params] = loadSessionData(meta, params, params.behav_only);
%   me  = loadMotionEnergy(obj, meta, params, datapth);
%   kin = getKinematics(obj, me, params);
%
% The export was made deliberately permissive (smooth = 1, lowFR = 1e-10, a
% -2.5..32 s window, three sample rates, the 10 standard conditions). This
% function applies each script's OWN settings on top, by the same code path the
% pipeline used (DataLoadingScripts/getSeq, removeLowFRClusters, findTrials,
% findClusters and utils/mySmooth):
%
%   params.alignEvent   'firstLick' -> objFL/kinFL,  'goCue' -> objGC/kinGC
%   params.dt           picks the exported rate: 1/300, 1/200 or 1/100 only
%   params.tmin/tmax    cropped to exactly the bins getSeq would have made
%   params.smooth       mySmooth(trialdat, params.smooth, params.bctype) AFTER
%                       cropping, so the reflect boundary sits at tmin exactly
%                       as it did in getSeq
%   params.condition    params.trialid = findTrials(obj, params.condition)
%   psth                per condition: mean over that condition's trials of the
%                       smoothed single-trial rates (== getSeq; smoothing is
%                       linear); a condition with no trials is all zeros, as in
%                       getSeq
%   params.lowFR        keep units whose MEAN RATE OVER ALL TRIALS of the
%                       session, in a goCue-aligned -2 .. 4 s window, is
%                       > lowFR. This ignores params.condition, params.alignEvent,
%                       params.tmin/tmax and params.smooth, so the same units
%                       survive in every figure. The goCue file is read for this
%                       even when the script aligns to firstLick.
%   params.quality      same list as the export -> exported units used directly.
%                       Different list -> findClusters on objFL.clu; units the
%                       export did not keep are re-binned from their spike
%                       times exactly as alignSpikes + getSeq do.
%
% OUTPUT SHAPE. obj.time is a row vector and obj.psth / obj.trialdat are
% double, as before. params.cluid has the per-probe shape the old scripts
% index into: a cell with one slot per probe, EMPTY for a probe the export
% dropped because no group map used it (so size(params.cluid{1,1},1) is 0 when
% probe 1 was dropped, and "probe 2 = units after cluid{1,1}" still holds).
% Single-probe sessions keep cluid as a plain array, as before. obj.probeOfUnit
% gives the original probe of every unit directly. The unit axis is sorted by
% probe, so "probe 1 = the first size(params.cluid{1,1},1) units, probe 2 = the
% rest" picks out exactly the units obj.probeOfUnit == p picks out, on every
% session including the ones the exporter treated as single-probe.
%
% kin.dat / kin.featLeg: getKinematics output at this rate and window.
% Kinematics do not depend on params.smooth or lowFR.
%
% KNOWN, SMALL DIFFERENCES from a fresh loadSessionData run:
%  * Spikes lying exactly on a bin edge can fall on the other side, because the
%    export's bin edges start at -2.5 s, not at params.tmin (floating point).
%  * getKinematics fills tongue x/y positions where the tongue is out of view
%    with the session mean start position computed within the load window; the
%    export's window is longer, so those filled values can differ slightly.
%  * If findTrials returns trial numbers beyond obj.bp.Ntrials (sessions with
%    obj.bp.fidx), those trials have no single-trial data and are left out of
%    the psth, as getSeq's trialdat also leaves them out.
%
% Uses the pipeline's own mySmooth, findTrials and findClusters (and the
% getStructVarNames / patternMatchCellArray they call) from
% <MATLAB Codes _ v2>\shared\pipelineCopies -- byte-identical copies of the
% uninstructedMovements_v2-main files, so the smoothing is exactly the
% pipeline's and nothing reads that folder.

    for f = {'mySmooth','findTrials','findClusters'}
        assert(exist(f{1}, 'file') == 2, ['slimToLegacy needs %s on the MATLAB path: ' ...
            'addpath(fullfile(<MATLAB Codes _ v2>, ''shared'', ''pipelineCopies'')).'], f{1});
    end
    if numel(params) > 1, params = params(1); end        % a struct array from a previous call
    if numel(meta) > 1
        % several sessions at once, as loadSessionData(all_meta, ...) did:
        % struct arrays out, one element per session
        O = cell(1, numel(meta));  P = O;  Kc = O;
        for i = 1:numel(meta)
            if nargout >= 3
                [O{i}, P{i}, Kc{i}] = slimToLegacy(meta(i), params);
            else
                [O{i}, P{i}] = slimToLegacy(meta(i), params);
            end
        end
        obj = [O{:}];  params = [P{:}];
        if nargout >= 3, kin = [Kc{:}]; end
        return
    end
    if isfield(params, 'timeWarp')
        assert(~params.timeWarp, 'slimToLegacy: time-warped data were not exported.');
    end
    assert(ismember(params.alignEvent, {'firstLick','goCue'}), ...
        'slimToLegacy: only firstLick and goCue alignments were exported, not %s.', params.alignEvent);

    rate = round(1/params.dt);
    assert(abs(1/rate - params.dt) < 1e-12, 'slimToLegacy: params.dt = %g is not 1/integer.', params.dt);

    if ~isfield(params, 'bctype'), bctype = 'none'; else, bctype = params.bctype; end   % as processData
    behavOnly = isfield(params, 'behav_only') && params.behav_only;

    % ---- read one alignment at one rate -------------------------------------
    dirPath = fileparts(meta.objFile);
    if nargout >= 3
        [S, K] = loadSlimSession(dirPath, meta.stem, params.alignEvent, rate);
    else
        S = loadSlimSession(dirPath, meta.stem, params.alignEvent, rate);   % kin file not read
    end
    behavOnly = behavOnly || S.behavOnly;

    % ---- crop to the script's bins ------------------------------------------
    edges = params.tmin : params.dt : params.tmax;       % exactly as getSeq
    tNew  = edges(1:end-1) + params.dt/2;
    n     = numel(tNew);
    t0    = S.time(1) - params.dt/2;                      % first exported bin edge
    i0    = round((params.tmin - t0) / params.dt) + 1;
    rows  = i0 : i0 + n - 1;
    assert(i0 >= 1 && rows(end) <= numel(S.time), ...
        ['slimToLegacy: window %g..%g s is outside the exported window %.3f..%.3f s for %s. ' ...
         'Narrow params.tmin/tmax.'], params.tmin, params.tmax, t0, ...
        S.time(end) + params.dt/2, meta.stem);
    assert(max(abs(S.time(rows)' - tNew)) < params.dt * 1e-3, ...
        'slimToLegacy: %s: params.tmin %g does not fall on the exported bin grid.', ...
        meta.stem, params.tmin);

    % ---- behaviour ------------------------------------------------------------
    obj = struct();
    obj.pth  = S.pth;
    obj.bp   = S.bp;
    obj.time = tNew;                                      % row, as loadSessionData
    trialid  = findTrials(struct('bp', S.bp, 'pth', S.pth), params.condition);

    if nargout >= 3
        kin = struct();
        kin.featLeg = K.featLeg;
        kin.dat     = double(K.dat(rows, :, :));
        kin.time    = tNew;
        clear K
    end

    params.trialid = trialid;
    params.probe   = S.keptProbes;

    if behavOnly
        params.cluid     = {};
        obj.probeOfUnit  = zeros(0,1);
        return
    end

    % ---- units: raw (unsmoothed) single-trial rates in the script's window ----
    nTr = size(S.trialdat, 3);
    [X, cluidK, probeOf] = selectUnits(S, meta, params, rows, edges, nTr);

    % ---- the rate the lowFR criterion tests ----------------------------------
    % [minFR] Not the psth averaged over params.condition, and not the script's
    % own window: each unit's mean rate over ALL trials of the session in a
    % goCue-aligned -2 .. 4 s window. Every script therefore keeps the same
    % units, whatever it aligns to and whichever conditions it asks for.
    if strcmp(params.alignEvent, 'goCue')
        [frMean, frIds, frProbes] = lowFRrate(S, meta, params);
        S.trialdat = [];  S.psth = [];                 % release the exported arrays
    else
        S.trialdat = [];  S.psth = [];                 % release before the goCue copy
        [frMean, frIds, frProbes] = lowFRrate([], meta, params);
    end
    assert(isequal(frIds(:), cluidK(:)) && isequal(frProbes(:), probeOf(:)), ...
        ['slimToLegacy: %s: the goCue export lists different units than the %s ' ...
         'export, so the firing-rate criterion cannot be matched to them. ' ...
         'Re-export this session.'], meta.stem, params.alignEvent);

    % [probe] Keep the unit axis in probe order, so the positional split the old
    % scripts use is the same set of units as selecting on obj.probeOfUnit.
    if ~issorted(probeOf)
        [probeOf, ord] = sort(probeOf(:), 'ascend');   % stable: order within a probe kept
        X = X(:, ord, :);  cluidK = cluidK(ord);  frMean = frMean(ord);
    else
        probeOf = probeOf(:);
    end

    % ---- smooth exactly as getSeq: mySmooth on each unit's time x trial block,
    % in place (one unit at a time keeps peak memory at one copy of trialdat)
    nU = size(X, 2);
    for u = 1:nU
        X(:, u, :) = reshape(mySmooth(reshape(X(:, u, :), n, nTr), params.smooth, bctype), n, 1, nTr);
    end

    % ---- psth per condition (getSeq) ------------------------------------------
    nC   = numel(trialid);
    psth = zeros(n, nU, nC);
    nBeyond = 0;
    for c = 1:nC
        trix = trialid{c};
        if isempty(trix), continue; end
        ok = trix(trix >= 1 & trix <= nTr);
        nBeyond = nBeyond + numel(trix) - numel(ok);
        if ~isempty(ok)
            psth(:, :, c) = sum(X(:, :, ok), 3) ./ numel(trix);
        end
    end
    if nBeyond > 0
        warning('slimToLegacy:trialsBeyondNtrials', ...
            '%s: %d trialid entries exceed the %d single trials and were left out of psth.', ...
            meta.stem, nBeyond, nTr);
    end

    % ---- removeLowFRClusters ------------------------------------------------
    use = frMean(:) > params.lowFR;                   % [minFR] goCue -2..4 s, all trials
    obj.psth        = psth(:, use, :);
    obj.trialdat    = X(:, use, :);
    obj.probeOfUnit = probeOf(use);
    cluidK          = cluidK(use);

    % ---- cluid in the shape the old code indexes into -------------------------
    if iscell(S.cluid)
        nSlots = max([2, S.keptProbes(:)']);
        cl = cell(1, nSlots);
        for p = 1:nSlots, cl{p} = zeros(0,1); end
        for p = S.keptProbes(:)'
            cl{p} = cluidK(obj.probeOfUnit == p);
            cl{p} = cl{p}(:);
        end
        params.cluid = cl;
    else
        params.cluid = cluidK(:);
    end
end


function [raw, cluidOut, probeOut] = selectUnits(S, meta, params, rows, edges, nTr, meanOnly)
% Raw single-trial rates (time x unit x trial) for the units params.quality
% selects, on the script's bins, with their cluster ids and original probes.
%
% meanOnly = true returns, instead of that block, one number per unit: its mean
% rate over every bin and every trial of the window. Same units, same order. It
% is built a unit at a time so a second full trialdat is never held in memory.
    if nargin < 7 || isempty(meanOnly), meanOnly = false; end
    exportQ = S.params.quality;
    sameQ   = isequal(sort(exportQ(:)), sort(params.quality(:)));

    probes = S.keptProbes(:)';
    if isempty(probes), probes = unique(S.probeOfUnit(:))'; end

    % Cluster ids in the SAME ORDER as the columns of trialdat. Usually cluid has
    % one cell per kept probe, so cell k belongs to probes(k). A session the
    % exporter treated as single-probe is different: it keeps the loaded probe's
    % cluid as it was (which can itself be a per-probe cell) and labels EVERY unit
    % with the mapped probe number, so cluid{k} does not line up with probes(k)
    % there (e.g. R16 TD25d 2025-07-14). Concatenating the cells is right in both
    % layouts, because the unit axis follows the cells in order.
    cluidAll  = allCluid(S.cluid);
    nUnitsExp = size(S.trialdat, 2);
    assert(numel(cluidAll) == nUnitsExp, ...
        ['%s: the export holds %d units but %d cluster ids -- the obj file is ' ...
         'inconsistent. Re-export this session.'], meta.stem, nUnitsExp, numel(cluidAll));

    if sameQ
        if meanOnly
            nU  = size(S.trialdat, 2);
            raw = zeros(nU, 1);
            for u = 1:nU
                raw(u) = mean(mean(double(S.trialdat(rows, u, :)), 3), 1);
            end
        else
            raw = double(S.trialdat(rows, :, :));
        end
        cluidOut = cluidAll;
        probeOut = S.probeOfUnit(:);
        return
    end

    % Different quality list: decide unit by unit from the full cluster list.
    % The unit axis follows the cluid cells in order, and clu has one entry per
    % cluid cell, so walk the cells rather than the kept probes -- the two are the
    % same list except on the sessions the exporter treated as single-probe.
    clu = getClu(S, meta);
    if iscell(S.cluid), cells = S.cluid(:)'; else, cells = {S.cluid}; end
    counts = cellfun(@numel, cells);
    offs   = [0 cumsum(counts)];
    n   = numel(edges) - 1;
    assert(isfield(S.bp.ev, params.alignEvent), ['%s: bp.ev.%s is missing, so units cannot ' ...
        'be re-binned to that event.'], meta.stem, params.alignEvent);
    ev  = S.bp.ev.(params.alignEvent);
    if meanOnly, raw = zeros(0, 1); else, raw = zeros(n, 0, nTr); end
    cluidOut = zeros(0,1);  probeOut = zeros(0,1);
    nRebuilt = 0;

    for k = 1:numel(cells)
        cols = offs(k)+1 : offs(k+1);               % this cell's columns of trialdat
        if isempty(cols), continue; end
        p    = S.probeOfUnit(cols(1));              % the probe the export labelled them with
        cK   = cluSegment(clu, k, numel(cells), meta.stem);
        have = cluidAll(cols);                      % their cluster ids, same order

        want = findClusters({cK.quality}', params.quality);
        want = want(:);
        if meanOnly
            Mk = zeros(numel(want), 1);
            for u = 1:numel(want)
                j = find(have == want(u), 1);
                if ~isempty(j)
                    col = double(S.trialdat(rows, cols(j), :));
                else
                    col = rebin(cK(want(u)), ev, edges, nTr);
                    nRebuilt = nRebuilt + 1;
                end
                Mk(u) = mean(mean(col, 3), 1);
            end
            raw = [raw; Mk];                       %#ok<AGROW>
        else
            Rk = zeros(n, numel(want), nTr);
            for u = 1:numel(want)
                j = find(have == want(u), 1);
                if ~isempty(j)
                    Rk(:, u, :) = double(S.trialdat(rows, cols(j), :));
                else
                    Rk(:, u, :) = rebin(cK(want(u)), ev, edges, nTr);
                    nRebuilt = nRebuilt + 1;
                end
            end
            raw = cat(2, raw, Rk);
        end
        cluidOut = [cluidOut; want];               %#ok<AGROW>
        probeOut = [probeOut; repmat(p, numel(want), 1)]; %#ok<AGROW>
    end
    if nRebuilt > 0 && ~meanOnly
        fprintf('  [slim] %s: %d units re-binned from spike times for quality {%s}\n', ...
            meta.stem, nRebuilt, strjoin(params.quality, ','));
    end
end


function [fr, ids, probes] = lowFRrate(S, meta, params)
% Each unit's mean firing rate over ALL trials of the session, in a goCue-aligned
% window from FRWIN(1) to FRWIN(2) seconds. This is the quantity params.lowFR is
% compared with, and it depends on nothing else the script sets -- not the
% alignment, not tmin/tmax, not smooth, not condition -- so a unit kept for one
% figure is kept for all of them.
%
% S is the already-loaded goCue session when the script aligns to goCue, and []
% when it does not, in which case the goCue file is read here.
    FRWIN = [-2 4];
    dt    = params.dt;

    if isempty(S)
        G = loadSlimSession(fileparts(meta.objFile), meta.stem, 'goCue', round(1/dt));
    else
        G = S;
    end

    edgesFR = FRWIN(1) : dt : FRWIN(2);
    nFR     = numel(edgesFR) - 1;
    t0      = G.time(1) - dt/2;
    i0      = round((FRWIN(1) - t0) / dt) + 1;
    rowsFR  = i0 : i0 + nFR - 1;
    assert(i0 >= 1 && rowsFR(end) <= numel(G.time), ...
        ['slimToLegacy: %s: the goCue %g .. %g s window the lowFR criterion uses ' ...
         'falls outside the exported window.'], meta.stem, FRWIN(1), FRWIN(2));

    pGC = params;
    pGC.alignEvent = 'goCue';
    [fr, ids, probes] = selectUnits(G, meta, pGC, rowsFR, edgesFR, size(G.trialdat, 3), true);
    fr = fr(:);
end


function R = rebin(c, ev, edges, nTr)
% alignSpikes (no time warp) + getSeq's histc, for one cluster, smooth = 1.
    ev  = ev(:);                                   % bp.ev.* may be stored as a row
    tr  = c.trial(:);
    ok  = tr >= 1 & tr <= min(nTr, numel(ev));
    tr  = tr(ok);
    tm  = c.trialtm(:);
    x   = tm(ok) - ev(tr);                         % both columns: element-wise
    [~, bin] = histc(x, edges);                    %#ok<HISTC> same binning as getSeq
    n   = numel(edges) - 1;
    keep = bin >= 1 & bin <= n;                   % histc's last slot (== edges(end)) is dropped by getSeq
    dt  = edges(2) - edges(1);
    R   = accumarray([bin(keep), tr(keep)], 1, [n, nTr]) ./ dt;
    R   = reshape(R, n, 1, nTr);
end


function clu = getClu(S, meta)
% objFL carries clu; for goCue-aligned loads read it from objFL.
    if isfield(S, 'clu') && ~isempty(S.clu)
        clu = S.clu;
        return
    end
    % current layout: clu on objFL; earlier layout (6 R1 sessions): on objFL300
    w = whos('-file', meta.objFile);  have = {w.name};
    extra = have(strncmp(have, 'objFL', 5) & ~strcmp(have, 'objFL'));
    extra = fliplr(sort(extra(:)'));                 % objFL300 before objFL200/objFL100
    cand  = [{'objFL'}, extra];
    clu = [];
    for v = cand
        if ~any(strcmp(have, v{1})), continue; end
        C = load(meta.objFile, v{1});
        if isfield(C.(v{1}), 'clu') && ~isempty(C.(v{1}).clu), clu = C.(v{1}).clu; break; end
    end
    assert(~isempty(clu), ['%s has no clu, so units with a different quality list ' ...
        'cannot be rebuilt.'], meta.stem);
end


function ids = allCluid(cluid)
% Cluster ids in the order of the exported unit axis, whatever shape cluid has:
% a cell (one entry per probe the export kept) or a plain array.
    if iscell(cluid)
        parts = cellfun(@(c) c(:), cluid(:), 'UniformOutput', false);
        ids   = vertcat(parts{:});
    else
        ids = cluid(:);
    end
    ids = ids(:);
end


function c = cluidForProbe(cluid, k)
    if iscell(cluid), c = cluid{k}(:); else, c = cluid(:); end
end


function cK = cluSegment(clu, k, nSeg, stem)
% The clu entry that goes with cluid cell k. clu is a cell with one entry per
% cluid cell on two-probe sessions, and the loaded probe's own cluster list on
% single-probe ones.
    if ~iscell(clu), cK = clu; return; end
    if numel(clu) == nSeg
        cK = clu{k};
    elseif numel(clu) == 1
        cK = clu{1};
    else
        error(['%s: clu has %d entries for %d cluid cells, so units cannot be ' ...
               'matched to clusters for a different quality list.'], stem, numel(clu), nSeg);
    end
end

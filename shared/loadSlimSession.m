function [obj, kin] = loadSlimSession(dirPath, stem, alignEvent, rate)
% LOADSLIMSESSION  One alignment of one exported session, at one sample rate,
% as a plain obj/kin pair.
%
%   [obj, kin] = loadSlimSession(dirPath, stem, alignEvent, rate)
%
%   dirPath     folder holding the files, e.g. fullfile('Data','R1')
%   stem        '<ANM>_<YYYY_MM_DD>', e.g. 'TD26d_2025_08_06'
%   alignEvent  'firstLick' (or 'FL') | 'goCue' (or 'GC')
%   rate        sample rate in Hz: 300, 200 or 100
%
% exportSlimObj writes objFL/objGC and kinFL/kinGC with ONE copy of everything
% that does not depend on sample rate (bp, cluid, probeOfUnit, trialid, params,
% featLeg ...) and one copy per rate of what does:
%
%   obj:  time300 time200 time100  psth300 ...  trialdat300 ...
%   kin:  time300 time200 time100  dat300 dat200 dat100
%
% This returns them with the chosen rate's arrays under the ordinary names
% (obj.time, obj.psth, obj.trialdat, kin.time, kin.dat) and every other rate's
% arrays removed, and sets obj.dt / obj.sampleRate / obj.params.dt to match --
% so it is a drop-in for what loadSessionData + getKinematics used to return:
%
%   [obj, kin] = loadSlimSession(fullfile('Data','R1'), 'TD26d_2025_08_06', 'goCue', 200);
%   units = find(obj.probeOfUnit == 2);          % select units by ORIGINAL probe
%   X     = obj.trialdat(:, units, :);           % time x unit x trial
%   tl    = kin.dat(:, :, strcmp(kin.featLeg, 'tongue_length'));
%
%   obj = loadSlimSession(...)      % one output: the kin file is not read at all
%
% Only the requested alignment is read from disk. The other rates of that
% alignment are read too (MATLAB cannot load part of a struct variable) and
% then dropped.
%
% Checks: obj and kin belong to the same session and alignment, the rate
% exists, obj and kin agree on the time base, and kin has no more trials than
% trialdat.

    switch lower(alignEvent)
        case {'firstlick','fl'}, tag = 'FL';  ev = 'firstLick';
        case {'gocue','gc'},     tag = 'GC';  ev = 'goCue';
        otherwise
            error('alignEvent must be ''firstLick'' or ''goCue'', not ''%s''.', alignEvent);
    end

    objF = fullfile(dirPath, [stem '_obj.mat']);
    kinF = fullfile(dirPath, [stem '_kin.mat']);
    assert(exist(objF, 'file') == 2, 'No such file: %s', objF);
    assert(exist(kinF, 'file') == 2, 'No such file: %s', kinF);

    obj = readOne(objF, 'obj', tag, rate, {'time','psth','trialdat'});
    if nargout < 2, return; end                 % obj only: the kin file is not read
    kin = readOne(kinF, 'kin', tag, rate, {'time','dat'});

    % ---- pairing checks ------------------------------------------------------
    assert(strcmp(obj.task, kin.task) && strcmp(obj.anm, kin.anm) && strcmp(obj.date, kin.date), ...
        'obj is %s %s %s but kin is %s %s %s.', obj.task, obj.anm, obj.date, kin.task, kin.anm, kin.date);
    assert(strcmp(obj.alignEvent, ev) && strcmp(kin.alignEvent, ev), ...
        'Asked for %s; obj is %s-aligned and kin is %s-aligned.', ev, obj.alignEvent, kin.alignEvent);
    assert(numel(obj.time) == numel(kin.time) && max(abs(obj.time - kin.time)) < 1e-9, ...
        '%s %s: obj and kin disagree about the %s time base at %d Hz.', obj.anm, obj.date, ev, rate);
    % getKinematicsFromVideo uses FEWER trials than bp.Ntrials on some sessions on
    % purpose (its per-session nTrials overrides, e.g. TD8d 2024-09-07), so kin
    % may be shorter than trialdat -- never longer.
    if isfield(obj, 'trialdat')
        assert(size(kin.dat, 2) <= size(obj.trialdat, 3), ...
            '%s %s: kin has %d trials, trialdat only %d.', obj.anm, obj.date, ...
            size(kin.dat, 2), size(obj.trialdat, 3));
    end
end


function x = readOne(fname, pre, tag, rate, names)
% Two file layouts exist in Data:
%   current  objFL / objGC, with time<r> psth<r> trialdat<r> inside  (kin: dat<r>)
%   earlier  objFL300 objFL200 ... one flat struct per rate           (6 R1 sessions:
%            TD10si 2024-07-09/10/11/13/14 and TD9si 2024-07-05 were exported with
%            the earlier version and kept when the export was re-run)
% Both hold the same data; this reads whichever the file has.
    w = whos('-file', fname);  have = {w.name};
    vNew = [pre tag];
    vOld = sprintf('%s%s%d', pre, tag, rate);
    if any(strcmp(have, vNew))
        L = load(fname, vNew);
        x = pickRate(L.(vNew), rate, names, fname);
    elseif any(strcmp(have, vOld))
        L = load(fname, vOld);
        x = L.(vOld);                                     % already one rate, flat
        assert(isfield(x, 'sampleRate') && x.sampleRate == rate, ...
            '%s: %s does not hold %d Hz data.', fname, vOld, rate);
        x.rates = rate;
        if isfield(x, 'params'), x.params.dt = 1/rate; end
    else
        error('%s has neither %s nor %s. Variables: %s', fname, vNew, vOld, strjoin(have, ', '));
    end
end


function x = pickRate(x, rate, names, fname)
% Put the chosen rate's arrays under the plain names and drop every rate copy.
    assert(isfield(x, 'rates'), ['%s has no .rates field -- it was written by an ' ...
        'older version of exportSlimObj.'], fname);
    assert(any(x.rates == rate), '%s: no %d Hz data. Available: %s Hz.', ...
        fname, rate, mat2str(x.rates));

    for n = 1:numel(names)
        f = sprintf('%s%d', names{n}, rate);
        if isfield(x, f), x.(names{n}) = x.(f); end
    end
    for r = x.rates
        for n = 1:numel(names)
            f = sprintf('%s%d', names{n}, r);
            if isfield(x, f), x = rmfield(x, f); end
        end
    end

    x.dt         = 1/rate;
    x.sampleRate = rate;
    if isfield(x, 'params'), x.params.dt = 1/rate; end
end

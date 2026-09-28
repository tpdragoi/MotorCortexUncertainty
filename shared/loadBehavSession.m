function [obj, kin, params] = loadBehavSession(dataDir, anmHint, dateStr, params)
% LOADBEHAVSESSION  Drop-in for the behavior-only pipeline chain
%
%     [obj, params] = loadSessionData(meta, params, params.behav_only);
%     me  = loadMotionEnergy(obj, meta, params, datapth);
%     kin = getKinematics(obj, me, params);
%
% reading ONLY the slim session files in dataDir (<ANM>_<YYYY_MM_DD>_obj.mat /
% _kin.mat).
%
% What it reproduces, and how:
%   alignment   params.alignEvent picks objGC/kinGC ('goCue') or objFL/kinFL ('firstLick')
%   dt          params.dt picks the stored rate (1/300, 1/200, 1/100). Any other
%               dt is an error -- nothing is resampled.
%   time        loadSessionData sets obj.time to the bin CENTERS of
%               params.tmin:params.dt:params.tmax. The files hold a longer
%               window on the same grid, so obj.time, kin.dat (and trialdat/psth)
%               are cropped to exactly those samples; a window outside the stored
%               one, or a grid that does not line up, is an error.
%   trialid     recomputed with the pipeline's own findTrials(obj, params.condition),
%               so it follows THIS script's conditions, not the stored ones.
%   bp, pth     as the pipeline returns them (full bp).
%   sglx        obj.sglx.laserTrigIX from <stem>_laser.mat when present (opto sessions).
%   kin         kin.dat (time x trials x feature) and kin.featLeg as getKinematics
%               returns them, on the cropped time base.
% With params.behav_only = 1 the neural arrays are dropped.

    rate = round(1 / params.dt);
    assert(abs(rate * params.dt - 1) < 1e-9, 'params.dt = %g is not 1/integer.', params.dt);

    stem = findSlimStem(dataDir, anmHint, dateStr);
    [obj, kin] = loadSlimSession(dataDir, stem, params.alignEvent, rate);

    % ---- the pipeline's time base, and the matching stored samples ----------
    edges = params.tmin:params.dt:params.tmax;
    tWant = edges + params.dt/2;
    tWant = tWant(1:end-1);                       % exactly loadSessionData
    tHave = obj.time(:)';
    idx   = round((tWant - tHave(1)) / params.dt) + 1;
    if any(idx < 1) || any(idx > numel(tHave))
        error(['%s: window [%g %g] s is outside the export (%.4f to %.4f s). ' ...
               'Narrow params.tmin/tmax.'], stem, params.tmin, params.tmax, ...
               tHave(1) - params.dt/2, tHave(end) + params.dt/2);
    end
    assert(max(abs(tHave(idx) - tWant)) < params.dt * 1e-3, ...
        '%s: the exported %d Hz grid does not line up with tmin = %g.', stem, rate, params.tmin);

    obj.time = tWant;                             % row, as loadSessionData
    kin.time = tWant;
    kin.dat  = kin.dat(idx, :, :);
    kin.dat  = double(kin.dat);
    if isfield(params, 'behav_only') && params.behav_only
        for f = {'psth', 'trialdat'}
            if isfield(obj, f{1}), obj = rmfield(obj, f{1}); end
        end
    else
        if isfield(obj, 'psth'),     obj.psth     = double(obj.psth(idx, :, :));     end
        if isfield(obj, 'trialdat'), obj.trialdat = double(obj.trialdat(idx, :, :)); end
    end

    % same nesting as the pipeline's params.trialid (obj.trialid is stored in
    % the shape loadSessionData returns, so copy that shape)
    tid = findTrials(obj, params.condition);
    if isfield(obj, 'trialid') && iscell(obj.trialid) && numel(obj.trialid) == 1 && iscell(obj.trialid{1})
        tid = {tid};
    end
    params.trialid = tid;
    obj.trialid    = tid;

    % laser triggers for the opto sessions:
    % obj.sglx.laserTrigIX exactly as the pipeline's obj carries it
    lf = fullfile(dataDir, [stem '_laser.mat']);
    if exist(lf, 'file') == 2
        L = load(lf, 'laserTrigIX', 'Ntrials');
        assert(L.Ntrials == obj.bp.Ntrials, '%s: laser file has %d trials, bp has %d.', stem, L.Ntrials, obj.bp.Ntrials);
        obj.sglx.laserTrigIX = L.laserTrigIX;
    end
    params.slimStem = stem;
end

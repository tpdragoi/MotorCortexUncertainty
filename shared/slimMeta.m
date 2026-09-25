function meta = slimMeta(loader, dateStr, varargin)
% SLIMMETA  Stand-in for the per-animal loader functions (loadTD10s_neur,
% loadTD13_neural, loadTD_inactivation_TD4f, ...). Instead of pointing at the
% raw data tree it points at the exported slim files in Data/.
%
%   meta = slimMeta(@loadTD26_neur, '2025-08-06')
%   meta = slimMeta('loadTD26_neur', '2025-08-06')
%   meta = slimMeta('TD26d', '2025-08-06')           % animal name also works
%
% The analysis scripts used to do
%       meta5 = loadTD26_neur(meta5, datapth, date);
% and now do
%       meta5 = slimMeta('loadTD26_neur', date);
% Any extra arguments are ignored, so a call in the old argument order can be
% converted by swapping only the function name.
%
% meta fields:  anm, date, loader, task, stem, objFile, kinFile, dataRoot, probe
% meta.probe is left empty; slimToLegacy fills params.probe from the file.
%
% Data folder: <folder containing this file>/Data, or <its parent>/Data (MATLAB Codes _ v2\shared ->
% MATLAB Codes _ v2\Data), unless the global variable
% SLIM_DATA_ROOT is set to another folder.

    if isa(loader, 'function_handle'), name = func2str(loader); else, name = char(loader); end
    name = regexprep(name, '^@', '');
    dateStr = char(dateStr);

    anm      = animalFromLoader(name);
    dataRoot = slimDataRoot();

    stemDate = strrep(dateStr, '-', '_');
    % '*' after the date: the one session whose recorded date carries a typo
    % (TDl3 2025-02-04 is stored as ..._2025_02_044_obj.mat) is still found.
    hits = dir(fullfile(dataRoot, '*', sprintf('%s_%s*_obj.mat', anm, stemDate)));
    % a longer date that merely starts with this one is not a match unless it is
    % the known typo form (one extra digit)
    keep = false(numel(hits), 1);
    for i = 1:numel(hits)
        rest = hits(i).name(numel(anm) + 2 + numel(stemDate) : end);   % after the date
        keep(i) = strcmp(rest, '_obj.mat') || ~isempty(regexp(rest, '^\d_obj\.mat$', 'once'));
    end
    hits = hits(keep);

    assert(~isempty(hits), ['slimMeta: no exported file for %s %s (loader %s) under %s.\n' ...
        'Expected <task>/%s_%s_obj.mat. Was this session exported?'], ...
        anm, dateStr, name, dataRoot, anm, stemDate);
    assert(numel(hits) == 1, 'slimMeta: %s %s matches %d files: %s', anm, dateStr, ...
        numel(hits), strjoin(fullfile({hits.folder}, {hits.name}), ', '));

    [~, task] = fileparts(hits.folder);
    meta = struct();
    meta.anm      = anm;
    meta.date     = dateStr;
    meta.loader   = name;
    meta.task     = task;
    meta.stem     = hits.name(1:end-8);                        % strip '_obj.mat'
    meta.objFile  = fullfile(hits.folder, hits.name);
    meta.kinFile  = fullfile(hits.folder, [meta.stem '_kin.mat']);
    meta.dataRoot = dataRoot;
    meta.probe    = [];
    assert(exist(meta.kinFile, 'file') == 2, 'slimMeta: %s has no matching kin file %s.', ...
        meta.objFile, meta.kinFile);
end


function anm = animalFromLoader(name)
% Loader function name (or animal name) -> animal name as stored in the files.
    % loadTD_inactivation_TD4f -> TD4f
    tok = regexp(name, '^loadTD_inactivation_(\w+)$', 'tokens', 'once');
    if ~isempty(tok), anm = tok{1}; return; end

    code = regexprep(name, '^load', '');
    code = regexprep(code, '_.*$', '');          % drop _neur, _neural, _neur222, _many ...
    table = { ...
        'TD10s','TD10si';  'TD9s','TD9si';   'TD27','TD27d';  'TD26','TD26d'; ...
        'TD1','TD1d';      'TD4','TD4d';     'TD13','TD13d';  'TD15','TD15d'; ...
        'TD8','TD8d';      'TD22','TD22d';   'TD23','TD23d';  'TD24','TD24d'; ...
        'TD25','TD25d';    'TD4f','TD4f';    'TD5f','TD5f';   'TD6f','TD6f'; ...
        'TD7f','TD7f';     'YH1','YH1';      'YH2','YH2';     'TDv1','TDv1'; ...
        'TDv4','TDv4';     'TDv5','TDv5';    'TDv6','TDv6'; ...
        'TD2l','TDl2';     'TD3l','TDl3';    'TD4l','TDl4';   'TD5l','TD5l'};
    r = find(strcmp(code, table(:,1)), 1);
    if ~isempty(r), anm = table{r,2}; return; end
    r = find(strcmp(name, table(:,2)), 1);          % already an animal name
    if ~isempty(r), anm = table{r,2}; return; end
    error('slimMeta: do not know which animal loader "%s" loads. Add it to the table in slimMeta.m.', name);
end


function d = slimDataRoot()
    global SLIM_DATA_ROOT %#ok<GVMIS>
    if ~isempty(SLIM_DATA_ROOT)
        d = SLIM_DATA_ROOT;
    else
        here = fileparts(mfilename('fullpath'));
        % [data] the obj/kin files live in Data, either beside this folder
        % or one level up from it.
        cand = { fullfile(here, 'Data'), fullfile(fileparts(here), 'Data') };
        d = cand{end};
        for ci = 1:numel(cand)
            if exist(cand{ci}, 'dir') == 7, d = cand{ci}; break; end
        end
    end
    assert(exist(d, 'dir') == 7, 'slimMeta: data folder %s does not exist.', d);
end

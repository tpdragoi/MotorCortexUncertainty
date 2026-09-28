function P = setPaths()
% SETPATHS  Tells the code where the data are. Edit the lines under
% "EDIT HERE" once, after downloading the data; nothing else needs changing.
%
% Every figure script calls this at the top. It returns a struct with the
% folders below and puts shared/ and shared/pipelineCopies on the MATLAB path.
%
%   dataRoot    the unzipped data release: the folder that holds the task
%               folders R1, R14, R16, VTA, Learning, C4Stim_R1, C4Stim_R16
%               and GCStim, each full of <ANM>_<DATE>_obj.mat / _kin.mat files
%   hmmRoot     the unzipped "Disengagement Times" folder (HMM-GLM fits, one
%               subfolder per session and probe). Only F2L, F2M, FS6B and
%               FS10_disengageTimeR14 read it.
%   cacheRoot   where the decoding and spike-rate scripts keep their caches
%               and the per-task summary files the cross-task comparisons
%               read. Can grow to several GB; safe to delete at any time.
%   outputRoot  where the scripts in utilities/ write their CSV tables.
%
% Leave a value empty ('') to use the default, a folder inside this
% repository: Data, Disengagement Times, cache and output. So there are two
% ways to set up: point dataRoot and hmmRoot at wherever you unzipped the
% data, or move the unzipped Data and Disengagement Times folders into the
% repository folder. All four default folders are listed in .gitignore.
%
% Paths may use / or \ on any system, and may start with ~ on macOS/Linux.

%% ============================ EDIT HERE ============================
dataRoot   = '';   % e.g. 'D:\CorticalDisengagement\Data'   or '~/data/CorticalDisengagement/Data'
hmmRoot    = '';   % e.g. 'D:\CorticalDisengagement\Disengagement Times'
cacheRoot  = '';   % optional; default <repository>/cache
outputRoot = '';   % optional; default <repository>/output
%% ===================================================================

repoRoot = fileparts(mfilename('fullpath'));

P.repoRoot   = repoRoot;
P.dataRoot   = resolveDir(dataRoot,   fullfile(repoRoot, 'Data'));
P.hmmRoot    = resolveDir(hmmRoot,    fullfile(repoRoot, 'Disengagement Times'));
P.cacheRoot  = resolveDir(cacheRoot,  fullfile(repoRoot, 'cache'));
P.outputRoot = resolveDir(outputRoot, fullfile(repoRoot, 'output'));

% shared functions (added once; addpath is skipped when already on the path)
onPath = [pathsep path pathsep];
for d = {fullfile(repoRoot, 'shared'), fullfile(repoRoot, 'shared', 'pipelineCopies')}
    if ~contains(onPath, [pathsep d{1} pathsep], 'IgnoreCase', ispc)
        addpath(d{1});
    end
end

if exist(P.dataRoot, 'dir') ~= 7
    error('setPaths:noData', ['Data folder not found:\n  %s\n' ...
        'Download the data release, unzip it, and set dataRoot at the top of\n  %s\n' ...
        'to the folder that contains R1, R14, R16, VTA, Learning, C4Stim_R1, C4Stim_R16 and GCStim.'], ...
        P.dataRoot, fullfile(repoRoot, 'setPaths.m'));
end
if exist(P.cacheRoot, 'dir') ~= 7
    [ok, msg] = mkdir(P.cacheRoot);
    assert(ok, 'setPaths: cannot create cache folder %s (%s). Set cacheRoot in setPaths.m.', P.cacheRoot, msg);
end
end


function d = resolveDir(d, default)
% Empty -> default; ~ -> home folder; separators -> this system's; no trailing separator.
    d = strtrim(char(d));
    if isempty(d), d = default; end
    if startsWith(d, '~') && ~ispc
        d = [getenv('HOME') d(2:end)];
    end
    d = strrep(strrep(d, '\', filesep), '/', filesep);
    while numel(d) > 1 && d(end) == filesep && ~(ispc && numel(d) == 3 && d(2) == ':')
        d(end) = [];
    end
end

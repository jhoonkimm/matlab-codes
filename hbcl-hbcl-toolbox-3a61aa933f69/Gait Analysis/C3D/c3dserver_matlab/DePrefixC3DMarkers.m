% GetRidOfC3DMarkerPrefix.m

itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package

% get all the files to remove prefixes from
filenames = {};
pathnames = {};
for i = 1:100
    fnms = {}; pnms = {};
    [fnms, pnms{1}] = uigetfile('*.c3d','Pick C3D files to Deprefix','MultiSelect','on');
    if ~iscell(fnms); if ~fnms; break; end; fnms = {fnms}; end
    pnms = repmat(pnms,length(fnms),1);
    filenames = cat(1,filenames,fnms{:});
    pathnames = cat(1,pathnames,pnms{:});
end
if ~filenames{end}; error('File Selection Cancelled'); end

% rework file names so the files go into the "Modified at IPS" folders

%% now open files one by one

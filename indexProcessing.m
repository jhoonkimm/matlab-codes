function output = indexProcessing(path,folder,varargin)
%  
% output = indexProcessing(path,folder,varargin)
% 
%  Outputs .mat, .txt and .c3d files after organizing.
%
%  Flags:
%       integrated
%           Find c3d files that are already integrated. Default is false
%
%
%  Created by Jay Kim (12/15/2017)

%% Default input arguments
integrated = 0;

%% Optional input arguments

for i = 1:2:length(varargin)
  switch varargin{i}
    case 'integrated'
    integrated = varargin{i + 1};
  end
end

%% Processing directory files
output = struct;

dirSBD = dir(fullfile(path,folder)); %get directory content from folder specified
dirSBDnames = {dirSBD.name}; % make directory names into cell array

%% If integrated files wanted
if integrated == 1
    idxDirSBD = find(not(cellfun('isempty',regexp(dirSBDnames,'.c3d')))); % find all c3d files in cell array
    idxCleaned = [];
    for i = idxDirSBD
        if regexp(lower(dirSBDnames{i}),'integrated') % INCLUDE all "integrated" c3d files
            idxCleaned =    [idxCleaned,   i];
        end
        if regexp(lower(dirSBDnames{i}),'static')
            idxStatic = i;
        end
    end
    idxPWSbase = []; idx150base = []; idx75base = []; idxSplit = []; idxPost = [];
    for i = idxCleaned
        if regexp(lower(dirSBDnames{i}),'pws')
        idxPWSbase =    [idxPWSbase,   i];
        elseif regexp(lower(dirSBDnames{i}),'150')
        idx150base =    [idx150base,   i];
        elseif regexp(lower(dirSBDnames{i}),'75')
        idx75base =     [idx75base, i];
        elseif regexp(lower(dirSBDnames{i}),'split')
        idxSplit =      [idxSplit,   i];
        elseif regexp(lower(dirSBDnames{i}),'post')
        idxPost =       [idxPost, i];
        end
    end
    if length(idxSplit)>1
        error('More than one split file detected')
    end
    if length(idxPost)>1
        error('More than one post file detected')
    end
    output.c3d.pws      = idxPWSbase;
    output.c3d.fast     = idx150base;
    output.c3d.slow     = idx75base;
    output.c3d.split    = idxSplit;
    output.c3d.post     = idxPost;
    
    output.static       = idxStatic;

end
%% No integrated files
if integrated==0
    idxDirSBD = find(not(cellfun('isempty',regexp(dirSBDnames,'.c3d')))); % find all c3d files in cell array
    idxDirLog = find(not(cellfun('isempty',regexp(dirSBDnames,'.txt')))); % find all txt files in cell array
    idxDirMat = find(not(cellfun('isempty',regexp(dirSBDnames,'.mat')))); % find all mat files in cell array

    idxCleaned = [];
    for i = idxDirSBD
        if regexp(lower(dirSBDnames{i}),'cleaned') % find all "cleaned" c3d files
            if isempty(regexp(lower(dirSBDnames{i}),'integrated')) %#ok<RGXP1> % remove all "integrated" c3d files
                if isempty(regexp(lower(dirSBDnames{i}),'.mat')) %#ok<RGXP1> % remove all "integrated" mat files
                idxCleaned =    [idxCleaned,   i];
                end
            end
        end
        if regexp(lower(dirSBDnames{i}),'static')
            idxStatic = i;
        end
    end

    idxPWSbase = []; idx150base = []; idx75base = []; idxSplit = []; idxPost = [];
    for i = idxCleaned
        if regexp(lower(dirSBDnames{i}),'pws')
        idxPWSbase =    [idxPWSbase,   i];
        elseif regexp(lower(dirSBDnames{i}),'150')
        idx150base =    [idx150base,   i];
        elseif regexp(lower(dirSBDnames{i}),'75')
        idx75base =     [idx75base, i];
        elseif regexp(lower(dirSBDnames{i}),'split')
        idxSplit =      [idxSplit,   i];
        elseif regexp(lower(dirSBDnames{i}),'post')
        idxPost =       [idxPost, i];
        end
    end

    idxPWSbaseMat = []; idx150baseMat = []; idx75baseMat = []; idxSplitMat = []; idxPostMat = [];
    for i = idxDirMat
        if isempty(regexp(lower(dirSBDnames{i}),'processed', 'once'))
            if regexp(lower(dirSBDnames{i}),'pws')
            idxPWSbaseMat =    [idxPWSbaseMat,   i];
            elseif regexp(lower(dirSBDnames{i}),'150')
            idx150baseMat =    [idx150baseMat,   i];
            elseif regexp(lower(dirSBDnames{i}),'75')
            idx75baseMat =     [idx75baseMat, i];
            elseif regexp(lower(dirSBDnames{i}),'split')
            idxSplitMat =      [idxSplitMat,   i];
            elseif regexp(lower(dirSBDnames{i}),'post')
            idxPostMat =       [idxPostMat, i];
            end
        end
    end

    idxPWSbaseLog = []; idx150baseLog = []; idx75baseLog = []; idxSplitLog = []; idxPostLog = [];
    for i = idxDirLog
        if regexp(lower(dirSBDnames{i}),'pws')
        idxPWSbaseLog =    [idxPWSbaseLog,   i];
        elseif regexp(lower(dirSBDnames{i}),'150')
        idx150baseLog =    [idx150baseLog,   i];
        elseif regexp(lower(dirSBDnames{i}),'75')
        idx75baseLog =     [idx75baseLog, i];
        elseif regexp(lower(dirSBDnames{i}),'split')
        idxSplitLog =      [idxSplitLog,   i];
        elseif regexp(lower(dirSBDnames{i}),'post')
        idxPostLog =       [idxPostLog, i];
        end
    end

    output.txt.pws      = idxPWSbaseLog;
    output.txt.fast     = idx150baseLog;
    output.txt.slow     = idx75baseLog;
    output.txt.split    = idxSplitLog;
    output.txt.post     = idxPostLog;

    output.mat.pws      = idxPWSbaseMat;
    output.mat.fast     = idx150baseMat;
    output.mat.slow     = idx75baseMat;
    output.mat.split    = idxSplitMat;
    output.mat.post     = idxPostMat;

    output.c3d.pws      = idxPWSbase;
    output.c3d.fast     = idx150base;
    output.c3d.slow     = idx75base;
    output.c3d.split    = idxSplit;
    output.c3d.post     = idxPost;
    
    output.static       = idxStatic;

end

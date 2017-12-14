clear all; close all; clc;

% Treadmill Integration for Jay's study
% November 3rd 2017
%---------------------------------------------------------------------------------------------------------------------------------------%
system = 'WIN';

if strcmp(system,'OS')
    path  = '\\Mac\Home\Google Drive\3-Research\2016\Trials';
    path2 = '\\Mac\Home\Google Drive\3-Research';
    path3 = '\\Mac\Home\Google Drive\3-Research\PaperComplianceStudy\FigureData';
elseif strcmp(system,'WIN')
    path  = 'E:\RESEARCH\MDSC 508';
end
%---------------------------------------------------------------------------------------------------------------------------------------%

folder = 'SBD_003';

customLog = 'phasespaceCustomLog1_11072017-PWS_01.txt';
treadmillData = '11072017-PWS_01.mat';
c3dFile = 'capture-20171107-PWS_01_CLEANED';

%---------------------------------------------------------------------------------------------------------------------------------------%
dirSBD = dir(fullfile(path,folder));
dirSBDnames = {dirSBD.name};
idxDirSBD = find(not(cellfun('isempty',regexp(dirSBDnames,'.c3d'))));
idxDirLog = find(not(cellfun('isempty',regexp(dirSBDnames,'.txt'))));
idxDirMat = find(not(cellfun('isempty',regexp(dirSBDnames,'.mat'))));

idxCleaned = [];
for i = idxDirSBD
    if regexp(lower(dirSBDnames{i}),'cleaned')
    idxCleaned =    [idxCleaned,   i];
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

[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2017(...
      fullfile(path,folder,dirSBDnames{idxPWSbaseLog}), fullfile(path,folder,dirSBDnames{idxPWSbaseMat}), fullfile(path,folder,dirSBDnames{idxPWSbase}), fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idxPWSbase})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2017(...
      fullfile(path,folder,dirSBDnames{idx150baseLog}), fullfile(path,folder,dirSBDnames{idx150baseMat}), fullfile(path,folder,dirSBDnames{idx150base}), fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx150base})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2017(...
      fullfile(path,folder,dirSBDnames{idx75baseLog}), fullfile(path,folder,dirSBDnames{idx75baseMat}), fullfile(path,folder,dirSBDnames{idx75base}), fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx75base})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
  
% for k = idxSplit
%     [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2017(...
%       fullfile(path,folder,dirSBDnames{idxSplitLog}), fullfile(path,folder,dirSBDnames{idxSplitMat}), fullfile(path,folder,dirSBDnames{k}), fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{k})),...
%       'c3dUnitsAreInMeters',1,...
%       'shouldFillGapsInMarkerData',0);
% end



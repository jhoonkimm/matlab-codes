function integrationC3dAndTreadmill(path,folder,varargin)
%  
% output = integrationC3dAndTreadmill(path,folder)
% 
%  Outputs integrated c3d files for the entire folder for MDSC 508 purposes
%
%  Flags:
%   exclude
%       excludes only 'pws','fast','slow','split' or 'post' c3d files. Can be cell array.
%
%
%  Created by Jay Kim (11/09/2017)

%% Default input arguments
exclude = [];


%% Optional input arguments

for i = 1:2:length(varargin)
  switch varargin{i}
    case 'exclude'
    exclude = varargin{i + 1};
  end
end
%% Recording command window log

diary(fullfile(path,folder,sprintf('integrationLog_%s.txt',datestr(now(),'yyyymmdd_HHMM'))))
%Clearing message window for diary
%clc;

%% Processing directories

dirSBD = dir(fullfile(path,folder));
dirSBDnames = {dirSBD.name};

idx = indexProcessing(path,folder);

%Finding missing indices - leading to error message
if mean(cellfun('isempty',struct2cell(idx.txt)))>0
    error('Function cannot run, one or more phasespaceCustomLog indices not found')
elseif mean(cellfun('isempty',struct2cell(idx.mat)))>0
    error('Function cannot run, one or more force mat file indices not found')
elseif mean(cellfun('isempty',struct2cell(idx.mat)))>0
    error('Function cannot run, one or more c3d file indices not found')
end

%% PWS Baseline Processing
if ~strcmp(exclude,'pws')
tstart = tic;
fprintf('Starting to integrate %s PWS Baseline...\n',folder)
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
      fullfile(path,folder,dirSBDnames{idx.txt.pws}), ...
      fullfile(path,folder,dirSBDnames{idx.mat.pws}),...
      fullfile(path,folder,dirSBDnames{idx.c3d.pws}),...
      fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.pws})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
t = toc(tstart);
fprintf('PWS Baseline took %.4f seconds to complete\n',t)
fprintf('Ending %s PWS Baseline integration...\n',folder)
fprintf('\n')
fprintf('\n')
end

%% Fast Baseline Processing
if ~strcmp(exclude,'fast')
tstart = tic;
fprintf('Starting to integrate %s Fast Baseline...\n',folder)
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
      fullfile(path,folder,dirSBDnames{idx.txt.fast}), ...
      fullfile(path,folder,dirSBDnames{idx.mat.fast}),...
      fullfile(path,folder,dirSBDnames{idx.c3d.fast}),...
      fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.fast})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
t = toc(tstart);
fprintf('Fast Baseline took %.4f seconds to complete\n',t)
fprintf('Ending %s Fast Baseline integration...\n',folder)
fprintf('\n')
fprintf('\n')
end

%% Slow Baseline Processing
if ~strcmp(exclude,'slow')
tstart = tic;
fprintf('Starting to integrate %s Slow Baseline...\n',folder)
[ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
      fullfile(path,folder,dirSBDnames{idx.txt.slow}), ...
      fullfile(path,folder,dirSBDnames{idx.mat.slow}),...
      fullfile(path,folder,dirSBDnames{idx.c3d.slow}),...
      fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.slow})),...
      'c3dUnitsAreInMeters',1,...
      'shouldFillGapsInMarkerData',0);
t = toc(tstart);
fprintf('Slow Baseline took %.4f seconds to complete\n',t)
fprintf('Ending %s Slow Baseline integration...\n',folder)
fprintf('\n')
fprintf('\n')
end

%% Split Processing
if ~strcmp(exclude,'split')
tstart = tic;
fprintf('Starting to integrate %s Split ...\n',folder)
%finding merged file name
firstSplit = dirSBDnames{idx.c3d.split(1)};
firstSplit = firstSplit(1,1:end-15); 
mergedFileName = sprintf('%sMERGED.c3d',firstSplit);
idxMerged = strfind(dirSBDnames, mergedFileName);

% Single phasespaceCustomLog for multilpe c3d files
    if length(idx.txt.split)==1
        if find(not(cellfun('isempty',idxMerged)))==[];
            if length(idx.c3d.split)==2
                mergeC3D(fullfile(path,folder,dirSBDnames{idx.c3d.split(1)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(2)}),...
                         fullfile(path,folder,mergedFileName));
            elseif length(idx.c3d.split)==3
                mergeC3D(fullfile(path,folder,dirSBDnames{idx.c3d.split(1)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(2)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(3)}),...
                         fullfile(path,folder,mergedFileName));
            elseif length(idx.c3d.split)==4
                mergeC3D(fullfile(path,folder,dirSBDnames{idx.c3d.split(1)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(2)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(3)}),...
                         fullfile(path,folder,dirSBDnames{idx.c3d.split(4)}),...
                         fullfile(path,folder,mergedFileName));
            end
        end
    end
% Multilpe phasespaceCustomLog for multiple c3d files
    if length(idx.txt.split)>1
        for k = 1:length(idx.txt.split)
    [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
          fullfile(path,folder,dirSBDnames{idx.txt.split(k)}), ...
          fullfile(path,folder,dirSBDnames{idx.mat.split(k)}),...
          fullfile(path,folder,dirSBDnames{idx.c3d.split(k)}),...
          fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(k)})),...
          'c3dUnitsAreInMeters',1,...
          'shouldFillGapsInMarkerData',0);
        end
        if length(idx.txt.split)==2
                mergeC3D(fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(1)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(2)})),...
                         fullfile(path,folder,sprintf('integrated_%s',mergedFileName)));
        elseif length(idx.txt.split)==3
                mergeC3D(fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(1)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(2)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(3)})),...
                         fullfile(path,folder,sprintf('integrated_%s',mergedFileName)));
        elseif length(idx.txt.split)==4
                mergeC3D(fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(1)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(2)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(3)})),...
                         fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.split(4)})),...
                         fullfile(path,folder,sprintf('integrated_%s',mergedFileName)));
        end

    elseif length(idx.txt.split)==1
    [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
          fullfile(path,folder,dirSBDnames{idx.txt.split}), ...
          fullfile(path,folder,dirSBDnames{idx.mat.split}),...
          fullfile(path,folder,mergedFileName),...
          fullfile(path,folder,sprintf('integrated_%s',mergedFileName)),...
          'c3dUnitsAreInMeters',1,...
          'shouldFillGapsInMarkerData',0);
    end

    t = toc(tstart);
    fprintf('Split took %.4f seconds to complete\n',t)
    fprintf('Ending %s Split integration...\n',folder)
    fprintf('\n')
    fprintf('\n')
end

%% Split Processing
if ~strcmp(exclude,'post')
    tstart = tic;
    fprintf('Starting to integrate %s Post ...\n',folder)
    %finding merged file name
    firstPost = dirSBDnames{idx.c3d.post(1)};
    firstPost = firstPost(1,1:end-15); 
    mergedFileName = sprintf('%sMERGED',firstPost);
    if length(idx.c3d.post)==2
        mergeC3D(fullfile(path,folder,dirSBDnames{idx.c3d.post(1)}),...
                 fullfile(path,folder,dirSBDnames{idx.c3d.post(2)}),...
                 fullfile(path,folder,mergedFileName));
    elseif length(idx.c3d.post)==3
        mergeC3D(fullfile(path,folder,dirSBDnames{idx.c3d.post(1)}),...
                 fullfile(path,folder,dirSBDnames{idx.c3d.post(2)}),...
                 fullfile(path,folder,dirSBDnames{idx.c3d.post(3)}),...
                 fullfile(path,folder,mergedFileName));
    end
    if length(idx.c3d.post)>1
    [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
          fullfile(path,folder,dirSBDnames{idx.txt.post}), ...
          fullfile(path,folder,dirSBDnames{idx.mat.post}),...
          fullfile(path,folder,mergedFileName),...
          fullfile(path,folder,sprintf('integrated_%s',mergedFileName)),...
          'c3dUnitsAreInMeters',1,...
          'shouldFillGapsInMarkerData',0);
    end
    if length(idx.c3d.post)==1
        [ phasespaceAndTreadmillStructure ] = savePhasespaceAndTreadmillDataInC3D2018(...
          fullfile(path,folder,dirSBDnames{idx.txt.post}), ...
          fullfile(path,folder,dirSBDnames{idx.mat.post}),...
          fullfile(path,folder,dirSBDnames{idx.c3d.post}),...
          fullfile(path,folder,sprintf('integrated_%s',dirSBDnames{idx.c3d.post})),...
          'c3dUnitsAreInMeters',1,...
          'shouldFillGapsInMarkerData',0);
    end
    t = toc(tstart);
    fprintf('Post took %.4f seconds to complete\n',t)
    fprintf('Ending %s Post integration...\n',folder)
    fprintf('\n')
    fprintf('\n')
end

%% Ending command window log recording

diary off
end


  

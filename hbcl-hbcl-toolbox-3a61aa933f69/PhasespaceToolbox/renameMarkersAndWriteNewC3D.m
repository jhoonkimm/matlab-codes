function renameMarkersAndWriteNewC3D(inputFilename,outputFilename,varargin)

% Modified by Amy Wu from John's code (Last Revised: 08/05/2014)
% Example:
%   inputFilename = 'P1234-02.c3d';
%   outputFilename = 'P1234-02_R.c3d';
%   renamingList = {...
%     {'M013','C7'},...
%     {'M015','CLAV'},...
%     {'M016','LSHO'};
%   renameMarkersWriteC3D(inputFilename,outputFilename,'markerRenaming',renamingList)
  
c3dUnitsAreInMeters = 1;
mocapSkipFactor = 1;
mocapFrequency = 480;
componentsToPlot = 3; %1; %2; %
shouldPlot = 0; %1; %
printFrequency = 50000;
markerRenaming = {};
templateFitErrorMaximum = 0.025; %0.015; %0.005; %0.01; %0.03;
maxPointIndexToLookAt = []; %1000

for i = 1:2:length(varargin)
  if strcmp('c3dUnitsAreInMeters', varargin{i})
    c3dUnitsAreInMeters = varargin{i + 1};
  end
  if strcmp('mocapSkipFactor', varargin{i})
    mocapSkipFactor = varargin{i + 1};
  end
  if strcmp('mocapFrequency', varargin{i})
    mocapFrequency = varargin{i + 1};
  end
  if strcmp('componentsToPlot', varargin{i})
    componentsToPlot = varargin{i + 1};
  end
  if strcmp('shouldPlot', varargin{i})
    shouldPlot = varargin{i + 1};
  end
  if strcmp('printFrequency', varargin{i})
    printFrequency = varargin{i + 1};
  end
  if strcmp('markerRenaming', varargin{i})
    markerRenaming = varargin{i + 1};
  end
  if strcmp('templateFitErrorMaximum', varargin{i})
    templateFitErrorMaximum = varargin{i + 1};
  end
  if strcmp('maxPointIndexToLookAt', varargin{i})
    maxPointIndexToLookAt = varargin{i + 1};
  end
end

%%
forceReprocessC3Ds = 0;
rawData = loadPhasespaceRecord(inputFilename, ...
  'c3dUnitsAreInMeters', c3dUnitsAreInMeters, 'forceReprocess', forceReprocessC3Ds);

[rawData] = renameMarkersInC3DData(rawData, markerRenaming);
allMarkerData = rmfield(rawData,'markerNames');

% ind = strfind(inputFilename,'c3d');
% getExt = inputFilename(ind:end);
% outputFilename = [inputFilename(1:ind-2) 'R.' getExt];
writeC3DFile(outputFilename, allMarkerData);

fprintf('finished renaming markers in file: %s\n', inputFilename);
fprintf('renamed markers in file: %s\n', outputFilename);




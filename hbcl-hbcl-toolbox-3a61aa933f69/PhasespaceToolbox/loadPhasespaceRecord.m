function [phasespaceRecord] = loadPhasespaceRecord(c3dFile, varargin)

addpath('./C3D');

verbose = 0;
forceReprocess = 0; %1; %
c3dUnitsAreInMeters = 0;
markersToRemove = [];

for i=1:length(varargin)
  if (strcmp(varargin{i}, 'verbose'))
    verbose = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'forceReprocess'))
    forceReprocess = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'c3dUnitsAreInMeters'))
    c3dUnitsAreInMeters = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'markersToRemove'))
    markersToRemove = varargin{i+1};
  end
end

%%

% dotIndex = findstr(c3dFile, '.');
% outputFilename = [c3dFile(1:(dotIndex(end) - 1)) '.mat'];
outputFilename = createMatOutputFilename(c3dFile);

if (forceReprocess || ~exist(outputFilename, 'file'))
  % c3dFile = 'proofTest.c3d';
  itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package
  openc3d(itf, 0, c3dFile);    % applies the correct open options to the file you choose
  
  try
    analogs = getanalogchannels(itf);
    phasespaceRecord.analogs = analogs;
  catch e
  end
  
  %%
  
  % frames = nframes(itf);  % calculates and displays the number of video frames
  frameRateIndex = itf.GetParameterIndex('POINT','RATE'); %% get Index for the Frame Rate
  frameRate = double(itf.GetParameterValue(frameRateIndex, 0));  % get value of the Frame Rate
  
  if (verbose > 0)
    fprintf('frame rate = %g Hz, loading marker positions...\n', frameRate);
  end
  allMocapMarkers = get3dtargets(itf);
  if (verbose > 0)
    fprintf('loaded marker positions\n');
  end
  closec3d(itf);
  
  %%
  allFields = fields(allMocapMarkers);
  %   markerPositions = [];
  phasespaceRecord.markerNames = {};
  for fieldNumber = 1:length(allFields)
    fieldName = allFields(fieldNumber);
    if (strcmp(fieldName, 'units'))
      continue;
    end
    
    thisField = getfield(allMocapMarkers, fieldName{:});
    if (~sum(sum(~isnan(thisField))))
      if (verbose > 0)
        fprintf('skipping %s ', fieldName{:});
      end
      continue
    end
    
    %     error('debug!');
    
    scaleFactor = 0.001;
    if (c3dUnitsAreInMeters)
      scaleFactor = 1.0;
    end
    
    phasespaceRecord.(fieldName{:}) = thisField * scaleFactor; % convert to meters
    phasespaceRecord.markerNames{length(phasespaceRecord.markerNames) + 1} = fieldName{:};
    
    if (verbose > 1)
      %     fprintf('found marker field %s\n', fieldName);
      plot3(thisField(:, 1), thisField(:, 2), thisField(:, 3))
      hold on;
      axis equal;
    end
  end
  
  %   phasespaceRecord.markerPositions = markerPositions;
  phasespaceRecord.frameRate = frameRate;
  phasespaceRecord.time = (0:(length(thisField)-1)) * (1/frameRate);
  save(outputFilename, 'phasespaceRecord');
else
  load(outputFilename);
end

if (~isempty(markersToRemove))
  markersToRemove = eval(markersToRemove);
  indecesToRemove = [];
  %   currentFieldnames = fields(phasespaceRecord);
  currentFieldnames = phasespaceRecord.markerNames;
  for i = 1:length(markersToRemove)
    for j = 1:length(currentFieldnames)
      if (strcmp(markersToRemove{i}, currentFieldnames{j}))
        phasespaceRecord = rmfield(phasespaceRecord, markersToRemove{i});
        indecesToRemove = [indecesToRemove j];
      end
    end
  end
  phasespaceRecord.markerNames(indecesToRemove) = [];
end

if (verbose > 1)
  visualizePhasespaceRecord(phasespaceRecord);
  %   figure;
  %   for fieldName = phasespaceRecord.markerNames
  %     thisField = phasespaceRecord.(fieldName{:});
  %     plot3(thisField(:, 1), thisField(:, 2), thisField(:, 3))
  %     hold on;
  %     axis equal;
  %   end
end



end



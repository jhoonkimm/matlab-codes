function [phasespaceScaleAndFrameAlignment] = ...
  getCartPhasespaceFrameAlignment(phasespaceRecord, cartMotion, varargin)
%% getCartPhasespaceFrameAlignment
% phasespaceRecord, cartMotion can be filenames or loaded data. Filenames
% are preferrable since it then saves a mat file with the calibration.

forceReprocess = 0;
forceReprocessSourceFiles = 0;
verbose = 0;
phasespaceCalibrationIndeces = [];
markerNameToUseInCalibration = [];
for i = 1:2:length(varargin)
  if (strcmp('forceReprocess', varargin{i}))
    forceReprocess = varargin{i + 1};
  end
  if (strcmp('forceReprocessSourceFiles', varargin{i}))
    forceReprocessSourceFiles = varargin{i + 1};
  end
  if (strcmp('verbose', varargin{i}))
    verbose = varargin{i + 1};
  end
  if (strcmp('phasespaceCalibrationIndeces', varargin{i}))
    phasespaceCalibrationIndeces = varargin{i + 1};
  end
  if (strcmp('markerNameToUseInCalibration', varargin{i}))
    markerNameToUseInCalibration = varargin{i + 1};
  end
end

%% check if we inputted data records or files with the data in them
if (ischar(phasespaceRecord))
  if (~(ischar(cartMotion)))
    error('phasespaceRecord and cartMotion must both be strings representing files, or must both be loaded data records');
  end
  cartSlashLocations = strfind(cartMotion, '\');
  cartFilenameRoot = cartMotion(cartSlashLocations(end)+1:end);
  outputFilename = [phasespaceRecord '_' cartFilenameRoot '-alignment.mat'];
  
  phasespaceRecord = loadPhasespaceRecord(phasespaceRecord, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles, ...
    'c3dUnitsAreInMeters', 0);
  alignmentMobileCartRecord = loadMobileCartRecord(cartMotion, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles);
  [cartMotion] = integrateMobileCartMotion(alignmentMobileCartRecord);
else
  if ((ischar(cartMotion)))
    error('phasespaceRecord and cartMotion must both be strings representing files, or must both be loaded data records');
  end
  outputFilename = [];
end

%%
numMarkerTicks = length(phasespaceRecord.(phasespaceRecord.markerNames{1}));
if (isempty(phasespaceCalibrationIndeces))
  phasespaceCalibrationIndeces = 1:(min(numMarkerTicks, 480 * 600));
  %   phasespaceCalibrationIndeces = 1:numMarkerTicks;
  %   phasespaceCalibrationIndeces = 1:(min(numMarkerTicks, 480 * 6));
  %   phasespaceCalibrationIndeces = 1:(min(numMarkerTicks, 480 * 12));
  
  %     phasespaceCalibrationIndeces = (480 * 30):(min(numMarkerTicks, 480 * 50));
  %     phasespaceCalibrationIndeces = (480 * 10):(min(numMarkerTicks, 480 * 40));
  %   phasespaceCalibrationIndeces = (480 * 5):(min(numMarkerTicks, 480 * 20));
  
  
end

%%
if (forceReprocess || isempty(outputFilename) || ~exist(outputFilename, 'file'))
  phasespaceCartLocationCalibration = cartMotion.positions;
  anglesCalibration = cartMotion.angles;
  
  phasespaceTimes = phasespaceRecord.time;
  %   linspace(0, (numMarkerTicks - 1)/phasespaceRecord.frameRate, numMarkerTicks);
  cartTimes = cartMotion.times;
  
  phasespaceCalibrationTimes = phasespaceTimes(phasespaceCalibrationIndeces);
  badCalibrationIndeces = phasespaceCalibrationTimes > cartTimes(end) | ...
    phasespaceCalibrationTimes < cartTimes(1);
  phasespaceCalibrationIndeces(badCalibrationIndeces) = [];
  %   phasespaceCalibrationIndeces = find(phasespaceCalibrationTimes < cartTimes(end));
  phasespaceCalibrationTimes = phasespaceTimes(phasespaceCalibrationIndeces);
  cartCalibrationIndeces = find(cartTimes > phasespaceCalibrationTimes(1), 1, 'first') : ...
    find(cartTimes < phasespaceCalibrationTimes(end), 1, 'last');
  
  % calibrationIndecesOffset = floor(calibrationIndecesOffset);
  % cartCalibrationIndecesToUse = cartCalibrationIndeces + calibrationIndecesOffset;
  cartCalibrationIndecesToUse = cartCalibrationIndeces;
  cartCalibrationIndecesToUse = cartCalibrationIndecesToUse(cartCalibrationIndecesToUse < length(anglesCalibration));
  cartCalibrationIndecesToUse = cartCalibrationIndecesToUse(cartCalibrationIndecesToUse > 0);
  
  if (isempty(markerNameToUseInCalibration))
    allMarkerData = zeros(numMarkerTicks, length(phasespaceRecord.markerNames));
    for i = 1:length(phasespaceRecord.markerNames)
      markerName = phasespaceRecord.markerNames{i};
      allMarkerData(:,i) = phasespaceRecord.(markerName)(:,3);
    end
    [tmp, highestMarkerHeightIndex] = max(nanmean(allMarkerData));
    markerNameToUseInCalibration = phasespaceRecord.markerNames{highestMarkerHeightIndex};
  end
  %%
  
  cartMotionArray = [phasespaceCartLocationCalibration(cartCalibrationIndecesToUse, :) ...
    anglesCalibration(cartCalibrationIndecesToUse)];
  markerDataToUse = phasespaceRecord.(markerNameToUseInCalibration);
  
  phasespaceCalibrationIndeces = phasespaceCalibrationIndeces(phasespaceCalibrationIndeces < length(markerDataToUse));
  markerDataToUse = markerDataToUse(phasespaceCalibrationIndeces, :);
  phasespaceTimesToUse = phasespaceTimes(phasespaceCalibrationIndeces);
  
  calibrationTimeOffset = 0.0;
  cartTimesToUse = cartTimes(cartCalibrationIndecesToUse) + calibrationTimeOffset;
  
  if (verbose > 1)
    %%
    figure;
    subplot(2,1,1);
    plot(cartTimesToUse, cartMotionArray);
    title('calibration cart data');
    subplot(2,1,2);
    plot(phasespaceTimesToUse, markerDataToUse);
    title('calibration marker data');
    %%
  end
  
  %%
  % figure;
  [transformationFromPhasespaceToCart, phasespaceScale, calibrationError] = ..., angleScale] = ...
    findAlignmentBetweenCartAndPhasespaceFrames( ...
    cartTimesToUse, cartMotionArray, ...
    phasespaceTimesToUse, markerDataToUse, ...
    'verbose', verbose, 'skipFrameRegistration', 0);
  %%
  
  phasespaceScaleAndFrameAlignment = [phasespaceScale; transformationFromPhasespaceToCart];
  
%   if (calibrationError > 0.1)
%     error('calibration error was greater than 0.1m, this is probably a bad thing, should be on the order of centimeters');
%   end
  if (calibrationError > 0.04)
%     error('calibration error was greater than 0.04m, this is probably a bad thing, should be on the order of small');
    fprintf('calibration error was greater than 0.04m, this is probably a bad thing, should be on the order of small\n');
    pause;
  end
  
  if (~isempty(outputFilename))
    save(outputFilename, 'phasespaceScaleAndFrameAlignment');
  end
  
else
  load(outputFilename, 'phasespaceScaleAndFrameAlignment');
end


end


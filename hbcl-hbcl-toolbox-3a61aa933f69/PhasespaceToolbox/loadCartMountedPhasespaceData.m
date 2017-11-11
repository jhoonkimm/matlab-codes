
function [markersInWorldFrame, calibration] = loadCartMountedPhasespaceData(...
  trialPhasespaceC3DFilename, trialCartRSTFilename, ...
  calibrationPhasespaceC3DFilename, calibrationCartRSTFilename, ...
  varargin)


phasespaceCalibrationIndeces = [];
forceReprocess = 0;
forceReprocessCalibration = 0;
calibrationToUse = [];
verbose = 0;

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch (option)
    case 'calibrationToUse'
      calibrationToUse = value;
    case 'verbose'
      verbose = value;
    case 'forceReprocessCalibration'
      forceReprocessCalibration = value;
    case 'forceReprocess'
      forceReprocessSourceFiles = value;
    case 'phasespaceCalibrationIndeces'
      phasespaceCalibrationIndeces = value;
  end
end


%%
forceReprocessSourceFiles = forceReprocess;
forceReprocessCalibration = forceReprocessCalibration || forceReprocess;

%%
if (isempty(calibrationToUse))
  
  calibration = getCartPhasespaceFrameAlignment( ...
    calibrationPhasespaceC3DFilename, ...
    calibrationCartRSTFilename, ...
    'forceReprocess', forceReprocess, ...
    'forceReprocessSourceFiles', forceReprocessSourceFiles, ...
    'verbose', verbose, ...
    'phasespaceCalibrationIndeces', phasespaceCalibrationIndeces);
else
  calibration = calibrationToUse;
end

if (verbose > 1)
  fprintf('using calibration [%s]\n', sprintf('%g, ', calibration));
end

%% load trial data
phasespaceData = loadPhasespaceRecord(trialPhasespaceC3DFilename, 'verbose', verbose, 'forceReprocess', forceReprocess);
markerFieldNames = phasespaceData.markerNames; %getMarkerFieldNames(processedPhasespaceTrial);

cartSensors = loadMobileCartRecord(trialCartRSTFilename, 'verbose', verbose, 'forceReprocess', forceReprocess);

%%
if (verbose > 2)
  figure
  colors = jet(length(markerFieldNames));
  for i = 1:length(markerFieldNames)
    markerData = phasespaceData.(markerFieldNames{i});
    plot3(markerData(:,1), markerData(:,2), markerData(:,3), '*', 'Color', colors(i,:));
    hold on;
    axis equal
  end
  title('all markers raw motion');
end

%%
if (~isfield(cartSensors, 'times'))
  cartMotion = integrateMobileCartMotion(cartSensors);
  cartTimes = cartMotion.times;
else
  cartTimes = cartSensors.times;
end
phasespaceTimes = phasespaceData.time;

%%
timeDifference = cartTimes(end) - phasespaceTimes(end);
if(verbose > 0)
  fprintf('trial cart data file has %g more seconds of data than the phasespace trial (%g phasespace frames)\n', ...
    timeDifference, timeDifference * 480);
end


%% transform phasespace marker to world frame
% [phasespaceCartLocation, angles] = integrateCartSensors(processedPhasespaceTrial.cartSensors);
% cartMotion = [phasespaceCartLocation angles];
% cartSpline = spline(cartTimes + trialTimeOffset, cartMotion', phasespaceTimes);

% transformedMarkerFile = [processedDataFolderName rootOfSavedFilename '_transformedMarkerData.mat'];

markersInWorldFrame = transformPhasespaceRecordIntoWorldWithMobileCartAlignment(...
  phasespaceData, ...
  cartSensors, calibration);
% save(transformedMarkerFile, 'markersInWorldFrame');

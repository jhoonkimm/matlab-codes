function [alignedPhasespaceRecord] = ...
  transformPhasespaceRecordIntoWorldWithMobileCartAlignment( ...
  phasespaceRecord, mobileCartMotion, phasespaceScaleAndAlignment, varargin)

forceReprocessSourceFiles = 0;
verbose = 0;
for i = 1:2:length(varargin)
  if (strcmp('forceReprocessSourceFiles', varargin{i}))
    forceReprocessSourceFiles = varargin{i + 1};
  end
  if (strcmp('verbose', varargin{i}))
    verbose = varargin{i + 1};
  end
end

%%
if (ischar(phasespaceRecord))
  phasespaceRecord = loadPhasespaceRecord(phasespaceRecord, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles);
end
if ((ischar(mobileCartMotion)))
  mobileCartRecord = loadMobileCartRecord(mobileCartMotion, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles);
  [mobileCartMotion] = integrateMobileCartMotion(mobileCartRecord);
end

if (~isfield(mobileCartMotion, 'positions'))
  [mobileCartMotion] = integrateMobileCartMotion(mobileCartMotion);
end

%%
trialTimeOffset = 0.0;

cartMotion = [mobileCartMotion.positions, mobileCartMotion.angles];

numMarkerTicks = length(phasespaceRecord.(phasespaceRecord.markerNames{1}));
phasespaceTimes = linspace(0, (numMarkerTicks - 1)/phasespaceRecord.frameRate, numMarkerTicks);
cartTimes = mobileCartMotion.times;

cartSpline = spline(cartTimes + trialTimeOffset, cartMotion', phasespaceTimes);

alignedPhasespaceRecord = phasespaceRecord;
markerNames = alignedPhasespaceRecord.markerNames;
for markerNumber = 1:length(markerNames)
  markerName = markerNames{markerNumber};
  markerDataToUse = alignedPhasespaceRecord.(markerName);
  markerDataInWorld = zeros(size(markerDataToUse));
  for i = 1:(length(phasespaceTimes) / 1)
    currentPhasespacePoint = markerDataToUse(i, 1:2)';
    cartPose = cartSpline(:,i);
    markerDataInWorld(i, 1:2) = ...
      transformPhasespacePointIntoWorld(currentPhasespacePoint, ...
      phasespaceScaleAndAlignment(2:4), phasespaceScaleAndAlignment(1), cartPose);
    markerDataInWorld(i, 3) = markerDataToUse(i, 3);
  end
  alignedPhasespaceRecord.(markerName) = markerDataInWorld;
end


end


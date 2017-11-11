function [] = makeVideoOfPhasespaceRecord(movieFilename, phasespaceRecord, varargin)

forceReprocessSourceFiles = 0;
mobileCartMotion = [];
frameSkip = floor(480/15);
recordName = '';
verbose = 0;
trackMarkers = 0;
for i = 1:2:length(varargin)
  %   if (strcmp('forceReprocess', varargin{i}))
  %     forceReprocess = varargin{i + 1};
  %   end
  if (strcmp('forceReprocessSourceFiles', varargin{i}))
    forceReprocessSourceFiles = varargin{i + 1};
  end
  if (strcmp('verbose', varargin{i}))
    verbose = varargin{i + 1};
  end
  if (strcmp('mobileCartMotion', varargin{i}))
    mobileCartMotion = varargin{i + 1};
  end
  if (strcmp('frameSkip', varargin{i}))
    frameSkip = varargin{i + 1};
  end
  if (strcmp('recordName', varargin{i}))
    recordName = varargin{i + 1};
  end
  if (strcmp('trackMarkers', varargin{i}))
    trackMarkers = varargin{i + 1};
  end
end

%% check if we inputted data records or files with the data in them
if (ischar(phasespaceRecord))
  phasespaceRecord = loadPhasespaceRecord(phasespaceRecord, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles);
end
if (~isempty(mobileCartMotion))
  if ((ischar(mobileCartMotion)))
    mobileCartRecord = loadMobileCartRecord(mobileCartMotion, 'verbose', verbose, 'forceReprocess', forceReprocessSourceFiles);
    [mobileCartMotion] = integrateMobileCartMotion(mobileCartRecord);
  end
  if (~isfield(mobileCartMotion, 'positions'))
    [mobileCartMotion] = integrateMobileCartMotion(mobileCartMotion);
  end
end


%% determine camera target

% for frameNumber = 1:frameSkip:length(phasespaceRecord.(phasespaceRecord.markerNames{1}))
cameraDistance = 2;
cameraTarget = zeros(3, 1);
numMarkers = length(phasespaceRecord.markerNames);

XLim = [Inf -Inf];
YLim = [Inf -Inf];
ZLim = [Inf -Inf];
if (~isempty(mobileCartMotion))
  XLim = [min(mobileCartMotion.positions(:, 1)) max(mobileCartMotion.positions(:, 1))];
  YLim = [min(mobileCartMotion.positions(:, 2)) max(mobileCartMotion.positions(:, 2))];
end

for markerNumber = 1:numMarkers;
  marker = phasespaceRecord.(phasespaceRecord.markerNames{markerNumber});
  %   cameraTarget = cameraTarget + (nanmean(marker)' / numMarkers);
  mins = nanmin(marker);
  maxs = nanmax(marker);
  XLim(1) = min([XLim(1) mins(1)]);
  YLim(1) = min([YLim(1) mins(2)]);
  ZLim(1) = min([ZLim(1) mins(3)]);
  XLim(2) = max([XLim(2) maxs(1)]);
  YLim(2) = max([YLim(2) maxs(2)]);
  ZLim(2) = max([ZLim(2) maxs(3)]);
end
cameraTarget = mean([XLim; YLim; ZLim]')';
cameraPosition = cameraTarget + [1, 1, 1]' * cameraDistance / norm([1, 1, 1]);

arrowLength = max(diff([XLim; YLim; ZLim]')) / 20.0;
% end

%% save out movie
fig = gcf;
aviobj = avifile(movieFilename);
tic
for frameNumber = 1:frameSkip:length(phasespaceRecord.(phasespaceRecord.markerNames{1}))
  cla;
    plotSingleFrameOfPhasespaceRecord(phasespaceRecord, frameNumber, ...
      'mobileCartMotion', mobileCartMotion, ...
      'cameraPosition', cameraPosition, 'cameraTarget', cameraTarget, ...
      'XLim', XLim, 'YLim', YLim, 'ZLim', ZLim, ...
      'arrowLength', arrowLength, 'trackMarkers', trackMarkers);
  
%   numMarkers = length(phasespaceRecord.markerNames);
%   colors = lines(numMarkers);
%   for markerNumber = 1:numMarkers
%     markerName = phasespaceRecord.markerNames{markerNumber};
%     %   thisMarker = phasespaceRecord.(markerName);
%     
%     plot3(phasespaceRecord.(markerName)(frameNumber, 1), phasespaceRecord.(markerName)(frameNumber, 2), phasespaceRecord.(markerName)(frameNumber, 3), ...
%       'o', 'MarkerSize', 12, 'MarkerFaceColor', colors(markerNumber, :)); %, ...
%     %     'color', colors(markerNumber, :));
%     hold on;
%     axis equal;
%   end
%   title(sprintf('phasespace record %s\n%g markers, %gHz, %g samples, t = %gs, frame %g', ...
%     recordName, numMarkers, phasespaceRecord.frameRate, length(phasespaceRecord.(markerName)), ...
%     length(phasespaceRecord.(markerName)) / phasespaceRecord.frameRate, frameNumber));
%   %   legend(phasespaceRecord.markerNames);
%   
%   fprintf('frameNumber = %g\n', frameNumber);
  
  frame = getframe(fig);
  aviobj = addframe(aviobj, frame);
end
toc
aviobj = close(aviobj);


end


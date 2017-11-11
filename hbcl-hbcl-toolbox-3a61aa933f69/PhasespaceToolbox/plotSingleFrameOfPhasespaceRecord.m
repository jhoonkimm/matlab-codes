function [] = plotSingleFrameOfPhasespaceRecord(phasespaceRecord, frameNumber, varargin)

mobileCartMotion = [];
recordName = '';
cameraTarget = [];
cameraPosition = [];
XLim = [];
YLim = [];
ZLim = [];
arrowLength = [];
trackMarkers = 0;
for i = 1:2:length(varargin)
  if (strcmp('mobileCartMotion', varargin{i}))
    mobileCartMotion = varargin{i + 1};
  end
  if (strcmp('recordName', varargin{i}))
    recordName = varargin{i + 1};
  end
  if (strcmp('cameraTarget', varargin{i}))
    cameraTarget = varargin{i + 1};
  end
  if (strcmp('cameraPosition', varargin{i}))
    cameraPosition = varargin{i + 1};
  end
  if (strcmp('XLim', varargin{i}))
    XLim = varargin{i + 1};
  end
  if (strcmp('YLim', varargin{i}))
    YLim = varargin{i + 1};
  end
  if (strcmp('ZLim', varargin{i}))
    ZLim = varargin{i + 1};
  end
  if (strcmp('arrowLength', varargin{i}))
    arrowLength = varargin{i + 1};
  end
  if (strcmp('trackMarkers', varargin{i}))
    trackMarkers = varargin{i + 1};
  end
end

numMarkers = length(phasespaceRecord.markerNames);
colors = lines(numMarkers);


positions = [];
for markerNumber = 1:numMarkers
  markerName = phasespaceRecord.markerNames{markerNumber};
  %   thisMarker = phasespaceRecord.(markerName);
  
  plot3(phasespaceRecord.(markerName)(frameNumber, 1), phasespaceRecord.(markerName)(frameNumber, 2), phasespaceRecord.(markerName)(frameNumber, 3), ...
    'o', 'MarkerSize', 12, 'MarkerFaceColor', colors(markerNumber, :)); %, ...
  hold on;
  
  positions = [positions; [phasespaceRecord.(markerName)(frameNumber, 1:3)]];
  
end

%%

if (~isempty(mobileCartMotion))
  
  if (isempty(arrowLength))
    motionBounds = max(max(mobileCartMotion.positions) - min(mobileCartMotion.positions));
    arrowLength = motionBounds / 20;
  end
  numMarkerTicks = length(phasespaceRecord.(phasespaceRecord.markerNames{1}));
  phasespaceTimes = linspace(0, (numMarkerTicks - 1)/phasespaceRecord.frameRate, numMarkerTicks);
  cartTimes = mobileCartMotion.times;
  
  position = interp1(cartTimes, mobileCartMotion.positions, phasespaceTimes(frameNumber));
  angle = interp1(cartTimes, mobileCartMotion.angles, phasespaceTimes(frameNumber));
  
  endPoint = position' + arrowLength * [cos(angle) sin(angle)]';
  line([position(1) endPoint(1)], [position(2) endPoint(2)], 'color', 'b', 'LineWidth', 4);
  plot(position(1), position(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

%%
axis equal;
if (~isempty(XLim))
  set(gca, 'XLim', XLim);
end
if (~isempty(YLim))
  set(gca, 'YLim', YLim);
end
if (~isempty(ZLim))
  set(gca, 'ZLim', ZLim);
end

% set(gca, 'XGrid', 'on');
% set(gca, 'YGrid', 'on');
% set(gca, 'Projection', 'perspective');


if (~isempty(cameraPosition))
  set(gca, 'CameraPosition', cameraPosition);
end
if (~isempty(cameraTarget))
  set(gca, 'CameraTarget', cameraTarget);
end
if (trackMarkers)
  camva(45);
  if (numel(positions) > 9)
    try
    set(gca, 'CameraPosition', nanmean(positions));
    set(gca, 'CameraTarget', nanmean(positions) + [1 1 1] * 4/norm([1,1,1]));
    catch e
      %       e
    end
  end
end

% set(gca, 'Projection', 'perspective');
% camproj('perspective')

title(sprintf('phasespace record %s\n%g markers, %gHz, %g samples, t = %gs, frame %g', ...
  recordName, numMarkers, phasespaceRecord.frameRate, length(phasespaceRecord.(markerName)), ...
  length(phasespaceRecord.(markerName)) / phasespaceRecord.frameRate, frameNumber));
% legend(phasespaceRecord.markerNames);

end


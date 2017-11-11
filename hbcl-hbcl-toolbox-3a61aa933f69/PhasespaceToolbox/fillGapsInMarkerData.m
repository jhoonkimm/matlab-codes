function [filledMarkerData] = fillGapsInMarkerData(markerName, markerData, varargin)

askUserToApproveFill = 1;
for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch option
    case 'askUserToApproveFill'
      askUserToApproveFill = value;
  end
end

%%
clf;
plot(markerData)
hold on;

xCoordinate = markerData(:, 1);
goodIndeces = ~isnan(xCoordinate);
allIndeces = 1:length(xCoordinate);
try
  splinedMarkerData = spline(allIndeces(goodIndeces), markerData(goodIndeces, :)', allIndeces)';
  
  firstGoodIndex = find(goodIndeces, 1, 'first');
  indecesToRelace = 1: firstGoodIndex;
  splinedMarkerData(indecesToRelace, :) = repmat(splinedMarkerData(firstGoodIndex, :), ...
    [length(indecesToRelace), 1]);
  
  lastGoodIndex = find(goodIndeces, 1, 'last');
  indecesToRelace = lastGoodIndex : size(splinedMarkerData, 1);
  splinedMarkerData(indecesToRelace, :) = repmat(splinedMarkerData(lastGoodIndex, :), ...
    [length(indecesToRelace), 1]);
  
  %   splinedMarkerData(splinedMarkerData > 10) = 10000;
  %   splinedMarkerData(splinedMarkerData < -10) = -10000;
catch e
  e
end

onlySplinedMarkerData = splinedMarkerData;
onlySplinedMarkerData(goodIndeces, :) = NaN;
color = [0.9 0.9 0.2];
plot(allIndeces, onlySplinedMarkerData, 'LineWidth', 3, 'color', color);

title(sprintf('marker: %s', markerName));
legend('x', 'y', 'z', 'proposed spline fills');
box off;

filledMarkerData = splinedMarkerData;

if (askUserToApproveFill)
  input(['Here are the proposed splines for this marker. If you don''t like them, you have to fix them yourself.\n' ...
    'If you don''t want to have the data splined, and wish to write out a c3d with missing data, add: ..., ''shouldFillGapsInMarkerData'', 0)\n' ...
    'to the call to savePhasespaceAndAnalogDataInC3D\n' ...
    'Press enter to accept']);
end


end


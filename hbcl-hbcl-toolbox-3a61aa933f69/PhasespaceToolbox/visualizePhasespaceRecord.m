function [] = visualizePhasespaceRecord(phasespaceRecord, varargin)

recordName = '';
plotInterval = 10;
plotDropoutRateHistogram = 1;
axisHandlesToPlotOn = [];
for i = 1:2:length(varargin)
  if (strcmp('recordName', varargin{i}))
    recordName = varargin{i + 1};
  end
  if (strcmp('plotInterval', varargin{i}))
    plotInterval = varargin{i + 1};
  end
  if (strcmp('plotDropoutRateHistogram', varargin{i}))
    plotDropoutRateHistogram = varargin{i + 1};
  end
  if (strcmp('axisHandlesToPlotOn', varargin{i}))
    axisHandlesToPlotOn = varargin{i + 1};
  end
end

%% first figure out or create the axes to plot on
if (plotDropoutRateHistogram)
  if (~isempty(axisHandlesToPlotOn))
    markerAxis = axisHandlesToPlotOn(1);
    dropoutHistogramAxis = axisHandlesToPlotOn(2);
  else
    markerAxis = subplot(3, 1, 1:2);
    dropoutHistogramAxis = subplot(3, 1, 3);
  end
else
  markerAxis = gca;
end

%%
axes(markerAxis);
numMarkers = length(phasespaceRecord.markerNames);
dropoutRates = zeros(1, numMarkers);
colors = lines(numMarkers);
for markerNumber = 1:numMarkers
  markerName = phasespaceRecord.markerNames{markerNumber};
  thisMarker = phasespaceRecord.(markerName);
  
  dropoutRates(markerNumber) = sum(isnan(thisMarker(:, 1))) / length(thisMarker);
  
  indecesToPlot = 1:plotInterval:length(thisMarker);
  plot3(thisMarker(indecesToPlot, 1), thisMarker(indecesToPlot, 2), thisMarker(indecesToPlot, 3), 'color', colors(markerNumber, :));
  hold on;
  axis equal;
end
title(sprintf('phasespace record %s plot interval = %g\n%g markers, %gHz, %g samples, t = %gs', ...
  recordName, plotInterval, ...
  numMarkers, phasespaceRecord.frameRate, length(thisMarker), length(thisMarker) / phasespaceRecord.frameRate));
legend(phasespaceRecord.markerNames);

%%
if (plotDropoutRateHistogram)
  axes(dropoutHistogramAxis);
  hist(dropoutRates);
  xlabel('dropout rate');
  ylabel('num markers in bucket');
  title('dropout histogram');
end



end
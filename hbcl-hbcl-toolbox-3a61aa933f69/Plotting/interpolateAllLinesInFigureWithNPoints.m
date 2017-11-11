function [] = interpolateAllLinesInFigureWithNPoints(varargin)
%%
% If, for example, you are exporting a figure to illustrator, and the lines
% in the plots are getting broken into separate illustrator paths (making it 
% annoying to edit), you can use this function to reduce the number of points 
% per line so matlab doesn't break the line into multiple paths.

figureHandle = gcf;
numberOfPointsToPlotPerLine = 100;

%%
for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch (option)
    case 'figureHandle'
      figureHandle = value;
    case 'numberOfPointsToPlotPerLine'
      numberOfPointsToPlotPerLine = value;
  end
end

%%

hLine = findobj(figureHandle, 'type', 'line');
for i = 1:length(hLine)
  h = hLine(i);
  xData = get(h, 'XData');
  yData = get(h, 'YData');
  
  if (length(xData) == 1)
    continue;
  end
  
  allIndeces = 1 : length(xData);
  xIndecesToPlot = linspace(1, allIndeces(end), numberOfPointsToPlotPerLine);
  xData = interp1(allIndeces, xData, xIndecesToPlot);
  yData = interp1(allIndeces, yData, xIndecesToPlot);
  
  set(h, 'XData', xData);
  set(h, 'YData', yData);
end

end


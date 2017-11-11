function [middlePointOfAxisLine] = drawAxes(axisDimension, axisLimits, yVal, tickLocations, tickHeights, varargin)

axisColor = 'k';
axisLineWidth = 1;

tickColor = 'k';
tickLineWidth = 0.5;

tickLabelTextOptions = {};
tickDirection = 1;

tickLabels = [];
centerFirstAndLastTickLabels = 0;

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch(option)
    case 'axisColor'
      axisColor = value;
    case 'tickLineWidth'
      tickLineWidth = value;
    case 'tickColor'
      tickColor = value;
    case 'axisLineWidth'
      axisLineWidth = value;
    case 'tickLabelTextOptions'
      tickLabelTextOptions = value;
    case 'tickDirection'
      tickDirection = value;
    case 'tickLabels'
      tickLabels = value;
    case 'centerFirstAndLastTickLabels'
      centerFirstAndLastTickLabels = value;
  end
end

%% expand out arguments if necessary
if (strcmp(axisDimension, 'x'))
  axisDimensionNumber = 1;
elseif (strcmp(axisDimension, 'y'))
  axisDimensionNumber = 2;
end

if (length(tickHeights) == 1)
  tickHeights = tickHeights * ones(size(tickLocations));
end

if (isempty(tickLabels))
  for i = 1:length(tickLocations)
    tickLabels{i} = num2str(tickLocations(i));
  end
end

%% draw the axis and ticks
middlePointOfAxisLine = drawLineInAxisDimension(axisLimits, [yVal yVal], axisDimensionNumber, 'color', axisColor, 'LineWidth', axisLineWidth);

for i = 1 : length(tickLocations)
  
  isFirstTick = i == 1;
  isLastTick = i == length(tickLocations);
  
  if (centerFirstAndLastTickLabels)
    isFirstTick = 0;
    isLastTick = 0;
  end
  
  if (tickDirection == 0)
    tickYVals = [yVal - tickHeights(i)/2 yVal + tickHeights(i)/2];
    %     drawLineInAxisDimension([yVal - tickHeights(i)/2 yVal + tickHeights(i)/2], tickLocations(i), ...
    %       axisDimensionNumber, 'color', tickColor, 'LineWidth', tickLineWidth);
  end
  if (tickDirection == 1)
    tickYVals = [yVal yVal + tickHeights(i)];
    %     drawLineInAxisDimension([yVal yVal + tickHeights(i)], tickLocations(i), ...
    %       axisDimensionNumber, 'color', tickColor, 'LineWidth', tickLineWidth);
  end
  if (tickDirection == -1)
    tickYVals = [yVal yVal - tickHeights(i)];
    %     drawLineInAxisDimension([yVal yVal - tickHeights(i)], tickLocations(i), ...
    %       axisDimensionNumber, 'color', tickColor, 'LineWidth', tickLineWidth);
  end
  
  drawLineInAxisDimension([tickLocations(i) tickLocations(i)], tickYVals, ...
    axisDimensionNumber, 'color', tickColor, 'LineWidth', tickLineWidth);
  
  
  if (length(tickLabelTextOptions) > 0)
    %     drawTextInAxisDimension(tickYVals(1), tickLocations(i), axisDimensionNumber, ...
    %       tickLabels{i}, isFirstTick, isLastTick, tickLabelTextOptions);
    drawTextInAxisDimension(tickLocations(i), tickYVals(1), axisDimensionNumber, ...
      tickLabels{i}, isFirstTick, isLastTick, tickLabelTextOptions{:});
  else
    drawTextInAxisDimension(tickLocations(i), tickYVals(1), axisDimensionNumber, ...
      tickLabels{i}, isFirstTick, isLastTick);
  end
  
  
end


end


function [middlePointOfLine] = drawLineInAxisDimension(xVals, yVals, axisDimensionNumber, varargin)
if (axisDimensionNumber == 1)
  line(xVals, yVals, 'clipping', 'off', varargin{:})
  middlePointOfLine = [mean(xVals), mean(yVals)];
else
  line(yVals, xVals, 'clipping', 'off', varargin{:})
  middlePointOfLine = [mean(yVals), mean(xVals)];
end
end

function [] = drawTextInAxisDimension(xVal, yVal, axisDimensionNumber, textToDraw, isFirstTick, isLastTick, varargin)

if length(varargin) == 1
  varargin = varargin{:};
end

fontSize = 10;
interpreter = 'tex';

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch(option)
    case 'fontSize'
      fontSize = value;
    case 'FontSize'
      fontSize = value;
    case 'interpreter'
      interpreter = value;
    case 'Interpreter'
      interpreter = value;
  end
end


if (axisDimensionNumber == 1)
  horizontalAlignment = 'center';
  verticalAlignment = 'top';
  if (isFirstTick)
    horizontalAlignment = 'left';
  elseif (isLastTick)
    horizontalAlignment = 'right';
  end
  
  
  for i = 1 : 2 : length(varargin)
    option = varargin{i};
    value = varargin{i + 1};
    switch(option)
      case 'horizontalAlignment'
        horizontalAlignment = value;
      case 'verticalAlignment'
        verticalAlignment = value;
    end
  end
  
  if (iscell(textToDraw))
    for i = 1 : length(textToDraw)
      h = text(xVal, yVal, textToDraw{i}, 'HorizontalAlignment', horizontalAlignment, 'VerticalAlignment', verticalAlignment, ...
        'FontSize', fontSize, 'interpreter', interpreter);
      height = get(h, 'Extent');
      height = height(4);
      height = height * 0.9;
      yVal = yVal - height;
    end
  else
    h = text(xVal, yVal, textToDraw, 'HorizontalAlignment', horizontalAlignment, 'VerticalAlignment', verticalAlignment, ...
      'FontSize', fontSize, 'interpreter', interpreter);
  end
  
else
  verticalAlignment = 'middle';
  horizontalAlignment = 'right';
  if (isFirstTick)
    verticalAlignment = 'baseline';
  elseif (isLastTick)
    verticalAlignment = 'cap';
  end
  
  for i = 1 : 2 : length(varargin)
    option = varargin{i};
    value = varargin{i + 1};
    switch(option)
      case 'horizontalAlignment'
        horizontalAlignment = value;
      case 'verticalAlignment'
        verticalAlignment = value;
    end
  end
  
  %   try
  h = text(double(yVal), double(xVal), textToDraw, ...
    'HorizontalAlignment', horizontalAlignment, ...
    'VerticalAlignment', verticalAlignment, ...
    'FontSize', fontSize, 'interpreter', interpreter);
  %   catch e
  %     e
  %   end
  
  %   textBoxPosition = dsxy2figxy([yVal xVal, 1e-3, 1e-3]);
  %   annotation('textbox', textBoxPosition, 'FitBoxToText', 'off', 'string', textToDraw, ...);
  %   'HorizontalAlignment', horizontalAlignment, 'VerticalAlignment', verticalAlignment, ...
  %     'FontSize', fontSize, 'interpreter', interpreter);
  
end
% if (length(varargin) > 0)
%   otherArgs = varargin{:};
%   set(h, otherArgs{1}, otherArgs{2});
% end

end

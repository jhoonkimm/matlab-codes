function axes3x3 = plot3x3( ...
  AnkleAngles, AnkleMoment, AnklePower, ...
  KneeAngles, KneeMoment, KneePower, ...
  HipAngles, HipMoment, HipPower, ...
  varargin)

% PLOT3X3 plot angles, moments and powers for ankle, knee and hip
%  Plot standard sagittal plane joint "3x3": angles, moments and powers
%  for ankle, knee and hip.
%  Convention Note (for angles and moments): EXTENSION is always POSITIVE
%
%  Options:
%   textFontSize: change default text font size
%   axesFontSize: change default axes font size
%   lineWidth: change plot line width
%   lineStyle: change line type ('-','--','.',etc)
%   axesColorOrder: change line color order; if only plotting one
%      curve/condition then this is same as 'Color'
%   xTickLabel: change x tick marks to label
%   xAxisLabel: change x axis name
%   angleUnits: specify angle units
%   momentUnits: specify moment units
%   powerUnits: specify power units
%   labelFlexAndExt: add 'Ext' and 'Flex' labels to plot? (Y/N)
%   plotZeroLine: add line at zero? (Y/N)
%   nonDimConstant: input [Cangle Cmoment Cpower] array to add
%      a 2nd Y axis with dimensionless units.
%      Angledimensionless = Angle*Cangle
%
%  Example calls
%   plot3x3(avgRAnkleAnglesR,avgRAnkleMomentR,avgRAnklePowerR,avgRKneeAnglesR,avgRKneeMomentR,avgRKneePowerR,avgRHipAnglesR,avgRHipMomentR,avgRHipPowerR)
%   plot3x3([avgRAnkleAnglesR avgLAnkleAnglesL],[avgRAnkleMomentR avgLAnkleMomentL],[avgRAnklePowerR avgLAnklePowerL],[avgRKneeAnglesR avgLKneeAnglesL],[avgRKneeMomentR avgLKneeMomentL],[avgRKneePowerR avgLKneePowerL],[avgRHipAnglesR avgLHipAnglesL],[avgRHipMomentR avgLHipMomentL],[avgRHipPowerR avgLHipPowerL])

%  Karl Zelik
%  11/19/09

%% don't look at this:
global gridOfAxes

%% Default plot options
textFontSize = 14;
axesFontSize = 14;
lineWidth = 3;
lineStyle = '-';
figureNum = [];
axesColorOrder = [];
xTickLabel = {'0','50','100'}; %% {'0','25','50','75','100'}; % default plot is heelstrike to heelstrike
xAxisLabel = 'Gait Cycle (%)';
angleUnits = 'degrees';
momentUnits = 'Nm/kg';
powerUnits = 'W/kg';
labelFlexAndExt = 1;
plotZeroLine = 1;
nonDimConstant = []; % input in order angle, moment, power
stdevOn = 0;
trialName = '';
calculateAndPlotStandardDeviations = 0;
precalculatedStandardDeviationsToPlot = [];
gridOfAxes = GridOfAxes(3, 3);
useMinimalAxes = 0; %1; %


plotStandardDeviationLines = 1;
plotStandardDeviationAreas = 1;

%%
gridOfAxes.normalizedRectangleOfAxisGrid = [0.1 0.2 0.8 0.65]; %[0.1 0.1 0.8 0.825]; %
gridOfAxes.horizontalSpacing = 0.05; %
gridOfAxes.verticalSpacing = -0.1; %-0.05; %0.0; %0.025; %0; %

%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'figure'
      figureNum = val;
    case 'textFontSize'
      textFontSize = val;
    case 'axesFontSize'
      axesFontSize = val;
    case 'lineWidth'
      lineWidth = val;
    case 'lineStyle'
      lineStyle = val;
    case 'axesColorOrder'
      axesColorOrder = val;
    case 'xTickLabel'
      xTickLabel = val;
    case 'xAxisLabel'
      xAxisLabel = val;
    case 'angleUnits'
      angleUnits = val;
    case 'momentUnits'
      momentUnits = val;
    case 'powerUnits'
      powerUnits = val;
    case 'labelFlexAndExt'
      labelFlexAndExt = val;
    case 'plotZeroLine'
      plotZeroLine = val;
    case 'nonDimConstant'
      nonDimConstant = val;
    case 'stdevOn'
      stdevOn = val;
    case 'stdevSizes'
      stdevSizes = val;
    case 'stdevColorOrder'
      stdevColorOrder = val;
    case 'trialName'
      trialName = val;
    case 'gridOfAxes'
      gridOfAxes = val;
    case 'calculateAndPlotStandardDeviations'
      calculateAndPlotStandardDeviations = val;
    case 'precalculatedStandardDeviationsToPlot'
      precalculatedStandardDeviationsToPlot = val;
    case 'plotStandardDeviationLines'
      plotStandardDeviationLines = val;
    case 'plotStandardDeviationAreas'
      plotStandardDeviationAreas = val;
    case 'useMinimalAxes'
      useMinimalAxes = val;
    otherwise
      error(['Error, incorrect varargin: ' opt]);
  end
end


%% Create figure
areaHandles = [];
if isempty(figureNum)
  figure();
else
  figure(figureNum);
end
set(gcf, 'color', 'w');
set(gcf, 'Name', trialName);
set(gcf, 'DefaultTextFontSize', textFontSize);
set(gcf, 'DefaultAxesFontSize', axesFontSize);
if isempty(axesColorOrder)
  axesColorOrder = get(0, 'DefaultAxesColorOrder');
else
  set(gca, 'ColorOrder', axesColorOrder);
  if (~isempty(gridOfAxes))
    clf;
  end
end
xTick = linspace(1, size(AnkleAngles,1), size(xTickLabel,2));

mySubplot = @(subplotNumber) activateSubplot(subplotNumber); %, gridOfAxes);

%% Store all joint data in 1 array
if ~stdevOn && ~calculateAndPlotStandardDeviations && isempty(precalculatedStandardDeviationsToPlot)
  data = [AnkleAngles KneeAngles HipAngles AnkleMoment KneeMoment HipMoment AnklePower KneePower HipPower];
  numConditions = size(AnkleAngles,2);
elseif stdevOn
  %   tempdata = [AnkleAngles KneeAngles HipAngles AnkleMoment KneeMoment HipMoment AnklePower KneePower HipPower];
  numConditions = size(stdevSizes,2);
  stdevSizes = [0 stdevSizes];
  for g = 1:numConditions
    stdevIndex = (stdevSizes(g)+1):stdevSizes(g+1);
    data(:,(1:9)*numConditions-(numConditions-g)) = [mean(AnkleAngles(:,stdevIndex), 2) mean(KneeAngles(:,stdevIndex),2) mean(HipAngles(:,stdevIndex),2) mean(AnkleMoment(:,stdevIndex),2) mean(KneeMoment(:,stdevIndex),2) mean(HipMoment(:,stdevIndex),2) mean(AnklePower(:,stdevIndex),2) mean(KneePower(:,stdevIndex),2) mean(HipPower(:,stdevIndex),2)];
    stdev(:,(1:9)*numConditions-(numConditions-g)) = [std(AnkleAngles(:,stdevIndex),0,2) std(KneeAngles(:,stdevIndex),0,2) std(HipAngles(:,stdevIndex),0,2) std(AnkleMoment(:,stdevIndex),0,2) std(KneeMoment(:,stdevIndex),0,2) std(HipMoment(:,stdevIndex),0,2) std(AnklePower(:,stdevIndex),0,2) std(KneePower(:,stdevIndex),0,2) std(HipPower(:,stdevIndex),0,2)];
  end
elseif calculateAndPlotStandardDeviations
  numConditions = size(AnkleAngles,2);
  data = [AnkleAngles KneeAngles HipAngles AnkleMoment KneeMoment HipMoment AnklePower KneePower HipPower];
  stdev = std(data, 0, 3);
  data = mean(data, 3);
elseif ~isempty(precalculatedStandardDeviationsToPlot)
  numConditions = size(AnkleAngles,2);
  data = [AnkleAngles KneeAngles HipAngles AnkleMoment KneeMoment HipMoment AnklePower KneePower HipPower];
  stdev = precalculatedStandardDeviationsToPlot;
  %   data = mean(data, 3);
end



%% Plot zero line
if plotZeroLine
  for j = 1:9
    mySubplot(j);
    ax = axis;
    hold on;
    plot([1 size(data, 1)], [0 0], ...
      'color', ones(3,1) * 0.5, ...
      'LineWidth', 0.5);
    axis(ax);
    hold on;
  end
end

%% Plotting loop
for i = 1 : numConditions
  lineColor = axesColorOrder(i,:); % set line color - different for each condition
  
  for j = 1:9
    eval(['h' num2str(j) '= mySubplot(j);']);
    
    hold on;
    indecesToPlot = j*numConditions-(numConditions-i);
    
    if ~stdevOn && ~calculateAndPlotStandardDeviations && isempty(precalculatedStandardDeviationsToPlot)
      %       plot(data(:,j*numConditions-(numConditions-i)), ...
      %         'LineWidth', lineWidth, ...
      %         'Color', lineColor, ...
      %         'LineStyle',lineStyle);
      
      
    elseif stdevOn
      plotMeanWithStd2(1:size(data,1), ...
        data(:,j*numConditions-(numConditions-i))', ...
        stdev(:,j*numConditions-(numConditions-i))',...
        lineStyle, lineWidth, 0.9-(i-1)*0.3, stdevColorOrder(i,:), lineColor);
      hold on;
      
    elseif calculateAndPlotStandardDeviations || ~isempty(precalculatedStandardDeviationsToPlot)

      if (plotStandardDeviationLines)
        plot(data(:, indecesToPlot) + stdev(:, indecesToPlot), ...
          'LineWidth', lineWidth / 2, ...
          'Color', lineColor, ...
          'LineStyle',lineStyle);
        hold on;
        plot(data(:, indecesToPlot) - stdev(:, indecesToPlot), ...
          'LineWidth', lineWidth / 2, ...
          'Color', lineColor, ...
          'LineStyle',lineStyle);
      end
      if (plotStandardDeviationAreas)
        xVals = [1 : size(data, 1) size(data, 1) : -1 : 1];
        areaColor = lineColor;
        areaHandles(length(areaHandles) + 1) = fill(xVals, ...
          [data(:, indecesToPlot) + stdev(:, indecesToPlot); ...
          flipud(data(:, indecesToPlot) - stdev(:, indecesToPlot))], ...
          areaColor, 'EdgeColor', 'none');
        hold on;
      end
    end
    
    plot(data(:, indecesToPlot), ...
      'LineWidth', lineWidth, ...
      'Color', lineColor, ...
      'LineStyle',lineStyle);
    
    xlim([1, size(data, 1)]);
    box off;
    setYAxisLimits(gca)
  end
end


for j = 1:9
  mySubplot(j);
  children = get(gca, 'children');
  areasHere = intersect(areaHandles, children);
  if (~isempty(areasHere))
    children = setdiff(children, areasHere);
    children = [children; areasHere];
    set(gca, 'children', children);
  end
end


%% Set equivalent axes
setEqualAxes([h1 h2 h3]); % angles
setEqualAxes([h4 h5 h6]); % moments
setEqualAxes([h7 h8 h9]); % powers

for j = [2 3 5 6 8 9]
  mySubplot(j);
  set(gca,'YTickLabel',{''});
end

%% labels, axes, ticks, oh my...
for i = 1 : numConditions
  
  for j = 1:9
    mySubplot(j);
    
    % Label joint names
    if j==1
      title('Ankle');
    end
    if j==2
      title('Knee');
    end
    if j==3
      title('Hip');
    end
    
    % x and y axes, labels
    if (useMinimalAxes)
      axis off;
      
      yAxisLocationAlongX = -25;
      yLabelLocationAlongX = -225; %yAxisLocationAlongX;
      tickLength = -25;
      yLabelProperties = {'rotation', 90, 'verticalAlignment', 'bottom', 'horizontalAlignment', 'center'};
      tickLabelTextOptions = {'fontSize', axesFontSize};
      if j==1
        middlePointOfAxisLine = drawAxes('y', [-45 0], yAxisLocationAlongX, [-45 0], tickLength, ...
          'tickLabelTextOptions', tickLabelTextOptions);
        labelString = sprintf(strcat('Angle\n(',angleUnits,')'));
        text(yLabelLocationAlongX, middlePointOfAxisLine(2), labelString, yLabelProperties{:});
      elseif j==4
        middlePointOfAxisLine = drawAxes('y', [0 100], yAxisLocationAlongX, [0 100], tickLength, ...
          'tickLabelTextOptions', tickLabelTextOptions);
        labelString = sprintf(strcat('Moment\n(',momentUnits,')'));
        text(yLabelLocationAlongX, middlePointOfAxisLine(2), labelString, yLabelProperties{:});
      elseif j==7
        middlePointOfAxisLine = drawAxes('y', [0 200], yAxisLocationAlongX, [0 200], tickLength, ...
          'tickLabelTextOptions', tickLabelTextOptions);
        labelString = sprintf(strcat('Power\n(',powerUnits,')'));
        text(yLabelLocationAlongX, middlePointOfAxisLine(2), labelString, yLabelProperties{:});
      end
      
      if ~isempty(nonDimConstant)
        
        yLabelProperties = {'rotation', 90, ...
          'verticalAlignment', 'top', ...
          'horizontalAlignment', 'center'};

        tickLabelTextOptions = {'fontSize', axesFontSize, ...
          ...'verticalAlignment', 'middle' ...
          'horizontalAlignment', 'left', ...
          };

        if j==3
          
          middlePointOfAxisLine = drawAxes('y', [-45 0], 1000 - yAxisLocationAlongX, [-45 0], -tickLength, ...
            'tickLabelTextOptions', tickLabelTextOptions, 'tickLabels', {'-pi/4', '0'});
          labelString = sprintf(strcat('Radians'));
          text(1000 - yLabelLocationAlongX, middlePointOfAxisLine(2), labelString, yLabelProperties{:});
        elseif j==6
          
          nonDimYAxisHeight = 0.1;
          middlePointOfAxisLine = drawAxes('y', ...
            [0 nonDimConstant(2) * nonDimYAxisHeight], ...
            1000 - yAxisLocationAlongX, ...
            [0 nonDimConstant(2) * nonDimYAxisHeight], ...
            -tickLength, ...
            'tickLabelTextOptions', tickLabelTextOptions, 'tickLabels', {'0', num2str(nonDimYAxisHeight)});
          labelString = sprintf('Nondim.\nTorque');
          text(1000 - yLabelLocationAlongX, double(middlePointOfAxisLine(2)), labelString, yLabelProperties{:});
        elseif j==9
          
          nonDimYAxisHeight = 0.05;
          middlePointOfAxisLine = drawAxes('y', ...
            [0 nonDimConstant(3) * nonDimYAxisHeight], ...
            1000 - yAxisLocationAlongX, ...
            [0 nonDimConstant(3) * nonDimYAxisHeight], ...
            -tickLength, ...
            'tickLabelTextOptions', tickLabelTextOptions, 'tickLabels', {'0', num2str(nonDimYAxisHeight)});
          labelString = sprintf('Nondim.\nPower');
          text(1000 - yLabelLocationAlongX, double(middlePointOfAxisLine(2)), labelString, yLabelProperties{:});
        end
      end
      
      tickLengthMultipleofRange = -0.025;
      labelDropDistanceMultipleofRange = 0.085;
      
      if j==7 | j==8 | j==9
        ylimits = ylim();
        yLimExcursion = diff(ylimits);
        tickLabelTextOptions = {'fontSize', axesFontSize};
        middlePointOfAxisLine = drawAxes('x', [0 1000], ylimits(1), [0 1000], yLimExcursion * tickLengthMultipleofRange, ...
          'tickLabelTextOptions', tickLabelTextOptions, 'tickLabels', {'0' '100'});
      end
      
      xLabelProperties = {'verticalAlignment', 'top', 'horizontalAlignment', 'center'};
      if j==8
        text(middlePointOfAxisLine(1), ...
          middlePointOfAxisLine(2) - yLimExcursion * labelDropDistanceMultipleofRange, ...
          xAxisLabel, xLabelProperties{:});
      end
      
    else
      % y-axis labels
      if j==1
        ylabel(strcat('Angle (',angleUnits,')'));
      elseif j==4
        ylabel(strcat('Moment (',momentUnits,')'));
      elseif j==7
        ylabel(strcat('Power (',powerUnits,')'));
      end
      
      % x axis label and ticks
      if j==8
        xlabel(xAxisLabel);
      end
      
      set(gca,'XTick',xTick);
      if j==7 | j==8 | j==9
        set(gca, 'XTickLabel', xTickLabel);
      else
        set(gca,'XTickLabel',{''});
      end
      
    end
  end
  axes3x3 = [h1 h2 h3 h4 h5 h6 h7 h8 h9];
end



%% Add Flexion & Extension text
% Angles
if ~isempty(nonDimConstant)
  mySubplot(1); hold on;
else
  mySubplot(3); hold on;
end
ax = axis;
ys = linspace(ax(3),ax(4),9);
if labelFlexAndExt
  text(ax(2)*1.05,ys(2),'Flex','FontSize',10);
  text(ax(2)*1.05,ys(8),'Ext','FontSize',10);
end

% Moments
if ~isempty(nonDimConstant)
  mySubplot(4); hold on;
else
  mySubplot(6); hold on;
end
ax = axis;
ys = linspace(ax(3),ax(4),9);
if labelFlexAndExt
  text(ax(2)*1.05,ys(2),'Flex','FontSize',10);
  text(ax(2)*1.05,ys(8),'Ext','FontSize',10);
end


%% Add dimensionless axes
if ~isempty(nonDimConstant)

  
  if (useMinimalAxes)
    %     axes(AX3(2));
    %     axis off;
    %     axes(AX6(2));
    %     axis off;
    %     axes(AX9(2));
    %     axis off;
    
%     mySubplot(3);
    

    
    
  else
    
    mySubplot(3);
    hold on;
    AX3 = addDimensionlessYAxis(gca,nonDimConstant(1)); % angle
    box off;
    
    mySubplot(6);
    hold on;
    AX6 = addDimensionlessYAxis(gca,nonDimConstant(2)); % moment
    box off;
    
    mySubplot(9);
    hold on;
    AX9 = addDimensionlessYAxis(gca,nonDimConstant(3)); % power
    box off;
    
    set(AX3(2),'XTickLabel',{''});
    set(AX6(2),'XTickLabel',{''});
    set(AX9(2),'XTickLabel',{''});
    set(gca,'XTickLabel',xTickLabel);
  end
  
end


end



function [handle] = activateSubplot(plotNumber)
%%
global gridOfAxes

if (isempty(gridOfAxes))
  handle = subplot(3, 3, plotNumber);
else
  [gridOfAxes, handle] = gridOfAxes.subplotFromIndex(plotNumber);
end


end


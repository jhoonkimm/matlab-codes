classdef BarPlot
  
  properties
    
    drawXAxis = 1;
    drawYAxis = 1;
    plotDatapointsWithBars = 1;
    
    barWidth = 0.5;
    
    xTickLocations = []; % default, 1:length(yData)
    xLabels = []; % default, xTickLocations
    xLabelsProperties = {'HorizontalAlignment', 'center'};
    
    barBaseline = 0;
    barColors = ones(3,1) * 0.2;
    
    whiskerColor = 'k';
    whiskerWidth = []; % default, some fraction of the barWidth
    
    legend = {};
    
    yAxisLabel = '';
    yAxisLabelProperties = {};
    yAxisTickLocations = []; % default, [min max] of all ydata
    yAxisTickStrings = {}; % default, [min max] of all ydata
    yAxisTickHeights = []; % default, some fraction of barWidth
    yAxisXPosition = []; % default, a fraction of barwidth to the left of the leftmost bar edge
    
    whiskerLengthsToUse = [];
    
    drawSignificantlyDifferentFromZeroStars = 0;
    
  end
  
  properties (Access = private)
    yData = [];
  end
  
  methods
    
    function [this] = BarPlot()
    end
    
    function [xTickLocations] = getXTickLocations(this)
      %%
      xTickLocations = this.xTickLocations;
      if (isempty(xTickLocations))
        xTickLocations = 1:length(this.yData);
      end
    end
    
    function [this] = plotBars(this, yData, varargin)
      %%
      this.yData = yData;
      [xTickLocations] = this.getXTickLocations();
      
      for i = 1 : length(yData)
        whiskerLength = nanstd(yData{i});
        
        if this.drawSignificantlyDifferentFromZeroStars
          if (ttest(yData{i}, 0, 0.01))
            text(xTickLocations(i), max(yData{i}) * 1.1, '*', 'HorizontalAlignment', 'center');
          end
        end
        
        if (~isempty(this.whiskerLengthsToUse))
          whiskerLength = this.whiskerLengthsToUse(i);
        end
        
        barColor = this.barColors(:, mod(i-1, size(this.barColors, 2)) + 1);
        
        barHandle = this.plotBar(xTickLocations(i), yData{i}, whiskerLength, barColor);
        
        %         %         hp = findobj(h2,'type','patch'); % findobj is highly customizable
        %         hh = hatchfill(barHandle,'cross',45,5,[0.8 0.8 0.8]);
        %         set(hh,'color','b','linewidth',2);

        
      end
    end
    
    function [this] = formatPlot(this, varargin)
      [xTickLocations] = this.getXTickLocations();
      
      if (isempty(this.xLabels))
        xLabels = xTickLocations;
      else
        xLabels = this.xLabels;
      end
      
      %%
      box off;
      axis off;
      
      if (this.drawXAxis)
        tickLocations = xTickLocations;
        tickHeights = 0;
        
        yLims = get(gca, 'ylim');
        axisHeight = min(yLims);
        
        xLabels = this.toStringCellArray(xLabels);
        
        drawAxes('x', [xTickLocations(1) xTickLocations(end)] + [-1 1] * this.barWidth*0.5, ...
          axisHeight, tickLocations, tickHeights, 'tickLabels', xLabels, 'centerFirstAndLastTickLabels', 1, ...
          'tickLabelTextOptions', this.xLabelsProperties); %, varargin)
      end
      
      this = this.plotTheYAxis();
      
    end
    
    function [this] = plotTheYAxis(this)
      %%
      
      if (this.drawYAxis)
        tickLocations = this.yAxisTickLocations;
        if (isempty(tickLocations))
          allYVals = [cat(2, this.yData{:}) this.barBaseline];
          tickLocations = [min(allYVals), max(allYVals)];
        end
        tickLabels = this.toStringCellArray(tickLocations);
        if (~isempty(this.yAxisTickStrings))
          tickLabels = this.yAxisTickStrings;
        end
        
        tickHeights = this.yAxisTickHeights;
        if (isempty(tickHeights))
          tickHeights = -this.barWidth * 0.1;
        end
        
        yAxisXPosition = this.yAxisXPosition;
        if (isempty(yAxisXPosition))
          %           xLims = get(gca, 'xlim');
          %           yAxisXPosition = min(xLims);
          yAxisXPosition = min(this.getXTickLocations()) - this.barWidth/2 - this.barWidth * 0.1;
        end
        
        middlePointOfAxisLine = drawAxes('y', [tickLocations(1) tickLocations(end)], ...
          yAxisXPosition, tickLocations, tickHeights, 'tickLabels', tickLabels);
        if (~isempty(this.yAxisLabel))
          allTextOptions = [{'rotation', 90, 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom'}, ...
            this.yAxisLabelProperties{:}];
          text(middlePointOfAxisLine(1) - this.barWidth * 0.1, middlePointOfAxisLine(2), this.yAxisLabel, allTextOptions{:});
        end
      end
      
      [xLims, yLims] = this.getXAndYLimitsOfPlotData();
      xlim(xLims);
      ylim(yLims);
    end
    
    
    function [xLims, yLims] = getXAndYLimitsOfPlotData(this)
      %%
      dataObjs = get(gca, 'Children');
      objTypes = get(dataObjs, 'Type');
      
      xData = [];
      yData = [];
      for i = 1 : length(objTypes)
        if (strcmp(objTypes{i}, 'line') || strcmp(objTypes{i}, 'patch'))
          
          newXs = get(dataObjs(i), 'XData');
          if (isempty(xData))
            xData = [min(min(newXs)) max(max(newXs))];
          end
          xData(1) = min(xData(1), min(min(newXs)));
          xData(2) = max(xData(2), max(max(newXs)));
          
          newYs = get(dataObjs(i), 'YData');
          if (isempty(yData))
            yData = [min(min(newYs)) max(max(newYs))];
          end
          yData(1) = min(yData(1), min(min(newYs)));
          yData(2) = max(yData(2), max(max(newYs)));
        end
      end
      
      xLims = xData;
      yLims = yData;
    end
    
    function [cellArray] = toStringCellArray(this, array)
      %%
      if (isnumeric(array))
        for i = 1 : length(array)
          cellArray{i} = num2str(array(i));
        end
      else
        cellArray = array;
      end
    end
    
    function [this] = plotLinesConnectingYData(this)
      %%
      [xTickLocations] = this.getXTickLocations();
      
      allYVals = cat(1, this.yData{:});
      lineHandles = plot(xTickLocations, allYVals);
      
      if (~isempty(this.legend))
        legend(lineHandles, this.legend);
      end
    end
    
    function [barHandle] = plotBar(this, x, ys, whiskerLength, barColor)
      %%
      xMin = x - (this.barWidth / 2);
      xMax = x + (this.barWidth / 2);
      verts = [ ...
        xMin this.barBaseline; ...
        xMin nanmean(ys); ...
        xMax nanmean(ys); ...
        xMax this.barBaseline; ...
        ];
      faces = [1 2 3 4];
      barHandle = patch('Faces', faces, 'Vertices', verts, 'FaceColor', barColor, 'EdgeColor', 'none');
      hold on;

      
      if (~isnan(whiskerLength))
        whiskerDisplacement = whiskerLength * sign(nanmean(ys));
        whiskerEndY = nanmean(ys) + whiskerDisplacement;
        
        whiskerWidth = this.whiskerWidth;
        if (isempty(whiskerWidth))
          whiskerWidth = this.barWidth * 0.25;
        end
        
        line([x x], [nanmean(ys), whiskerEndY], 'color', this.whiskerColor);
        line([x - whiskerWidth, x + whiskerWidth], [whiskerEndY whiskerEndY], 'color', this.whiskerColor)
      end
      
      if (this.plotDatapointsWithBars)
        plot(x, ys, '.k');
      end
      
    end
    
  end
  
end


classdef LineLabelling
  %LINELABELLING places labels and arrows on plotted lines
  
  properties
    drawConnectingLine = 0;
    horizontalAlignment = 'left';
    verticalAlignment = 'middle';
    labelLineProperties = struct('lineStyle', '-', 'lineWidth', 0.5, 'color', 'k');
    labelLineShrinkFraction = 0;
    markerProperties = struct('markerStyle', 'o', 'markerSize', 8, 'markerFaceColor', 'k', 'markerEdgeColor', 'k');
    
    offsetCurveArrowProperties = struct('lineStyle', '-', 'lineWidth', 2, 'color', ones(3, 1) * 0.3);
    
    arrowheadProperties = {}; % these get passed through to the arrow function
    labelTextProperties = {};
    preString = '';
  end
  
  methods (Static)
    
    function [] = test()
      figure('color', 'w');
      ts = 0:0.01:2*pi;
      xs = sin(ts);
      ys = cos(ts);
      handle = plot(xs, ys);
      lineLabelling = LineLabelling();
      lineLabelling.labelLineAtIndex(handle, 'top point', find(ys == max(ys)));
      
      lineLabelling.drawConnectingLine = 1;
      lineLabelling.labelLineAtIndex(handle, 'right point', find(xs == max(xs)), 'offsetOfLabel', [0.1, 0.1]);

      lineLabelling.horizontalAlignment = 'right';
      lineLabelling.labelLineShrinkFraction = 0.5;
      lineLabelling.labelLineAtIndex(handle, 'bottom point', find(ys == min(ys)), 'offsetOfLabel', [-0.1, -0.1]);
      
      lineLabelling.horizontalAlignment = 'right';
      lineLabelling.verticalAlignment = 'cap';
      lineLabelling.labelLineShrinkFraction = 0.3;
      lineLabelling.labelLineAtIndex(handle, 'left point', find(xs == min(xs)), 'offsetOfLabel', [-0.2, -0.1], ...
        'fractionAlongLabelEdgesToDrawLineTo', [-0.2, -0.2]);
      lineLabelling.plotMarkerAtIndex(handle, find(xs == min(xs)));
      
      lineLabelling.arrowheadProperties = {'length', 10};
      lineLabelling.drawOffsetCurveArrow(handle, 125:150, [-0.075 0]);
      
      %% test a broken case
      annotationsStringStyles = {'FontSize', 10, 'interpreter', 'latex'};
      lineLabeller = LineLabelling();
      lineLabeller.labelTextProperties = annotationsStringStyles;
      lineLabeller.horizontalAlignment = 'center';
      lineLabeller.verticalAlignment = 'top';
      lineLabeller.drawConnectingLine = 1;
      lineLabeller.labelLineShrinkFraction = 0;
      [textHandle] = lineLabeller.labelLineAtIndex(handle, '$\alpha = 0$, impulsive ', 100, 'offsetOfLabel', [-0.3,-0.05], ...
        'fractionAlongLabelEdgesToDrawLineTo', [0, 0]);
      thisBox = get(textHandle, 'extent');
      thisBox
      
      axis equal;
      axis off;
    end
  end
  
  methods
    
    function [this] = LineLabelling()
    end
    
    
    function [] = drawOffsetCurveArrow(this, plotLineHandle, indecesToDraw, offsetVector)
      %%
      xs = get(plotLineHandle, 'XData');
      ys = get(plotLineHandle, 'YData');
      
      xs = xs(indecesToDraw) + offsetVector(1);
      ys = ys(indecesToDraw) + offsetVector(2);
      
      
      lastTwoXs = xs(end-1:end);
      lastTwoYs = ys(end-1:end);
      
      if (isempty(this.arrowheadProperties))
        %         arrow([lastTwoXs(1) lastTwoYs(1)], [lastTwoXs(2) lastTwoYs(2)]);
        propertiesToUse = {'FaceColor', this.offsetCurveArrowProperties.color};
      else
        %         arrow([lastTwoXs(1) lastTwoYs(1)], [lastTwoXs(2) lastTwoYs(2)], this.arrowheadProperties{:});
        propertiesToUse = {'FaceColor', this.offsetCurveArrowProperties.color, this.arrowheadProperties{:}};
      end
      
      arrow([lastTwoXs(1) lastTwoYs(1)], [lastTwoXs(2) lastTwoYs(2)], propertiesToUse{:});

      
      hold on;
      plot(xs(1:(end-1)), ys(1:(end-1)), 'lineStyle', this.offsetCurveArrowProperties.lineStyle, ...
        'LineWidth', this.offsetCurveArrowProperties.lineWidth, ...
        'color', this.offsetCurveArrowProperties.color);
      
    end
    
    function [] = plotMarkerAtIndex(this, plotLineHandle, indexToLabelAt)
      %%
      xs = get(plotLineHandle, 'XData');
      ys = get(plotLineHandle, 'YData');
      
      thisX = xs(indexToLabelAt);
      thisY = ys(indexToLabelAt);
      
      hold on;
      plot(thisX, thisY, this.markerProperties.markerStyle, ...
        'MarkerSize', this.markerProperties.markerSize, ...
        'MarkerFaceColor', this.markerProperties.markerFaceColor, ...
        'MarkerEdgeColor', this.markerProperties.markerEdgeColor);
    end
    
    
    function [textHandle] = labelLineAtIndex(this, plotLineHandle, label, indexToLabelAt, varargin)
      %%
      offsetOfLabel = [0, 0];
      fractionAlongLabelEdgesToDrawLineTo = [0, 0];
      
      for i = 1 : 2 : length(varargin)
        option = varargin{i};
        value = varargin{i + 1};
        switch option
          case 'fractionAlongLabelEdgesToDrawLineTo'
            fractionAlongLabelEdgesToDrawLineTo = value;
          case 'offsetOfLabel'
            offsetOfLabel = value;
        end
      end
      
      xs = get(plotLineHandle, 'XData');
      ys = get(plotLineHandle, 'YData');
      
      thisX = xs(indexToLabelAt);
      thisY = ys(indexToLabelAt);
      
      if (isempty(this.labelTextProperties))
        textHandle = text(thisX + offsetOfLabel(1), thisY + offsetOfLabel(2), [this.preString label], ...
          'HorizontalAlignment', this.horizontalAlignment, 'VerticalAlignment', this.verticalAlignment);
      else
        textHandle = text(thisX + offsetOfLabel(1), thisY + offsetOfLabel(2), [this.preString label], ...
          'HorizontalAlignment', this.horizontalAlignment, 'VerticalAlignment', this.verticalAlignment, this.labelTextProperties{:});
      end
      
      if (this.drawConnectingLine)
        box = get(textHandle, 'Extent');
        %         lineOrigin = [0, 0];
        %         switch this.horizontalAlignment
        %           case 'left'
        %             lineOrigin(1) = box(1);
        %           case 'center'
        %             lineOrigin(1) = box(1) + box(3)/2;
        %           case 'right'
        %             lineOrigin(1) = box(1) + box(3);
        %         end
        %         switch this.verticalAlignment
        %           case 'bottom'
        %             lineOrigin(2) = box(2);
        %           case 'baseline'
        %             lineOrigin(2) = box(2);
        %           case 'middle'
        %             lineOrigin(2) = box(2) + box(4)/2;
        %           case 'top'
        %             lineOrigin(2) = box(2) + box(4);
        %           case 'cap'
        %             lineOrigin(2) = box(2) + box(4);
        %         end
        lineOrigin = get(textHandle, 'Position');
        lineOrigin(3) = [];
        
        lineOrigin(1) = lineOrigin(1) + fractionAlongLabelEdgesToDrawLineTo(1) * box(3);
        lineOrigin(2) = lineOrigin(2) + fractionAlongLabelEdgesToDrawLineTo(2) * box(4);
        
        lineXs = [thisX lineOrigin(1)];
        lineYs = [thisY lineOrigin(2)];
        
        midPointOfLine = [mean(lineXs), mean(lineYs)];
        
        lineXs = (lineXs - midPointOfLine(1)) * (1 - this.labelLineShrinkFraction);
        lineXs = lineXs + midPointOfLine(1);
        
        lineYs = (lineYs - midPointOfLine(2)) * (1 - this.labelLineShrinkFraction);
        lineYs = lineYs + midPointOfLine(2);
        
        line(lineXs, lineYs, ...
          'lineStyle', this.labelLineProperties.lineStyle, ...
          'LineWidth', this.labelLineProperties.lineWidth, ...
          'color', this.labelLineProperties.color);
      end
      
      
    end
    
  end
  
end


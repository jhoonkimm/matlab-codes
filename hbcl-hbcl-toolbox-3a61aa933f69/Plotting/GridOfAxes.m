classdef GridOfAxes
  
  properties
    
    numRows;
    numCols;
    
    verticalSpacing = 0.1;
    horizontalSpacing = 0.1;
    
    normalizedRectangleOfAxisGrid = [0 0 1 1];
    
    xSpacings = [];
    ySpacings = [];
    plotWidth = [];
    plotHeight = [];
    
    axesHandles = [];
    
  end
  
  methods (Static)
    
    function [] = test()
      %%
      
      figure('color', 'w');
      gridAx = GridOfAxes(2, 2);
      gridAx.verticalSpacing = 0.03;
      gridAx.horizontalSpacing = 0.03;
      
      gridAx.normalizedRectangleOfAxisGrid = [0.2 0 0.8 1];
      
      gridAx = gridAx.subplot(1,1);
      peaks;
      grid off; axis off; title('');
      gridAx = gridAx.subplot(2,1);
      peaks;
      grid off; axis off; title('');
      gridAx = gridAx.subplot(1,2);
      peaks;
      grid off; axis off; title('');
      gridAx = gridAx.subplot(2,2);
      peaks;
      grid off; axis off; title('');
    end
    
  end
  
  methods
    
    function [this] = GridOfAxes(numRows, numCols)
      %%
      this.numRows = numRows;
      this.numCols = numCols;
      this.axesHandles = ones(numRows, numCols) * NaN;
    end
    
    function [this, handle] = subplotFromIndex(this, index)
      %%
%       [row, col] = ind2sub([this.numRows, this.numCols], index);
      [col, row] = ind2sub([this.numCols, this.numRows], index);
      [this, handle] = subplot(this, row, col);
    end
    
    function [this, handle] = subplot(this, row, col)
      %%
      this = updateAxesSpacing(this);
      if (~isnan(this.axesHandles(row, col)))
        axes(this.axesHandles(row, col));
      else
      this.axesHandles(row, col) = axes('position', ...
        [this.xSpacings(col), this.ySpacings(this.numRows - row + 1), ...
        this.plotWidth, this.plotHeight]);
      end
      handle = this.axesHandles(row, col);
    end
    
  end
  
  methods (Access = private)
    
    function [this] = updateAxesSpacing(this)
      %%
      if (isempty(this.xSpacings))
        [this.xSpacings, this.plotWidth] = this.updateAxisSpacing( ...
          this.normalizedRectangleOfAxisGrid(1), this.normalizedRectangleOfAxisGrid(1) + this.normalizedRectangleOfAxisGrid(3), ...
          this.horizontalSpacing, this.numCols);
        
        [this.ySpacings, this.plotHeight] = this.updateAxisSpacing( ...
          this.normalizedRectangleOfAxisGrid(2), this.normalizedRectangleOfAxisGrid(2) + this.normalizedRectangleOfAxisGrid(4), ...
          this.verticalSpacing, this.numRows);
      end
    end
    
    function [spacings, lengthOfPlot] = updateAxisSpacing(this, start, stop, spacing, numPlotsInDimension)
      %%
      lengthOfPlot = ((stop - start) - (numPlotsInDimension + 1) * spacing) / numPlotsInDimension;
      
      spacings = linspace( ...
        start + spacing, ...
        stop - spacing - lengthOfPlot, ...
        numPlotsInDimension);
    end
    
  end
  
end


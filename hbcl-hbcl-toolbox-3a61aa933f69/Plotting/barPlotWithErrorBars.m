function handles = barPlotWithErrorBars(x,y,errors,varargin)

% BARPLOTWITHERRORBARS
%  Similar to BARWEB, but less obnoxious about input parameters
%  Inputs:
%   x: x-axis values; if you don't care, then input (1:length(y))
%   y: bar values; input column vector, add columns to create groups of bars
%   errors: standard dev/error bars to plot; erros should be same size as y
%
%  Examples:
%   x=[1 2 4 8]'; y=[1 2; 2 3; 3 4; 5 6]; errors=ones(4,2);
%   h = barPlotWithErrorBars(x,y,errors);
%   h = barPlotWithErrorBars(x,y,errors,barLineWidth,1);

%  Karl Zelik
%  11/8/09


%% Default variables
figureNum = [];
titleName = '';
errorBarWidth = 80;
errorBarLineWidth = 2;
barWidth = 1;
barLineWidth = 2;
barLineColor = 'k';
yAxisLabel = '';
xAxisLabel = '';
textFontSize = 14;
axesFontSize = 14;
errorBarType = ''; % eventually build in ability to input different upper and lower error bars

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
    case 'barLineWidth'
      barLineWidth = val;
    case 'barWidth'
      barWidth = val;  
    case 'barLineColor'
      barLineColor = val; 
    case 'xAxisLabel'
      xAxisLabel = val;   
    case 'yAxisLabel'
        yAxisLabel = val;
    case 'title'
        titleName = val;
    case 'errorBarWidth'
        errorBarWidth = val;
    case 'errorBarLineWidth'
        errorBarLineWidth = val;
    case 'errorBarType'
        errorBarType = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


%% Bar function does not like single inputs
if size(x,1)==1 & size(x,2)==1
    x = [0; x];
    y = [zeros(1,size(y,2)); y];
    e = [zeros(1,size(e,2)); e];
end


%% Create figure
if isempty(figureNum)
    figure(); 
else figure(figureNum); 
end
set(gcf,'DefaultTextFontSize',textFontSize);
set(gcf,'DefaultAxesFontSize',axesFontSize);


%% Bar plot
hold on;
handles.bars = bar(x,y,barWidth,'edgecolor',barLineColor, 'linewidth', barLineWidth);
barsPerGroup = size(y,2);

%% Add error bars
for j=1:size(y,1)
    for i=1:barsPerGroup
        if barsPerGroup <= 5
            
            % Define x-offset value to plot error bars in proper locations
            % I approximated these values b/c not sure how/where to find precise bar width values
            if barsPerGroup==1, xoffset = 0;
            elseif barsPerGroup==2, halfBarWidth = 0.022; xoffset = [-halfBarWidth halfBarWidth];
            elseif barsPerGroup==3, halfBarWidth = 0.019; xoffset = [-2*halfBarWidth 0 2*halfBarWidth];
            elseif barsPerGroup==4, halfBarWidth = 0.0135; xoffset = [-3*halfBarWidth -halfBarWidth halfBarWidth 3*halfBarWidth];
            elseif barsPerGroup==5, halfBarWidth = 0.0115; xoffset = [-4*halfBarWidth -2*halfBarWidth 0 2*halfBarWidth 4*halfBarWidth];
            end
            
            xerrors(j,i) = x(j) + xoffset(i);
            hold on; handles.errorbars(j,i) = errorbar(xerrors(j,i), y(j,i), errors(j,i), 'k', 'linestyle', 'none', 'linewidth', errorBarLineWidth);           
            hold on; errorbar_tick(handles.errorbars(j,i),errorBarWidth);
            
        else error('Need to generalize code for more than 5 inputs. X-offset for errorbars is unknown.');
        end
    end
end


%% Label axes, title
xlabel(xAxisLabel); ylabel(yAxisLabel); title(titleName);





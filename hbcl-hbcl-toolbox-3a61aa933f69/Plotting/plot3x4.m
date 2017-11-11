function plot3x4(AnkleAngles,AnkleMoment,AnklePower,KneeAngles,KneeMoment,KneePower,HipAngles,HipMoment,HipPower,TrunkAngles,TrunkMoment,TrunkPower,varargin)

% PLOT3X4 plot angles, moments and powers for ankle, knee, hip & trunk
%  Plot standard sagittal plane joint "3x4": angles, moments and powers 
%  for ankle, knee hip & trunk. 
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
%   plot3x3(avgRAnkleAnglesR,avgRAnkleMomentR,avgRAnklePowerR,avgRKneeAnglesR,avgRKneeMomentR,avgRKneePowerR,avgRHipAnglesR,avgRHipMomentR,avgRHipPowerR,avgTrunkAnglesR,avgTrunkMomentR,avgTrunkPowerR)

%  Karl Zelik
%  11/19/09


%% Default plot options
textFontSize = 14;
axesFontSize = 14;
lineWidth = 3;
lineStyle = '-';
figureNum = [];
axesColorOrder = [];
xTickLabel = {'0','25','50','75','100'}; % default plot is heelstrike to heelstrike
xAxisLabel = 'Gait Cycle (%)';
angleUnits = 'degrees';
momentUnits = 'Nm/kg';
powerUnits = 'W/kg';
labelFlexAndExt = 1;
plotZeroLine = 1;
nonDimConstant = []; % input in order angle, moment, power

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
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


%% Create figure
if isempty(figureNum)
    figure(); 
else figure(figureNum); 
end
set(gcf,'DefaultTextFontSize',textFontSize);
set(gcf,'DefaultAxesFontSize',axesFontSize);
if isempty(axesColorOrder)
    axesColorOrder = get(0,'DefaultAxesColorOrder'); 
else
    set(gca,'ColorOrder',axesColorOrder);
end
xTick = linspace(1,size(AnkleAngles,1),size(xTickLabel,2));


%% Store all joint data in 1 array
data = [AnkleAngles KneeAngles HipAngles TrunkAngles AnkleMoment KneeMoment HipMoment TrunkMoment AnklePower KneePower HipPower TrunkPower];
numConditions = size(AnkleAngles,2);


%% Plotting loop
for i=1:numConditions
    lineColor = axesColorOrder(i,:); % set line color - different for each condition
    
    for j=1:12
        eval(['h' num2str(j) '= subplot(3,4,j);']);
        set(gca,'Box','on');
        
        % Default Labeling of Axes
        % Column titles
        if j==1, title('Ankle'); end
        if j==2, title('Knee'); end
        if j==3, title('Hip'); end
        if j==4, title('Trunk'); end
        
        % y-axis labels
        if j==1, ylabel(strcat('Angle (',angleUnits,')'));
        elseif j==5, ylabel(strcat('Moment (',momentUnits,')'));
        elseif j==9, ylabel(strcat('Power (',powerUnits,')')); 
        end
        
        % x-axis label
        if j==10, xlabel(xAxisLabel); end
       
        % x-axis ticks
        set(gca,'XTick',xTick); 
        if j==9 | j==10 | j==11 | j==12
            set(gca,'XTickLabel',xTickLabel); 
        else
            set(gca,'XTickLabel',{''});
        end
        
        % Plot data
        hold on; plot(data(:,j*numConditions-(numConditions-i)),'LineWidth',lineWidth,'Color',lineColor,'LineStyle',lineStyle);
        setYAxisLimits(gca)
    end
end
   

%% Set equivalent axes
setEqualAxes([h1 h2 h3 h4]); % angles
setEqualAxes([h5 h6 h7 h8]); % moments
setEqualAxes([h9 h10 h11 h12]); % powers

for j = [2 3 4 6 7 8 10 11 12]
   subplot(3,4,j); set(gca,'YTickLabel',{''});
end

%% Plot zero line
if plotZeroLine
    for j=1:12
        subplot(3,4,j); 
        ax=axis; 
        hold on; plot([ax(1) ax(2)],[0 0],'k','LineWidth',1); axis(ax);
    end
end
    
%% Add Flexion & Extension text
% Angles
if ~isempty(nonDimConstant)
    subplot(341); hold on;
else
    subplot(344); hold on;
end
ax = axis;
ys = linspace(ax(3),ax(4),9);
if labelFlexAndExt
    text(ax(2)*1.05,ys(2),'Flex','FontSize',10);
    text(ax(2)*1.05,ys(8),'Ext','FontSize',10);
end

% Moments
if ~isempty(nonDimConstant)
    subplot(345); hold on;
else
    subplot(348); hold on;
end
ax = axis;
ys = linspace(ax(3),ax(4),9);
if labelFlexAndExt
    text(ax(2)*1.05,ys(2),'Flex','FontSize',10);
    text(ax(2)*1.05,ys(8),'Ext','FontSize',10);
end


%% Add dimensionless axes
if ~isempty(nonDimConstant)
    subplot(344); hold on; AX3 = addDimensionlessYAxis(gca,nonDimConstant(1)); % angle
    set(AX3(2),'XTickLabel',{''});
    subplot(348); hold on; AX6 = addDimensionlessYAxis(gca,nonDimConstant(2)); % moment
    set(AX6(2),'XTickLabel',{''});
    subplot(3,4,12); hold on; AX9 = addDimensionlessYAxis(gca,nonDimConstant(3)); % power
    set(AX9(2),'XTickLabel',{''}); set(gca,'XTickLabel',xTickLabel); 
end


end


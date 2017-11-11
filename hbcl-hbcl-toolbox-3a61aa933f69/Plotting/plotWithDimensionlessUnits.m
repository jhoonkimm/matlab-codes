function plotWithDimensionlessUnits(X,Y,nonDimConstant,varargin)

% PLOTWITHDIMENSIONLESSUNITS plot 2nd y-axis
%  Plot equivalent dimensiolness units on right side y-axis.
%  Input X,Y data the same way you would for PLOT.
%  Input nondimensionalization constant (nonDimConstant) where
%   Ydimensionless = Ydimensional*nonDimConstant
%
%  Examples:
%   X=0:10:360; Y=sind(X); nonDimConstant=0.1;
%   plotWithDimensionlessUnits(X,Y,nonDimConstant);

%  Karl Zelik
%  11/21/09


%% Attempt 3
% Default variables
axesColorOrder = get(0,'DefaultAxesColorOrder');
lineWidth = 2;
lineColor = 'b';
lineStyle = '-';
marker = 'none';
markerSize = 10;
markerEdgeColor = 'auto';
markerFaceColor = 'none';
addAxisOnly = 0;

%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  
  switch opt 
    case 'LineWidth'
      lineWidth = val;
    case 'Color'
      lineColor = val;
    case 'LineStyle'
      lineStyle = val;
    case 'Marker'
      marker = val;
    case 'MarkerSize'
      markerSize = val;
    case 'MarkerEdgeColor'
      markerEdgeColor = val;
    case 'MarkerFaceColor'
      markerFaceColor = val;
    case 'axesColorOrder'
      axesColorOrder= val;
    case 'addAxisOnly'
      addAxisOnly= val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


if size(Y,2)==1
    hold on; h = plot(X,Y,'LineWidth',lineWidth,'Color',lineColor,'LineStyle',lineStyle,'Marker',marker,'MarkerSize',markerSize,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor);
else
    hold on; set(gcf,'DefaultAxesColorOrder',axesColorOrder);
    h = plot(X,Y,'LineWidth',lineWidth,'LineStyle',lineStyle,'Marker',marker,'MarkerSize',markerSize,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor);
end


addDimensionlessYAxis(gca,nonDimConstant);








%% Attempt 2
% % %% Default variables
% % axesColorOrder = get(0,'DefaultAxesColorOrder');
% % lineColor = 'b';
% % lineWidth = 2;
% % 
% % 
% % %% Optional input arguments
% % opt_argin = varargin;
% % while length(opt_argin) >= 2,
% %   opt = opt_argin{1};
% %   val = opt_argin{2};
% %   opt_argin = opt_argin(3:end);
% %   
% %   switch opt 
% %     case 'LineWidth'
% %       lineWidth = val;
% %     case 'Color'
% %       lineColor= val;
% %     case 'axesColorOrder'
% %       axesColorOrder= val;
% %     otherwise
% %       error('\nError: incorrect varargin\n')
% %   end
% % end
% % 
% % %% Plot x & y
% % hold on; 
% % plot(X,Y) %,'Color',lineColor,'LineWidth',lineWidth);
% % ax1 = axis; ymin1=ax1(3); ymax1=ax1(4);
% % 
% % %% Add second (dimensionless) axis and scale new axis to align properly with axis 1
% % addaxisreset(gca); % remove old right y-axis if it exists to prepare for new right y-ais
% % hold on; addaxis(X,Y*nonDimConstant,'k','linewidth',0.0002)
% % ymin2=ymin1*nonDimConstant; ymax2=ymax1*nonDimConstant;
% % addaxisset([ymin2 ymax2],2);
% % 
% % %% Replot x % y to cover up plot from addaxis
% % set(gca,'ColorOrder',axesColorOrder)
% % hold on; plot(X,Y) %,'Color',lineColor,'LineWidth',lineWidth);
% % 
% % 
% % 




% % %%
% % %% Attempt 1
% % % Default variables
% % axesColorOrder = get(0,'DefaultAxesColorOrder');
% % lineWidth = 2;
% % lineColor = 'b';
% % lineStyle = '-';
% % marker = 'none';
% % markerSize = 10;
% % markerEdgeColor = 'auto';
% % markerFaceColor = 'none';
% % addAxisOnly = 0;
% % 
% % %% Optional input arguments
% % opt_argin = varargin;
% % while length(opt_argin) >= 2,
% %   opt = opt_argin{1};
% %   val = opt_argin{2};
% %   opt_argin = opt_argin(3:end);
% %   
% %   switch opt 
% %     case 'LineWidth'
% %       lineWidth = val;
% %     case 'Color'
% %       lineColor = val;
% %     case 'LineStyle'
% %       lineStyle = val;
% %     case 'Marker'
% %       marker = val;
% %     case 'MarkerSize'
% %       markerSize = val;
% %     case 'MarkerEdgeColor'
% %       markerEdgeColor = val;
% %     case 'MarkerFaceColor'
% %       markerFaceColor = val;
% %     case 'axesColorOrder'
% %       axesColorOrder= val;
% %     case 'addAxisOnly'
% %       addAxisOnly= val;
% %     otherwise
% %       error('\nError: incorrect varargin\n')
% %   end
% % end
% % 
% % [AX,H1,H2] = plotyy(X,Y,X,Y*nonDimConstant);
% % delete(H1);
% % delete(H2);
% % 
% % set(AX(1),'YColor','k')
% % set(AX(2),'YColor','k')
% % 
% % ax1 = axis(AX(1)); 
% % ax2 = [ax1(1) ax1(2) ax1(3)*nonDimConstant ax1(4)*nonDimConstant];
% % 
% % axes(AX(2)); % set to current axis
% % set(AX(2), 'YTickMode','auto')
% % axis(ax2);
% % axes(AX(1));
% % 
% % if ~addAxisOnly
% %     if size(Y,2)==1
% %         hold on; h = plot(X,Y,'LineWidth',lineWidth,'Color',lineColor,'LineStyle',lineStyle,'Marker',marker,'MarkerSize',markerSize,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor);
% %     else
% %         hold on; set(gcf,'DefaultAxesColorOrder',axesColorOrder);
% %         h = plot(X,Y,'LineWidth',lineWidth,'LineStyle',lineStyle,'Marker',marker,'MarkerSize',markerSize,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor);
% %     end
% % end





function PlotAxisAtOrigin(ax);
%PlotAxisAtOrigin Plot 2D axes through the origin
%   This is a 2D version of Plot3AxisAtOrigin written by Michael Robbins
%   File exchange ID: 3245. 
%
%   Have hun! 
%
%   Example:
%   x = -2*pi:pi/10:2*pi;
%   y = sin(x);
%   PlotAxisAtOrigin(x,y)
%

% PLOT
% if nargin == 2 
%     plot(x,y);
%     hold on;
% else
%     display('   Not 2D Data set !')
% end;

% GET TICKS
X=get(ax,'Xtick');
Y=get(ax,'Ytick');

% GET LABELS
XL=get(ax,'XtickLabel');
YL=get(ax,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(ax,'XLim'))./100;
Yoff=diff(get(ax,'YLim'))./100;

% DRAW AXIS LINEs
plot(get(ax,'XLim'),[0 0],'k');
plot([0 0],get(ax,'YLim'),'k');

% Plot new ticks  
for i=1:length(X)
    plot([X(i) X(i)],[0 Yoff],'-k');
end;
for i=1:length(Y)
   plot([Xoff, 0],[Y(i) Y(i)],'-k');
end;

% % ADD LABELS
% text(X,zeros(size(X))-2.*Yoff,XL);
%  text(zeros(size(Y))-3.*Xoff,Y,YL);

box off;
% axis square;
axis off;
set(gcf,'color','w');


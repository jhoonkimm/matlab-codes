function setYAxisLimits(ax,varargin)

% SETYAXISLIMITS set y-axis limits
%  Similar to command AXIS, but allows user to specify amount of vertical white space
%  above and below the data. The command AXIS TIGHT is usually too zoomed in on data 
%  and command AXIS NORMAL is too zoomed out. SETYAXISLIMITS allows user to define 
%  the vertical excursion.
%
%  SETYAXISLIMITS(AX) expands y-limits by 5% over 'axis tight' for axis AX
%
%  SETYAXISLIMITS(AX,VAL) expands y-limits by VAL*100% over 'axis tight' for axis AX
%
%  Examples:
%    x=1:10; y=1:10; figure; plot(x,y); setYAxisLimits(gca,0.1)

% Karl Zelik
% updated 10/29/09

if nargin==2 & ~isempty(varargin{1})
    scaleFactor = varargin{1};
else
    scaleFactor = 0.1; % default
end


axes(ax); % make this axis current
oldxlimits = xlim(ax);

axis tight % sets the axis limits to the range of the data
oldylimits = ylim(ax);
yrange = range(oldylimits);
newylimits = [oldylimits(1)-yrange*scaleFactor/2 oldylimits(2)+yrange*scaleFactor/2];
% example: for scaleFactor = 0.05 (5%), add 2.5% at top ad 2.5% at bottom


axis([oldxlimits newylimits])
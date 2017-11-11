function setXAxisLimits(ax,varargin)

% SETXAXISLIMITS set x-axis limits
%  Similar to command AXIS, but allows user to specify amount of horizontal white space
%  above and below the data. The command AXIS TIGHT is usually too zoomed in on data 
%  and command AXIS NORMAL is too zoomed out. SETXAXISLIMITS allows user to define 
%  the horizontal excursion.
%
%  SETYAXISLIMITS(AX) expands x-limits by 5% over 'axis tight' for axis AX
%
%  SETYAXISLIMITS(AX,VAL) expands x-limits by VAL*100% over 'axis tight' for axis AX
%
%  Examples:
%    x=1:10; y=1:10; figure; plot(x,y); setXAxisLimits(gca,0.1)

% Karl Zelik
% updated 11/19/09

if nargin==2 & ~isempty(varargin{1})
    scaleFactor = varargin{1};
else
    scaleFactor = 0.1; % default
end

axes(ax); % make this axis current
oldylimits = ylim(ax);

axis tight % sets the axis limits to the range of the data
oldxlimits = xlim(ax);
xrange = range(oldxlimits);
newxlimits = [oldxlimits(1)-xrange*scaleFactor/2 oldxlimits(2)+xrange*scaleFactor/2];
% example: for scaleFactor = 0.05 (5%), add 2.5% at top ad 2.5% at bottom


axis([newxlimits oldylimits])
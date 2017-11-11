function [AX] = addDimensionlessYAxis(ax,nonDimConstant,varargin)

% ADDDIMENSIONLESSYAIS add 2nd y-axis with dimensionless units
%  Add equivalent dimensiolness units on right side y-axis.
%  Input nondimensionalization constant (nonDimConstant) where
%   Ydimensionless = Ydimensional*nonDimConstant
%
%  CAUTION: Axis are not permantently linked, so if you
%   rescale 1st xis you need to rerun this
%   ADDDIMENSIONLESSYAXIS script.
%
%  Examples:
%   X=0:10:360; Y=sind(X); nonDimConstant=0.1;
%   figure; plot(X,Y); addDimensionlessYAxis(gca,nonDimConstant);

%  Karl Zelik
%  11/21/09


% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'yAxisScaleFactor'
      yAxisScaleFactor = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end

axes(ax); % set input axis as current
hold on; oldAX = axis;
[AX,H1,H2] = plotyy(0,0,0,0);
axis(oldAX);
delete(H1); % delete plots, but keep both sets of axes
delete(H2);

% Set color of axes
set(AX(1),'YColor','k')
set(AX(2),'YColor','k')

% Calculate equivalent dimensionless axis limits
ax1 = axis(AX(1)); 
ax2 = [ax1(1) ax1(2) ax1(3)*nonDimConstant ax1(4)*nonDimConstant];

axes(AX(2)); % set dimensionless axis to current axis
set(AX(2), 'YTickMode','auto')
axis(ax2); % rescale dimensionless axis to correspond to dimensional axis limits
axes(AX(1)); % set dimensional axis to current axis



function animationexample(showorsave)
% Demonstrates Matlab animation, with option to save to avi file.
% Use animationexample(1) to save the file.

if nargin == 0    % no input argument given
  showorsave = 0;
end

% generate some sample data
l1 = 0.4; l2 = 0.4;
w = 2*pi;
ts = 0:1/15:2;
theta1s = 10*pi/180*cos(w*ts);
theta2s = 20*pi/180*sin(w*ts);
x1 = -l1*sin(theta1s(1)); y1 = l1*cos(theta2s(1));
x2 = x1 - l2*sin(theta2s(1)); y2 = y1 + l2*cos(theta2s(1));

% Draw the arm as a single line with joints
clf;
hArm = line([0 x1 x2],[0 y1 y2],'LineWidth', 2, 'LineStyle', '-', ...
  'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
  'MarkerEdgeColor', 'r');
axis([-0.5 0.5 0 1]);

Mov(1:length(ts)) = getframe; % This pre-allocates memory for a movie

for i=1:length(ts)
  x1 = -l1*sin(theta1s(i)); y1 = l1*cos(theta2s(i));
  x2 = x1 - l2*sin(theta2s(i)); y2 = y1 + l2*cos(theta2s(i));
  xinfo = [0 x1 x2]; yinfo = [0 y1 y2];
  set(hArm,'Xdata', xinfo, 'Ydata', yinfo);
  drawnow; 
  if showorsave == 0 % just show the movie
    pause(1/15); % wait before continuing
  else               % save the movie
    Mov(i) = getframe;
  end % if
end % for

if showorsave == 1 % save the movie to a file
  movie2avi(Mov, 'twolinkAnimation.avi');
end

end % animationexample
function out = circle(x,y,r)
% CIRCLE(x,y,r)  draws a circle centered at (x,y) with radius r
% 'Circle' is a constrained call to the 'Rectangle' function.
% I think it's silly that MATLAB does not have
% a 'Circle' function already.

out = rectangle('Curvature',[1 1],'Position',[x-r y-r 2*r 2*r]);

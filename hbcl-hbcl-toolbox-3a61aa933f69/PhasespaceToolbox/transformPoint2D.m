function [p] = transformPoint2D(p, xyTheta)
theta = xyTheta(3);
p = [cos(theta) -sin(theta);
 sin(theta) cos(theta)] * p;
p = p + xyTheta(1:2);
end

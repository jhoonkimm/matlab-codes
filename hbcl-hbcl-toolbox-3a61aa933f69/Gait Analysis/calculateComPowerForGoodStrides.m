function [RComP LComP RComPAverage LComPAverage RComPStd LComPStd velocity velocityAverage] = ...
  calculateComPowerForGoodStrides(rightForces, leftForces, t, goodStrides, treadmillSpeed, varargin)
%% calculateComPowerForGoodStrides given left/right leg forces, COM velocity, and indeces of good strides
%
% rightForces and leftForces must be structs with fields: {'x', 'y', 'z'}
% each of dimension: [n, numberOfStrides]
%
% t is the stride time, dimension [numberOfStrides]
%
% goodStrides lists the indeces of good strides


Rx = rightForces.x;
Ry = rightForces.y;
Rz = rightForces.z;

Lx = leftForces.x;
Ly = leftForces.y;
Lz = leftForces.z;

RComP = zeros(size(Rx,1), length(goodStrides));
LComP = zeros(size(Lx,1), length(goodStrides));
velocity = zeros(size(Lx,1), 3, length(goodStrides));

for j = 1 : size(goodStrides, 2)
  RGRF = [Rx(:,goodStrides(j)) Ry(:,goodStrides(j)) Rz(:,goodStrides(j))];
  LGRF = [Lx(:,goodStrides(j)) Ly(:,goodStrides(j)) Lz(:,goodStrides(j))];
  [RComP(:,j) ,LComP(:,j), velocity(:, :, j)] = calculateComWorkRateForSteadyStateGait( ...
    RGRF, LGRF, treadmillSpeed, t(goodStrides(j)), 'subtractVelocitySlope', 1);
end

RComPAverage = mean(RComP, 2);
RComPStd = std(RComP, 0, 2);

LComPAverage = mean(LComP, 2);
LComPStd = std(LComP, 0, 2);

velocityAverage = mean(velocity, 3);

end
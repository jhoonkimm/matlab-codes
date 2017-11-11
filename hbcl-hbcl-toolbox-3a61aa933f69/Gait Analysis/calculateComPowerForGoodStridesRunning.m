function [comPower comPowerAverage comPowerStandardDeviation velocity velocityAverage] = ...
  calculateComPowerForGoodStridesRunning(forces, t, goodStrides, treadmillSpeed, varargin)
%% calculateComPowerForGoodStrides given left/right leg forces, COM velocity, and indeces of good strides
%
% rightForces and leftForces must be structs with fields: {'x', 'y', 'z'}
% each of dimension: [n, numberOfStrides]
%
% t is the stride time, dimension [numberOfStrides]
%
% goodStrides lists the indeces of good strides


comPower = zeros(size(forces.x,1), length(goodStrides));
velocity = zeros(size(forces.x,1), 3, length(goodStrides));

for j = 1 : size(goodStrides, 2)
  GRF = [forces.x(:, goodStrides(j)) forces.y(:, goodStrides(j)) forces.z(:, goodStrides(j))];
  [comPower(:,j), velocity(:, :, j)] = calculateComWorkRateForSteadyStateGaitRunning( ...
    GRF, treadmillSpeed, t(goodStrides(j)), 'subtractVelocitySlope', 1);
end

comPowerAverage = mean(comPower, 2);
comPowerStandardDeviation = std(comPower, 0, 2);
velocityAverage = mean(velocity, 3);

end
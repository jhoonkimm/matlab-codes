

%%
% A test of using the RigidBody class to track markers on a rigid body.
% There are a few use cases that we want to handle with this code.
%
% - First, there are times when we have a rigid body with more than 3
% markers, but during some frames some markers are missing in such a way
% that we still have 3 or more markers. We want to fill in the missing
% markers with using the existing markers on the body. (currently we have
% this capability).
%
% - Second, there are times when we have a rigid body with 3 or more
% markers, and there are periods where some markers are missing such that
% we have less than 3 markers remaining. If this period is short enough, we
% want to fill in the missing markers by using an interpolation of the
% rigid body data across the gap, comined with any existing information
% during the gap. This might require assumptions on how fast the body is
% expected to move (or might not?)
%

%
% how do we do 2? Eventually a function would exist in RigidBody, e.g.:
% fillMarkerDataAcrossGap(measuredMarkerData)
%
% There are sort of two cases. If the gap contains no marker data, then the
% best we might be able to do is spline the rigid body data across the gap,
% then recreate the markers from that. If there are markers in a frame in
% the gap, the best we can do at that point is some compromise between the
% rigid body frame before the point, the limited data we have in this
% frame, and the rigid body frame after the point. This is a good problem
% to solve with a kalman filter type solution, where:
%
% estimated state is the rigid body state in each frame where we don't have a
% rigid body state yet.
%
% Measurements are:
% - the locations of the trackers in each of the unknown frames.
% - an estimate of the state in each frame based on the average of the one
%   estimated state before and after. (note that this won't spline!)
%
% The simplest thing at the beginning is probably just to be lazy and
% estimate every state in the trajectory, whether we have 3 markers in the
% frame or not.
%
% The two major parameters here will be
% - ratio of trust in markers to trust in smooth rigid body states
% - how fast to change states in the filter iteration process with respect
%   to the measurement error.
% there would also be a few others to normalize orientations and positions
% with respect to each other.
%
% If we want to change this so that it would spline rigid body states, we
% could add some estimate of dynamics to the estimated state. So you would
% estimate the rigid body state, as well as the velocity and angular
% velocity of the rigid body at each tick. This would add a few more
% normalization parameters.
%
% At the moment, there might be a lazier approach, which is just to use
% fminunc again, with all of the states thrown in, and estimation error set
% to the measurement error.
%
%

clear all;
close all;
format compact;

r = RigidBody;

r = r.addTrackerPositionInBodyFrame([0 0 0]');
r = r.addTrackerPositionInBodyFrame([1 0 0]');
r = r.addTrackerPositionInBodyFrame([0 2 0]');
r = r.addTrackerPositionInBodyFrame([0 0 3]');

r = r.addTrackerPositionInBodyFrame([3 0 3]');

% for i = 1:10
%   r = r.addTrackerPositionInBodyFrame([10 10 10]');
%   r = r.addTrackerPositionInBodyFrame([1 2 3]');
%   r = r.addTrackerPositionInBodyFrame([4 -6 5]');
% end

numTimePoints = 20;
markerNoiseStandardDeviation = 0.3; %05; %

%%
linearArray = linspace(1, 2, numTimePoints);
distanceMultiple = 3;
positionTrajectory = [linearArray*10; sin(linearArray*5); linearArray*2]*distanceMultiple;
orientationTrajectory = [linearArray*pi; linearArray*pi*0.5 * 0; linearArray*pi/16*0];

pointsInWorld = cell(numTimePoints, 1);
noisyPointsInWorld = cell(numTimePoints, 1);
noisyPointsInWorldMissingMarker = cell(numTimePoints, 1);
for i = 1:numTimePoints
  r.position = positionTrajectory(:,i);
  r.orientation = orientationTrajectory(:,i);
  
  pointsInWorld{i} = r.getTrackerPositionsInWorld();
  noisyPointsInWorld{i} = pointsInWorld{i} + ...
    randn(size(pointsInWorld{i})) .* markerNoiseStandardDeviation;
  
  noisyPointsInWorldMissingMarker{i} = noisyPointsInWorld{i};
%   noisyPointsInWorldMissingMarker{i}(:, 1) = NaN;
end

%%
figure();
subplot(2,2,[1]);
for i = 1:numTimePoints
  r.plotTheseWorldPoints(pointsInWorld{i}, {'ko'});
  if (0)
    r.plotTheseWorldPoints(noisyPointsInWorld{i}, {'go'});
    r.plotTheseWorldPoints(noisyPointsInWorldMissingMarker{i}, {'bo'});
  end
end
axis equal;

subplot(2, 2, 3:4);
plot(positionTrajectory', 'k');
hold on;

%%
noisyPositionTrajectory = positionTrajectory;
noisyMissingPositionTrajectory = positionTrajectory;
noisyMissingPositionTrajectorySVD = positionTrajectory;

noisyOrientationTrajectory = orientationTrajectory;
noisyMissingOrientationTrajectory = orientationTrajectory;
noisyMissingOrientationTrajectorySVD = orientationTrajectory;

r.position = [0 0 0]';
r.orientation = [0 0 0]';
totalTime = 0;
totalTimeSVDMethod = 0;

for i = 1:numTimePoints
  actualPosition = positionTrajectory(:,i);
  actualOrientation = orientationTrajectory(:,i);
  
  %   [position, orientation, ...
  %     trackersFromTemplate, templateFitError] = ...
  %     r.calculatePositionAndOrientationFromTrackersInWorld(pointsInWorld{i});
  
  [position, orientation, ...
    trackersFromTemplate, templateFitError] = ...
    r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(pointsInWorld{i});
  r.position = position;
  r.orientation = orientation;
  
  if (templateFitError > 1e-3)
    error('recreated markers look good, but error is bad!');
  end
  
  %   if ((norm(position - actualPosition) + norm(orientation - actualOrientation)) > 1e-2)
  %     norm(position - actualPosition)
  %     norm(orientation - actualOrientation)
  %     error('we were unable to recreate the position and orientation exactly from perfect data!');
  %   end
  %   if (templateFitError > 1e-3)
  %     error('recreated markers look good, but error is bad!');
  %   end
  
  [positionFromNoisyData, orientationFromNoisyData, ...
    noisyTrackersFromTemplate, templateFitErrorToNoisyData] = ...
    r.calculatePositionAndOrientationFromTrackersInWorld(noisyPointsInWorld{i});
  %   if ((norm(position - actualPosition) + norm(orientation - actualOrientation)) > ...
  %       markerNoiseStandardDeviation * 8)
  %     error(['recreation of marker data is more than some number of standard deviations ' ...
  %       'from correct... might be bad?']);
  %   end
  
  noisyPositionTrajectory(:,i) = positionFromNoisyData;
  noisyOrientationTrajectory(:,i) = orientationFromNoisyData;
  
  tic
  [positionFromNoisyMissingData, orientationFromNoisyMissingData, ...
    noisyMissingTrackersFromTemplate, templateFitErrorToNoisyMissingData] = ...
    r.calculatePositionAndOrientationFromTrackersInWorld(noisyPointsInWorldMissingMarker{i});
  totalTime = totalTime + toc;
  
  %   if ((norm(position - actualPosition) + norm(orientation - actualOrientation)) > ...
  %       markerNoiseStandardDeviation * 8)
  %     error(['recreation of marker data is more than some number of standard deviations ' ...
  %       'from correct... might be bad?']);
  %   end
  
  noisyMissingPositionTrajectory(:,i) = positionFromNoisyMissingData;
  noisyMissingOrientationTrajectory(:,i) = orientationFromNoisyMissingData;
  
  tic
  [positionFromNoisyMissingDataSVD, orientationFromNoisyMissingDataSVD, ...
    noisyMissingTrackersFromTemplateSVD, templateFitErrorToNoisyMissingDataSVD] = ...
    r.calculatePositionAndOrientationFromTrackersInWorldLeastSquares(noisyPointsInWorldMissingMarker{i});
  totalTimeSVDMethod = totalTimeSVDMethod + toc;
  
  noisyMissingPositionTrajectorySVD(:,i) = positionFromNoisyMissingDataSVD;
  noisyMissingOrientationTrajectorySVD(:,i) = orientationFromNoisyMissingDataSVD;
  
  
  subplot(2, 2, 2);
  r.plotTheseWorldPoints(pointsInWorld{i}, {'ko'});
  r.plotTheseWorldPoints(trackersFromTemplate, {'ro'});
  r.plotTheseWorldPoints(noisyTrackersFromTemplate, {'go'});
  r.plotTheseWorldPoints(noisyMissingTrackersFromTemplate, {'bo'});
  r.plotTheseWorldPoints(noisyMissingTrackersFromTemplateSVD, {'mo'});
  axis equal;
end

fprintf('time taken per frame = %gs\n', totalTime / numTimePoints);
fprintf('time taken per frame for SVD method = %gs\n', totalTimeSVDMethod / numTimePoints);

subplot(2, 2, 3:4);
plot(noisyPositionTrajectory', 'g');
plot(noisyMissingPositionTrajectory', 'b');
plot(noisyMissingPositionTrajectorySVD', 'm');


return;

%% interpolate gaps in rigid body estimates:

proofMarkerTrajectory = pointsInWorld;

markerTrajectoryWithGaps = noisyPointsInWorld;
for i = 1:length(markerTrajectoryWithGaps)
  indecesToKnockOut = randperm(4);
  markerTrajectoryWithGaps{i}(:, indecesToKnockOut(1:2)) = NaN;
end

tic
[interpolatedPositions, interpolatedOrientations, interpolatedMarkersEstimatedfromTemplate] = ...
  r.interpolateRigidBodyFromMarkers(markerTrajectoryWithGaps);
toc

%%
for i = 1:numTimePoints
  subplot(2, 2, [1]);
  r.plotTheseWorldPoints(markerTrajectoryWithGaps{i}, {'ro'});
  r.position = positionTrajectory(:, i);
  r.orientation = orientationTrajectory(:, i);
  r.plotWireframesInWorld();
  subplot(2, 2, 2);
  r.plotTheseWorldPoints(pointsInWorld{i}, {'ko'});
  r.plotTheseWorldPoints(interpolatedMarkersEstimatedfromTemplate{i}, {'ro'});
  axis equal;
end

subplot(2, 2, 3:4);
plot(interpolatedPositions', 'r');










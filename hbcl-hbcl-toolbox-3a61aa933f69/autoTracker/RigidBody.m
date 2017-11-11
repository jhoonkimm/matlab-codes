classdef RigidBody
  
  properties
    position = [0 0 0]';
    orientation = [0 0 0]';
    trackersInBodyFrame = [];
  end
  
  methods
    
    %     function [this] = normalizeTrackers(this)
    %       this.trackersInBodyFrame = this.trackersInBodyFrame - ...
    %         repmat(mean(this.trackersInBodyFrame, 2), [1, size(this.trackersInBodyFrame, 2)]);
    %     end
    
    function [this] = setStateFromPositionAndRStructure(this, positionAndRStructure)
      %%
      this.position = positionAndRStructure.position;
      this.position = reshape(this.position, [3, 1]);
      
      if any(size(positionAndRStructure.R) ~= [3, 3])
        error('given R must be a 3x3 rotation matrix');
      end
      
      this.orientation = rodrigues(positionAndRStructure.R);
    end
    
    function ret = transformPointFromBodyFrameToWorldFrame(this, point)
      %       this = this.normalizeTrackers();
      rotatedPoint = rodrigues(this.orientation) * point;
      ret = rotatedPoint + this.position;
    end
    
    function ret = transformPointFromWorldFrameToBodyFrame(this, point)
      %       this = this.normalizeTrackers();
      translatedPoint = point - this.position;
      ret = rodrigues(this.orientation)' * translatedPoint;
    end
    
    function [this] = addTrackerPositionInBodyFrame(this, point)
      %       this = this.normalizeTrackers();
      this.trackersInBodyFrame(:, size(this.trackersInBodyFrame, 2) + 1) = point;
      %       this = this.normalizeTrackers();
    end
    
    function [this] = addTrackerPositionInWorldFrame(this, point)
      %       this = this.normalizeTrackers();
      this = this.addTrackerPositionInBodyFrame( ...
        this.transformPointFromWorldFrameToBodyFrame(point));
    end
    
    function ret = getTrackerPositionsInWorld(this)
      %       this = this.normalizeTrackers();
      normalizedTrackersInBodyFrame = this.trackersInBodyFrame - ...
        repmat(mean(this.trackersInBodyFrame, 2), [1, size(this.trackersInBodyFrame, 2)]);
      
      ret = normalizedTrackersInBodyFrame;
      R = rodrigues(this.orientation);
      for i = 1:size(this.trackersInBodyFrame, 2)
        %         ret(:, i) = ...
        %           this.transformPointFromBodyFrameToWorldFrame(this.trackersInBodyFrame(:, i));
        rotatedPoint = R * normalizedTrackersInBodyFrame(:, i);
        ret(:, i) = rotatedPoint + this.position;
      end
    end
    
    %% data filling and template fitting methods:
    function [position, orientation, trackersFromTemplateToReturn, fitError] = ...
        calculatePositionAndOrientationFromTrackersInWorldLeastSquares(this, measuredTrackerLocations, varargin)

      shouldPlotFit = 0;
      
      for i = 1:2:length(varargin)
        option = varargin{i};
        value = varargin{i + 1};
        switch (option)
          case 'shouldPlotFit'
            shouldPlotFit = value;
        end
      end
      
      %%
      % well, crap, look at:
      % http://kwon3d.com/theory/jkinem/rotmat.html
      %
      % let's look at equation 22, which forward calculates this,
      % and work backward to see if we can implement it.
      % (spoiler, it works...)

      normalizedTrackersInBodyFrame = this.trackersInBodyFrame - ...
        repmat(mean(this.trackersInBodyFrame, 2), [1, size(this.trackersInBodyFrame, 2)]);

      
      this.position = [0 0 0]';
      this.orientation = [0 0 0]';
      
      templateTrackerLocations = normalizedTrackersInBodyFrame; %getTrackerPositionsInWorld();
      
      markerIndecesToUse = ~isnan(measuredTrackerLocations(1,:));
      
      try
        x1 = templateTrackerLocations(:, markerIndecesToUse);
      catch e
        e
      end
      x2 = measuredTrackerLocations(:, markerIndecesToUse);
      
      x1 = x1 - repmat(mean(x1, 2), [1, size(x1, 2)]);
      x2 = x2 - repmat(mean(x2, 2), [1, size(x2, 2)]);
      
      % eq 14:
      c = zeros(3);
      for i = 1:size(x1, 2)
        c = c + x2(:, i) * x1(:, i)';
        %         c = c + x1(:, i) * x2(:, i)';
      end
      c = c / (size(x1, 2));
      
      % eq 16:
      try
        [u, w, v] = svd(c);
      catch e
        e
      end
      
      
      % eq 22:
      R = u * [1 0 0; 0 1 0; 0 0 det(u*v')] * v';
      %       R = u * v';
      %       R = R';
      
      orientation = rodrigues(R);
      %       position = mean(measuredTrackerLocations(:, markerIndecesToUse), 2) - ...
      %         mean(this.trackersInBodyFrame(:, markerIndecesToUse), 2);
      %       position = mean(this.trackersInBodyFrame(:, markerIndecesToUse), 2) - ...
      %         mean(measuredTrackerLocations(:, markerIndecesToUse), 2);
      
      % inlined for efficiency (bleh)
      trackersFromTemplate = this.trackersInBodyFrame;
      %       trackersFromTemplate = trackersFromTemplate - mean(trackersFromTemplate, 2);
      trackersFromTemplate = trackersFromTemplate - ...
        repmat(mean(trackersFromTemplate(:, markerIndecesToUse), 2), [1, size(trackersFromTemplate, 2)]);
      
      if (shouldPlotFit)
        figure
        %         this.plotTheseWorldPoints(mean(trackersFromTemplate(:, markerIndecesToUse), 2), {'*k'}, 'plotWireframes', 0);
        %         this.plotTheseWorldPoints(mean(measuredTrackerLocations(:, markerIndecesToUse), 2), {'*r'}, 'plotWireframes', 0);
      end

      % for this to work, we have to rotate the markers around the mean of
      % only the points we are using
      originOfTheFoundRotation = mean(trackersFromTemplate(:, markerIndecesToUse), 2);
      position = mean(measuredTrackerLocations(:, markerIndecesToUse), 2) - ...
        originOfTheFoundRotation; % - mean(trackersFromTemplate, 2);
      %       trackersFromTemplateOriginal = trackersFromTemplate;
      
      subsetOfTemplateMarkersForTransform = markerIndecesToUse;
      
      for i = 1:size(this.trackersInBodyFrame, 2)
        %         rotatedPoint = R * (this.trackersInBodyFrame(:, i) - mean(trackersFromTemplate, 2));
        rotatedPoint = R * trackersFromTemplate(:, i);
        trackersFromTemplate(:, i) = rotatedPoint + position;
      end
      trackersFromTemplateToReturn = trackersFromTemplate;
      
      %       % but to regenerate the points later, we need to rotate about the
      %       % midpoint of all the markers?!?! I suspect this is bad, and we need to fix the rotation matrix
      %       % since we are rotating about a new point in that case. Damnit?
      %       position = mean(measuredTrackerLocations(:, markerIndecesToUse), 2) - mean(trackersFromTemplateOriginal, 2);
      
      
      % At the moment... This is slow...
      %       fitSquareError = this.getTotalSquaredTrackerErrorOfStateEstimate( ...
      %         measuredTrackerLocations, position, orientation);
      
      %       trackersFromTemplate
      %       measuredTrackerLocations
      %       nanmean(nanmean(max(abs(trackersFromTemplate - measuredTrackerLocations))))
      
      
      fitError = sqrt(nansum(nansum((trackersFromTemplate - measuredTrackerLocations).^2)));
      numNonNaNMarkers = sum(~isnan(sum(trackersFromTemplate - measuredTrackerLocations)));
      fitError = fitError / numNonNaNMarkers;
      
      %       if (numNonNaNMarkers ~= 6)
      %         fprintf('test\n');
      %       end
      
      if (shouldPlotFit)
        this.plotTheseWorldPoints(trackersFromTemplate, {'ok'}, 'plotWireframes', 1);
        this.plotTheseWorldPoints(measuredTrackerLocations, {'or'}, 'plotWireframes', 1);
        axis equal;
      end
      
      % ok, now for an important bit... we need to change the translation and rotation that we return to be 
      % with respect to the mean point of the template markers, not the
      % mean point of the template markers that were used for the fit.
      % There are probably great, fast ways of doing this, which we will
      % not attempt or even look for. Here's a slow way:
      
      % we have the transform that went between the points, which is
      % about the origin of a subset of the points:
      rotationOriginal = R;
      originOfOriginalRotation = [0 0 0]'; %originOfTheFoundRotation; WRT used subset mean
      translationOriginal = position;
      
      % we have a stack of transformed points
      transformedTrackersFromTemplate = trackersFromTemplate;
      
      % and the same points, before transform
      trackersFromTemplate = this.trackersInBodyFrame;
      %       trackersFromTemplate = trackersFromTemplate - ...
      %         repmat(mean(trackersFromTemplate(:, markerIndecesToUse), 2), [1, size(trackersFromTemplate, 2)]);
      originOfNewTransform = mean(trackersFromTemplate, 2);
      
      position = mean(transformedTrackersFromTemplate, 2); % - originOfNewTransform;
      
      % we want to find and T2
      % xTransformed = R2 * (xRaw - originOfNewTransform) + newTranslation
      
      %       xTransformed = reshape(transformedTrackersFromTemplate(:, 1:3), [9 1]); % x,y,z of 3 transformed points.
      
      %       xRaw = reshape(trackersFromTemplate(:, 1:3), [9 1]); % x,y,z of 3 raw points.
      
      % R2: [col1 col2 col3]
      
      %       x1 = trackersFromTemplate;
      %       x2 = transformedTrackersFromTemplate;
      %
      %       x1 = x1 - repmat(mean(x1, 2), [1, size(x1, 2)]);
      %       x2 = x2 - repmat(mean(x2, 2), [1, size(x2, 2)]);
      %
      %       % eq 14:
      %       c = zeros(3);
      %       for i = 1:size(x1, 2)
      %         c = c + x2(:, i) * x1(:, i)';
      %       end
      %       c = c / (size(x1, 2));
      %       % eq 16:
      %       try
      %         [u, w, v] = svd(c);
      %       catch e
      %         e
      %       end
      %       % eq 22:
      %       R = u * [1 0 0; 0 1 0; 0 0 det(u*v')] * v';
      %       orientation = rodrigues(R);
      
    end
    
    function [position, orientation, trackersFromTemplate, fitSquareError] = ...
        calculatePositionAndOrientationFromTrackersInWorld(this, trackerLocations)
      %%
      %       this = this.normalizeTrackers();
      
      error('have to check this to make sure trackers are normalized correctly');
      functionToMinimize = @(x) this.getTotalSquaredTrackerErrorOfStateEstimate( ...
        trackerLocations, x(1:3), x(4:6));
      %       initialStateGuess = [0 0 0 0 0 0]';
      
      initialStateGuess = [this.position; this.orientation];
      functionToMinimize(initialStateGuess);
      
      %       optimizationSettings = optimset(...
      %         'Display', 'off', ...
      %         'LargeScale', 'off'); %'algorithm', 'sqp');
      
      %       optimizationSettings = optimset(...
      %         'TolFun', 2e-3, ...
      %         'TolX', 2e-3, ...
      %         'Display', 'off', ...
      %         'LargeScale', 'off'); %'algorithm', 'sqp');
      
      optimizationSettings = optimset( ...
        'TolFun', 2e-5, ...
        'TolX', 2e-5, ...
        ...'MaxFunEvals', 2000, ...
        'Display', 'off', ...
        'LargeScale', 'off'); %'algorithm', 'sqp');
      
      
      [optimalState, fitSquareError] = fminunc( ...
        functionToMinimize, ...
        initialStateGuess, ...
        optimizationSettings);
      
      position = optimalState(1:3);
      orientation = optimalState(4:6);
      
      this.position = position;
      this.orientation = orientation;
      
      trackersFromTemplate = this.getTrackerPositionsInWorld();
    end
    
    function [err] = getTotalSquaredTrackerErrorOfStateEstimate(this, ...
        measuredTrackerPositionsInWorld, estimatedPosition, estimatedOrientation)
      %%
      %       this = this.normalizeTrackers();
      error('check for normalization');
      this.position = estimatedPosition;
      this.orientation = estimatedOrientation;
      estimatedTrackerLocations = this.getTrackerPositionsInWorld();
      
      try
        err = nansum(nansum((measuredTrackerPositionsInWorld - estimatedTrackerLocations).^2));
      catch e
        e
      end
    end
    
    %     function [interpolatedPositions, interpolatedOrientations, markersEstimatedfromTemplate] = ...
    %         interpolateRigidBodyFromMarkersKalmanFilter(this, measuredMarkerTrajectory)
    %
    %
    %       %       function [ orientations, angularVelocities, ...
    %       %           stateEstimationCovariances, estimatedMeasurements ] = ...
    %       %           imuKalmanFilterV4( ...
    %       %           initialOrientation, initialOmega, ...
    %       %           initialOrientationEstimateVariance, initialOmegaEstimateVariance, ...
    %       %           accelerometerReadings, gyroScopeReadings, ...
    %       %           accelerometerNoiseVariance, gyroScopeNoiseVariance, ...
    %       %           angularAccelerationVariance, ...
    %       %           deltaT, g, angleCovariance)
    %
    %       %%
    %       % estimated state is the rigid body state in each frame where we don't have a
    %       % rigid body state yet.
    %       %
    %       % Measurements are:
    %       % - the locations of the trackers in each of the unknown frames.
    %       % - an estimate of the state in each frame based on the average of the one
    %       %   estimated state before and after. (note that this won't spline!)
    %       %
    %       % The simplest thing at the beginning is probably just to be lazy and
    %       % estimate every state in the trajectory, whether we have 3 markers in the
    %       % frame or not.
    %
    %
    %       rigidBodyInterpolation = zeros(6, length(measuredMarkerTrajectory));
    %       for i = 1:length(measuredMarkerTrajectory)
    %         [position, orientation, trackersFromTemplate, fitSquareError] = ...
    %           this.calculatePositionAndOrientationFromTrackersInWorld(measuredMarkerTrajectory{i});
    %         this.position = position;
    %         this.orientation = orientation;
    %         rigidBodyInterpolation(1:3, i) = position;
    %         rigidBodyInterpolation(4:6, i) = orientation;
    %       end
    %
    %       rigidBodyInterpolation = reshape(numel(rigidBodyInterpolation), 1);
    %
    %       %g = 9.81;
    %       %       dt = deltaT;
    %
    %       %         if (~initialOrientationEstimateVariance)
    %       %           %   fprintf('setting default initialOrientationEstimateVariance');
    %       %           initialOrientationEstimateVariance = (((dt^2) / 2)^2) * angularAccelerationVariance; %0.01; %
    %       %           %     initialOrientationEstimateVariance = angleCovariance * angularAccelerationVariance; %(((dt^2) / 2)^2) * angularAccelerationVariance;
    %       %         end
    %       %
    %       %         if (~initialOmegaEstimateVariance)
    %       %           %   fprintf('setting default initialOmegaEstimateVariance');
    %       %           initialOmegaEstimateVariance = (dt^2) * angularAccelerationVariance; %
    %       %           %     initialOmegaEstimateVariance = angularAccelerationVariance; %(dt^2) * angularAccelerationVariance; %
    %       %         end
    %
    %
    %       initialCovarianceEstimate = blkdiag(eye(4,4) .* initialOrientationEstimateVariance, ...
    %         eye(3,3) .* initialOmegaEstimateVariance);
    %
    %       initialState = [initialOrientation; initialOmega];
    %
    %       % processNoiseCovarianceMatrix = blkdiag((((dt^2) / 2)^2) * angularAccelerationVariance .* eye(4,4), ...
    %       %   (dt^2) * angularAccelerationVariance .* eye(3,3));
    %
    %       processNoiseCovarianceMatrix = blkdiag(((dt^2) / 2) * angularAccelerationVariance .* eye(4,4), ...
    %         (dt) * angularAccelerationVariance .* eye(3,3));
    %
    %       % processNoiseCovarianceMatrix = blkdiag(angleCovariance * angularAccelerationVariance .* eye(4,4), ...
    %       %     angularAccelerationVariance .* eye(3,3));
    %
    %       measurementNoiseCovarianceMatrix = blkdiag(eye(3,3) .* accelerometerNoiseVariance, ...
    %         eye(3,3) .* gyroScopeNoiseVariance);
    %
    %       %%
    %       % State Vector is:
    %       %
    %       % $$x = \left[ \begin[array][c]
    %       % \vec[\epsilon] \\
    %       % \vec[\omega]
    %       % \end[array] \right]$$
    %       %
    %
    %       %% State Transition Function:
    %       %
    %       %
    %
    %       function [xNew] = stateTransitionFunction (x, u)
    %         q0 = x(1);
    %         q1 = x(2);
    %         q2 = x(3);
    %         q3 = x(4);
    %
    %         wx = x(5);
    %         wy = x(6);
    %         wz = x(7);
    %
    %         magOmega = sqrt(wx*wx + wy*wy + wz*wz);
    %
    %         if(isnan(magOmega))
    %           error('magnitude of omega is NaN');
    %         end
    %
    %         S = sin(0.5 *  deltaT * magOmega);
    %         C = cos(0.5 *  deltaT * magOmega);
    %         magQ = q0^2 + q1^2 + q2^2 + q3^2;
    %
    %         xNew = [(C*magOmega*q0 + S*(q1*wx + q2*wy + q3*wz))/(magOmega* ...
    %           sqrt(magQ)), (C*magOmega*q1 -  ...
    %           S*(q0*wx - q3*wy + q2*wz))/(magOmega* ...
    %           sqrt(magQ)), (C*magOmega*q2 -  ...
    %           S*(q3*wx + q0*wy - q1*wz))/(magOmega*sqrt(magQ)),  ...
    %           (C*magOmega*q3 + S*(q2*wx - q1*wy - q0*wz))/(magOmega* ...
    %           sqrt(magQ)), wx, wy, wz]';
    %
    %
    %       end
    %
    %       %%
    %       %
    %       function [jacobian] = stateTransitionJacobian(x, u)
    %         q0 = x(1);
    %         q1 = x(2);
    %         q2 = x(3);
    %         q3 = x(4);
    %
    %         wx = x(5);
    %         wy = x(6);
    %         wz = x(7);
    %
    %         magOmega = sqrt(wx*wx + wy*wy + wz*wz);
    %
    %         S = sin(0.5 * deltaT * magOmega);
    %         C = cos(0.5 * deltaT * magOmega);
    %         magQ = q0^2 + q1^2 + q2^2 + q3^2;
    %
    %
    %         % this came from the mathematica file: kalmanFilterEquations.nb
    %
    %         %           jacobian = [[]];
    %
    %
    %         jacobian = reshape(jacobian, 7,7)';
    %
    %       end
    %
    %       stateTransitionFunctionSet.valueAt = @stateTransitionFunction;
    %       stateTransitionFunctionSet.jacobianWithRespectToState = @stateTransitionJacobian;
    %
    %       function [z] = measurementFunction(x)
    %         q0 = x(1);
    %         q1 = x(2);
    %         q2 = x(3);
    %         q3 = x(4);
    %         wx = x(5);
    %         wy = x(6);
    %         wz = x(7);
    %
    %         qMag = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    %
    %         q0 = q0 / qMag;
    %         q1 = q1 / qMag;
    %         q2 = q2 / qMag;
    %         q3 = q3 / qMag;
    %
    %         magOmega = sqrt(wx*wx + wy*wy + wz*wz);
    %
    %         S = sin(0.5 * deltaT * magOmega);
    %         C = cos(0.5 * deltaT * magOmega);
    %         magQ = q0^2 + q1^2 + q2^2 + q3^2;
    %
    %         %           z = [
    %         %             ]';
    %       end
    %
    %
    %       function [jacobian] = measurementFunctionJacobian(x)
    %         q0 = x(1);
    %         q1 = x(2);
    %         q2 = x(3);
    %         q3 = x(4);
    %
    %         wx = x(5);
    %         wy = x(6);
    %         wz = x(7);
    %
    %         %           jacobian =
    %         ...[NaN, NaN, NaN, NaN, NaN, 5, 10]]; %
    %
    %
    %       jacobian = reshape(jacobian, 7, 6)';
    %       end
    %
    %       measurementFunctionSet.valueAt = @measurementFunction;
    %       measurementFunctionSet.jacobianWithRespectToState = @measurementFunctionJacobian;
    %
    %       measurements = [accelerometerReadings; gyroScopeReadings];
    %       inputs = ones(size(measurements)) .* 1;
    %
    %       [ estimatedStates, stateEstimationCovariances, estimatedMeasurements] = kalmanFilter( ...
    %         initialState, initialCovarianceEstimate, ...
    %         stateTransitionFunctionSet, measurementFunctionSet, ...
    %         inputs, measurements, processNoiseCovarianceMatrix, ...
    %         measurementNoiseCovarianceMatrix);
    %
    %
    %       estimatedMeasurements = estimatedMeasurements';
    %       orientations = estimatedStates(:,1:4)';
    %       angularVelocities = estimatedStates(:,5:7)';
    %       % covarianceMatrices = stateEstimationCovariances;
    %
    %     end
    
    
    function [interpolatedPositions, interpolatedOrientations, markersEstimatedfromTemplate] = ...
        interpolateRigidBodyFromMarkers(this, measuredMarkerTrajectory)
      %%
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
      
      %       this = this.normalizeTrackers();
      error('check for normalization');

      
      
      rigidBodyInterpolation = zeros(6, length(measuredMarkerTrajectory));
      numConstraints = 0;
      for i = 1:length(measuredMarkerTrajectory)
        [position, orientation, trackersFromTemplate, fitSquareError] = ...
          this.calculatePositionAndOrientationFromTrackersInWorld(measuredMarkerTrajectory{i});
        this.position = position;
        this.orientation = orientation;
        rigidBodyInterpolation(1:3, i) = position;
        rigidBodyInterpolation(4:6, i) = orientation;
        
        %         if (sum(sum(isnan(measuredMarkerTrajectory{i}))) == 0)
        %           numConstraints = numConstraints + 1;
        %
        %           %           awwww, fuck it.
        %           %           Add the locations of the points to the estimated state, then the existing points have to be met exactly. Or, probably better, add a
        %           %           slowly changing rigid body to the state. So I guess a rigid body at each state? Or for every N states.
        %
        %           r = RigidBody;
        %           r.position = wholeStateEstimate(1:3, indecesIntoAllData(1)) / 1000.0;
        %           r.orientation = wholeStateEstimate(4:6, indecesIntoAllData(1));
        %           for j = 1:length(markerNames)
        %             r = r.addTrackerPositionInWorldFrame(rawMarkers{indecesIntoAllData(1)}(:,j) / 1000.0);
        %           end
        %         end
        
      end
      
      functionToMinimize = @(x) this.rigidBodyInterpolationError(measuredMarkerTrajectory, x);
      
      %       optimizationSettings = optimset(...
      %         ...'TolFun', 2e-5, ...
      %         ...'TolX', 2e-5, ...
      %         'TolFun', 2e-3, ...
      %         'TolX', 2e-3, ...
      %         ...'MaxFunEvals', 2000, ...
      %         'Display', 'Iter', ...
      %         'LargeScale', 'off'); %'algorithm', 'sqp');
      
      %       [rigidBodyInterpolation, fitSquareError] = fminunc( ...
      %         functionToMinimize, ...
      %         rigidBodyInterpolation, ...
      %         optimizationSettings);
      %
      
      %       optimizationSettings = optimset(...
      %         ...'TolFun', 2e-5, ...
      %         ...'TolX', 2e-5, ...
      %         'TolFun', 2e-3, ...
      %         'TolX', 2e-3, ...
      %         ...'MaxFunEvals', 2000, ...
      %         'Display', 'Iter', ...
      %         'algorithm', 'sqp');
      %
      %       constraintMatrix = zeros(1, numel(rigidBodyInterpolation)); %[]; %
      %       constraintVector = 0; %[]; %
      %       [rigidBodyInterpolation, fitSquareError] = fmincon( ...
      %         functionToMinimize, ...
      %         rigidBodyInterpolation, ...
      %         constraintMatrix, constraintVector, ...
      %         [], [], ...
      %         [], [], ...
      %         [], ...
      %         optimizationSettings);
      
      optimizationSettings = optimset(...
        ...'TolFun', 2e-5, ...
        ...'TolX', 2e-5, ...
        'TolFun', 2e-3, ...
        'TolX', 2e-3, ...
        ...'MaxFunEvals', 2000, ...
        'Display', 'Iter', ...
        'LargeScale', 'off'); %'algorithm', 'sqp');
      [rigidBodyInterpolation, fitSquareError] = fminsearch( ...
        functionToMinimize, ...
        rigidBodyInterpolation, ...
        optimizationSettings);
      
      interpolatedPositions = rigidBodyInterpolation(1:3, :);
      interpolatedOrientations = rigidBodyInterpolation(4:6, :);
      
      markersEstimatedfromTemplate = measuredMarkerTrajectory;
      for i = 1:length(measuredMarkerTrajectory)
        this.position = interpolatedPositions(:,i);
        this.orientation = interpolatedOrientations(:,i);
        markersEstimatedfromTemplate{i} = this.getTrackerPositionsInWorld();
      end
      
      
    end
    
    
    function [err] = rigidBodyInterpolationError(this, measuredMarkerTrajectory, rigidBodyInterpolation)
      err = 0;
      markerToRigidBodySmoothnessRatio = 0.01; %0.5; %
      
      for i = 1:length(measuredMarkerTrajectory)
        %         this.position = rigidBodyInterpolation(1:3, i);
        %         this.orientation = rigidBodyInterpolation(4:6, i);
        markerEstimationError = getTotalSquaredTrackerErrorOfStateEstimate(this, ...
          measuredMarkerTrajectory{i}, rigidBodyInterpolation(1:3, i), rigidBodyInterpolation(4:6, i));
        err = err + (markerEstimationError * (markerToRigidBodySmoothnessRatio));
      end
      
      for i = 2:(length(measuredMarkerTrajectory) - 1)
        smoothedRigidBodyState = mean([rigidBodyInterpolation(:, i-1) rigidBodyInterpolation(:, i+1)], 2);
        smoothRigidBodyError = sum((smoothedRigidBodyState - rigidBodyInterpolation(:, i)).^2);
        err = err + (smoothRigidBodyError * (1 - markerToRigidBodySmoothnessRatio));
      end
    end
    
    %% find corresponding trackers:
    function [correspondences, rootMeanSquaredError] = findCorrespondecesOfThisBodyToGivenPoints(this, pointsInCloud, varargin)
      
      verbose = 0;
      minSquaredErrorLimit = 0.1;
      constrainedCorrespondences = [];
      for i = 1:2:length(varargin)
        if (strcmp('minSquaredErrorLimit', varargin{i}))
          minSquaredErrorLimit = varargin{i + 1};
        end
        if (strcmp('verbose', varargin{i}))
          verbose = varargin{i + 1};
        end
        if (strcmp('constrainedCorrespondences', varargin{i}))
          constrainedCorrespondences = varargin{i + 1};
        end
      end
      
      done = 0;
      cloudBody = RigidBody;
      cloudBody.trackersInBodyFrame = pointsInCloud;
      templateBody = this;
      pointsInTemplate = templateBody.trackersInBodyFrame;
      for numberOfPointsToUseFromTemplate = size(pointsInTemplate, 2):-1:3
        for numberOfPointsToUseFromCloud = numberOfPointsToUseFromTemplate:-1:3
          
          %           fprintf('size(pointsInTemplate, 2) = %g\n', size(pointsInTemplate, 2));
          
          setOfIndecesToTryFromTemplate = ...
            getCombinationsOfListElements(1:size(pointsInTemplate, 2), numberOfPointsToUseFromTemplate, ...
            'getSortedUniqueCombinations', 1);
          
          possibleIndecesInCloud = 1:size(pointsInCloud, 2);
          % if there are fewer cloud points then template points, we have to be
          % able to add in nan points corresponding to the template points. These
          % will be signified with negative indeces:
          possibleIndecesInCloud = [possibleIndecesInCloud ...
            -1:-1:(numberOfPointsToUseFromCloud - numberOfPointsToUseFromTemplate)];
          
          setOfIndecesToTryFromCloud = ...
            getCombinationsOfListElements(possibleIndecesInCloud, numberOfPointsToUseFromTemplate, ...
            'getSortedUniqueCombinations', 0);
          
          allResidualErrors = ...
            zeros(size(setOfIndecesToTryFromTemplate, 2), size(setOfIndecesToTryFromCloud, 2)) * NaN;
          shouldBreakOut = 0;
          for templatesIndecesNumber = 1:size(setOfIndecesToTryFromTemplate, 2);
            
            templateIndeces = setOfIndecesToTryFromTemplate(:, templatesIndecesNumber);
            
            subTemplateRigidBody = RigidBody;
            subTemplateRigidBody.trackersInBodyFrame = templateBody.trackersInBodyFrame(:, templateIndeces);
            for cloudIndecesNumber = 1:size(setOfIndecesToTryFromCloud, 2)
              
              cloudIndeces = setOfIndecesToTryFromCloud(:, cloudIndecesNumber);
              
              thisPotentialCorrespondenceSatisfiesConstraints = 1;
              for constraintNumber = 1 : size(constrainedCorrespondences, 1)
                thisConstraint = constrainedCorrespondences(constraintNumber, :);
                %                 if (sum(templateIndeces == thisConstraint(1)) == 0)
                inds = templateIndeces == thisConstraint(1);
                if (~any(inds))
                  thisPotentialCorrespondenceSatisfiesConstraints = 0;
                  break;
                end
                if (cloudIndeces(inds) ~= thisConstraint(2))
                  thisPotentialCorrespondenceSatisfiesConstraints = 0;
                  break;
                end
              end
              if (~thisPotentialCorrespondenceSatisfiesConstraints)
                continue;
              end
              
              %               fprintf('this set of [template : cloud] indeces satisfy the constraints... \n');
              %               [templateIndeces cloudIndeces]
              
              if (mod(cloudIndecesNumber, 10000) == 0)
                if (verbose > 0)
                  fprintf('attempting fit for template indeces [%s], %g%% done\n', ...
                    sprintf('%g, ', templateIndeces), ...
                    100 * cloudIndecesNumber /size(setOfIndecesToTryFromCloud, 2));
                end
              end
              
              subCloudRigidBody = RigidBody;
              
              subCloudRigidBody.trackersInBodyFrame = NaN * zeros(3, length(cloudIndeces));
              subCloudRigidBody.trackersInBodyFrame(:, find(cloudIndeces > 0)) = ...
                cloudBody.trackersInBodyFrame(:, cloudIndeces(cloudIndeces > 0));
              
              [position, orientation, trackersFromTemplate, fitSquareError] = ...
                subTemplateRigidBody.calculatePositionAndOrientationFromTrackersInWorldLeastSquares( ...
                subCloudRigidBody.trackersInBodyFrame);
              
              if (isnan(fitSquareError))
                error('pretty sure we shouldn''t have a nan error');
              end
              
              % fprintf('fitError = %gmm\n', fitSquareError);
              
              allResidualErrors(templatesIndecesNumber, cloudIndecesNumber) = fitSquareError;
              

              %                 shouldBreakOut = 1;
              %                 break;
              %               end
            end
            %             if (shouldBreakOut)
            %               break;
            %             end
          end
          
          if (size(allResidualErrors, 1) > 1)
            [minVal, minCol] = min(min(allResidualErrors));
            [minVal, minRow] = min(allResidualErrors(:, minCol));
          else
            minRow = 1;
            [minVal, minCol] = min(allResidualErrors);
          end
          
          bestTemplateIndeces = setOfIndecesToTryFromTemplate(:, minRow);
          bestCloudIndeces = setOfIndecesToTryFromCloud(:, minCol);
          
          %           minVal = sqrt(minVal) / (min(sum(bestCloudIndeces > 0), sum(bestTemplateIndeces > 0)));
          
          if (verbose > 0)
            fprintf('found best correspondence error of %g\n', minVal);
            [bestTemplateIndeces bestCloudIndeces]'
          end
          
          %           figure;
          %           plot(allResidualErrors');
          %           title('residual squared fit errors for correspondence calculation');
          fprintf('best error = %gmm\n', minVal);
          
          if (isnan(minVal))
            error('pretty sure we shouldn''t have a nan error');
          end
          
          
          if (minVal < minSquaredErrorLimit)
            bestTemplateIndeces = setOfIndecesToTryFromTemplate(:, minRow);
            bestCloudIndeces = setOfIndecesToTryFromCloud(:, minCol);
            done = 1;
          end
          if (done)
            break;
          end
        end
        
        if (done)
          break;
        end
      end
      
      goodCorrespondenceIndeces = find(bestCloudIndeces > 0);
      correspondences = [bestTemplateIndeces(goodCorrespondenceIndeces) ...
        bestCloudIndeces(goodCorrespondenceIndeces)];
      
      rootMeanSquaredError = minVal;
      
      %       if (verbose > 0)
      %         fprintf('found min error of %g, correspondence indeces are:\n', minVal);
      %         correspondences
      %       end
      
      if (verbose > 1)
        figure;
        plot(allResidualErrors');
        title('residual squared fit errors for correspondence calculation');
      end
      
    end
    
    %%
    function [closestPointIndexInCloud, closestPointIndexInTemplate] = ...
        findClosestPointInCloudToAnyPointInBody(this, cloud)
      
      thesePoints = this.trackersInBodyFrame;
      
      expandedCloud = repmat(cloud, [1 size(thesePoints, 2)]);
      expandedThesePoints = ...
        reshape(repmat(reshape(thesePoints', [1 numel(thesePoints)]), [size(cloud, 2) 1]), ...
        size(expandedCloud'))';
      
      differences = sum((expandedCloud - expandedThesePoints).^2);
      differences = reshape(differences, [size(cloud, 2), size(thesePoints, 2)])';
      [minValues, minCols] = min(differences);
      
      [minVal, minRow] = min(minValues);
      
      closestPointIndexInCloud = minRow;
      closestPointIndexInTemplate = minCols(minRow);
    end
    
    %%
    function [closestIndeces] = ...
        getClosestNPointsToPointsInThisBody( ...
        this, cloudPoints, numberOfClosestPointsToFind, varargin)
      
      shouldPlot = 0;
      for i = 1:2:length(varargin)
        if (strcmp('shouldPlot', varargin{i}))
          shouldPlot = varargin{i + 1};
        end
      end
      
      %       rigidBody = rigidBodyChain.rigidBodySegments.pelvis;
      rigidBody = this;
      %       numberOfClosestPointsToFind = size(rigidBody.trackersInBodyFrame, 2);
      closeCloudPoints = zeros(3, numberOfClosestPointsToFind);
      closeTemplatePoints = zeros(3, numberOfClosestPointsToFind);
      possibleCorrespondingPoints = cloudPoints;
      templatePoints = rigidBody.trackersInBodyFrame;
      
      if (shouldPlot)
        figure;
        plotLine = @(a, b, color, lineWidth) line([a(1) b(1)], ...
          [a(2) b(2)], [a(3) b(3)], 'color', color, 'LineWidth', lineWidth);
      end
      
      closestIndeces = false(1, size(possibleCorrespondingPoints, 2));
      for i = 1:numberOfClosestPointsToFind
        
        allIndeces = 1:size(possibleCorrespondingPoints, 2);
        allIndeces(closestIndeces) = [];
        
        if (isempty(allIndeces))
          break;
        end
        
        [closestPointIndexInCloud, closestPointIndexInTemplate] = ...
          rigidBody.findClosestPointInCloudToAnyPointInBody( ...
          possibleCorrespondingPoints(:, allIndeces));
        
        closestPointIndexInCloud = allIndeces(closestPointIndexInCloud);
        closeCloudPoints(:, i) = possibleCorrespondingPoints(:, closestPointIndexInCloud);
        closeTemplatePoints(:, i) = templatePoints(:, closestPointIndexInTemplate);
        %         possibleCorrespondingPoints(:, closestPointIndexInCloud) = [];
        closestIndeces(closestPointIndexInCloud) = 1;
        
        if (shouldPlot)
          plotLine(closeTemplatePoints(:, i), closeCloudPoints(:, i), 'k', 2);
          hold on;
        end
      end
      
      closestIndeces = find(closestIndeces);
      
      if (shouldPlot)
        rigidBody.plotTheseWorldPoints(templatePoints, {'ok'});
        rigidBody.plotTheseWorldPoints(cloudPoints, {'or'});
        axis equal
      end
    end
    
    %% plotting methods:
    function [] = plotTheseWorldPoints(this, points, styles, varargin)
      %       this = this.normalizeTrackers();
      plotWireframes = 0;
      
      for i = 1:2:length(varargin)
        option = varargin{i};
        value = varargin{i + 1};
        switch (option)
          case 'plotWireframes'
            plotWireframes = value;
        end
      end
      
      for i = 1:size(points, 2)
        p = points(:, i);
        plot3(p(1), p(2), p(3), styles{mod(i, length(styles)) + 1});
        hold on;
      end
      
      if (plotWireframes)
        lastPosition = points(:, 1); %mean(points, 2);
        for i = 2:size(points, 2)
          p = points(:, i);
          line([lastPosition(1) p(1)], [lastPosition(2) p(2)], [lastPosition(3) p(3)], ...
            'color', 'k');
          hold on;
          lastPosition = p;
        end
      end
      
    end
    
    
    function [] = plotTrackersInWorldStyled(this, styles)
      %       this = this.normalizeTrackers();
      
      trackersInWorld = getTrackerPositionsInWorld(this);
      this.plotTheseWorldPoints(trackersInWorld, styles);
    end
    
    function [] = plotTrackersInWorld(this)
      %       this = this.normalizeTrackers();
      
      this.plotTrackersInWorldStyled({'.r', '.g', '.b', '.c', '.m', '.y', '.k'});
    end
    
    function [] = plotWireframesInWorld(this)
      %       this = this.normalizeTrackers();
      
      trackersInWorld = getTrackerPositionsInWorld(this);
      meanPosition = mean(trackersInWorld, 2);
      for i = 1:size(trackersInWorld, 2)
        p = trackersInWorld(:, i);
        line([meanPosition(1) p(1)], [meanPosition(2) p(2)], [meanPosition(3) p(3)], ...
          'color', 'k');
        hold on;
      end
    end
    
    
  end
  
end


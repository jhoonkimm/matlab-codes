function [transformationFromPhasespaceToCart, phasespaceScale, calibrationError] = ..., angleScale] = ...
  findAlignmentBetweenCartAndPhasespaceFrames(cartTimes, cartMotion, markerTimes, ...
  markerMotion, varargin)

global markerTimesGlobal markerMotionGlobal cartSpline findAlignmentVerbose

markerTimesGlobal = markerTimes;
markerMotionGlobal = markerMotion;

verbose = 0;
skipFrameRegistration = 0; %1; %
for i = 1:2:length(varargin)
  if(strcmp('verbose', varargin{i}))
    verbose = varargin{i + 1};
  end
  if (strcmp('skipFrameRegistration', varargin{i}))
    skipFrameRegistration = varargin{i + 1};
  end
end

%%

cartSpline = spline(cartTimes, cartMotion', markerTimes);

% initialScaleAndTransformation = [0.2046, 0, 0, 0];
% initialScaleAndTransformation = [-0.2046, -2000, 0, 0];
% initialScaleAndTransformation = [0.2041   -0.2940    2.2269    1.5048];
% initialScaleAndTransformation = [0.2046, 0, 0, -pi/4];
% initialScaleAndTransformation = [1, 0, 0, 0];
% initialScaleAndTransformation = [0.2046, 2000, 0, 0];

% using gyro for rad:
% initialScaleAndTransformation = 1.0e+003 * [0.000203764491361   2.292691480875553   0.000010881381698  -0.000017113278423];
% using encoders for rad:
% initialScaleAndTransformation = 1.0e+003 * [0.000202958497687   2.312438354922517   0.000001137713308  -0.000030358606314, 1.1e-3];
% initialScaleAndTransformation = 1.0e+003 * [0.000201356218096   2.136905410196281   0.000001111127455  -0.000032799955892   0.001061100609960];

% using encoders for rad and second half of peter integration test
% initialScaleAndTransformation = 1.0e+003 * [0.000187006780900   2.000598809655708   0.000001609803220  -0.000024149252845]; %   0.001031673633622];

% John walk test:
% initialScaleAndTransformation = 1.0e+003 * [0.000216559171677   2.000598754882812   0.000001609803177  -0.000019874637946];

% eye candy:
% initialScaleAndTransformation = 1.0e+003 * [0.000210137426591   3.577345548450465  -0.000008581601465   0.000002352079625];

% yufanl variability pilot
% initialScaleAndTransformation = 1.0e+003 * [0.000220130387821   3.654871134337918  -0.000010717088318  -0.000027201647440];
% initialScaleAndTransformation = 1.0e+003 * [0.0002206   3.6071421  -0.0000111  -0.0000233];

% now that we are in meters instead of encoder ticks and mm:
initialScaleAndTransformation = [1.0 3.6071421 0 0];
% initialScaleAndTransformation = [1.0 0 0 pi];


if (skipFrameRegistration)
  scaleAndTransformationFromPhasespaceToCart = initialScaleAndTransformation;
  phasespaceScale = scaleAndTransformationFromPhasespaceToCart(1);
  transformationFromPhasespaceToCart = scaleAndTransformationFromPhasespaceToCart(2:4)';
  % angleScale = scaleAndTransformationFromPhasespaceToCart(5);
  fprintf('skipping the glorious frame registration routine!\n');
  calibrationError = NaN;
  return;
end

%%
findAlignmentVerbose = verbose;
if (findAlignmentVerbose > 2)
  figure;
end

% initialScaleAndTransformation = [0.1, 0, 0, pi];
% initialScaleAndTransformationOffset = [0.1, 0, 0, pi, 0.8];

% initialScaleAndTransformation = [0.2046, 0, 0, 0, 0.8 ];
% initialScaleAndTransformation = [0, 0, 0];

options = optimset( ...
  'Display', 'off', ...
  'FinDiffType', 'central', ...
...  'Algorithm', 'lm-line-search' ...'interior-point' ... 'active-set', ...
  'LargeScale', 'off' ...
  );

initialScaleAndTransformation(1) = [];

getErrorOfTransformationFromPhasespaceToCart(initialScaleAndTransformation);
% scaleAndTransformationFromPhasespaceToCart = fminsearch(@getErrorOfTransformationFromPhasespaceToCart, initialScaleAndTransformation);
scaleAndTransformationFromPhasespaceToCart = ...
  fminunc(@getErrorOfTransformationFromPhasespaceToCart, initialScaleAndTransformation, options);
scaleAndTransformationFromPhasespaceToCart = double(scaleAndTransformationFromPhasespaceToCart);
scaleAndTransformationFromPhasespaceToCart = fminsearch(@getErrorOfTransformationFromPhasespaceToCart, ...
  scaleAndTransformationFromPhasespaceToCart, options);
[scaleAndTransformationFromPhasespaceToCart calibrationError] = ...
  fminunc(@getErrorOfTransformationFromPhasespaceToCart, scaleAndTransformationFromPhasespaceToCart, options);
% scaleAndTransformationFromPhasespaceToCart = fminsearch(@getErrorOfTransformationFromPhasespaceToCart, scaleAndTransformationFromPhasespaceToCart);
% scaleAndTransformationFromPhasespaceToCart
% transformationFromPhasespaceToCart = scaleAndTransformationFromPhasespaceToCart(2:end)';
% phasespaceScale = scaleAndTransformationFromPhasespaceToCart(1);

% getErrorOfTransformationFromPhasespaceToCart(initialScaleAndTransformationOffset);
% scaleAndTransformationOffsetFromPhasespaceToCart = ...
%   fminunc(@getErrorOfTransformationFromPhasespaceToCart, initialScaleAndTransformationOffset);
% scaleAndTransformationOffsetFromPhasespaceToCart = double(scaleAndTransformationOffsetFromPhasespaceToCart);
% scaleAndTransformationOffsetFromPhasespaceToCart = fminsearch(@getErrorOfTransformationFromPhasespaceToCart, scaleAndTransformationOffsetFromPhasespaceToCart);
% scaleAndTransformationOffsetFromPhasespaceToCart = fminunc(@getErrorOfTransformationFromPhasespaceToCart, scaleAndTransformationOffsetFromPhasespaceToCart);
% scaleAndTransformationOffsetFromPhasespaceToCart

% phasespaceScale = scaleAndTransformationFromPhasespaceToCart(1);
% transformationFromPhasespaceToCart = scaleAndTransformationFromPhasespaceToCart(2:4)';

phasespaceScale = 1;
transformationFromPhasespaceToCart = scaleAndTransformationFromPhasespaceToCart';


if (verbose > 0)
  fprintf('calculated calibration: [%s] with error %g\n', ...
    sprintf('%g, ', scaleAndTransformationFromPhasespaceToCart), calibrationError);
end

%%

end


function [errVal] = getErrorOfTransformationFromPhasespaceToCart(scaleAndTransformation)

frameSkip = 50; %20; %100; %10; %200; %

global markerTimesGlobal markerMotionGlobal cartSpline findAlignmentVerbose

markerTimes = markerTimesGlobal;
markerMotion = markerMotionGlobal;

% phasespaceScale = scaleAndTransformation(1);
% transformationFromPhasespaceToCart = scaleAndTransformation(2:4)';
phasespaceScale = 1;
transformationFromPhasespaceToCart = scaleAndTransformation';


% transformationFromPhasespaceToCart = scaleAndTransformation(2:5)';
phasespaceOffset = 0.0; %scaleAndTransformation(5);% 0.8; %
angleScale = 1.0; %scaleAndTransformation(5); 

% phasespaceScale = 0.2041;
% transformationFromPhasespaceToCart = scaleAndTransformation(1:3)';

initialPhasespacePoint = markerMotion(1, 1:2)';
cartPose = cartSpline(:,1);
phasespaceOffset = max(0.0, phasespaceOffset);
% cartPose = interp1(markerTimes, cartSpline', markerTimes(1) + phasespaceOffset)';

initialPointInWorld = transformPhasespacePointIntoWorld(initialPhasespacePoint, ...
  transformationFromPhasespaceToCart, phasespaceScale, cartPose);

calculateErrorBasedOnPosition = 0; %1; % else calculate based on velocities
if (findAlignmentVerbose > 2)
  clf;
  subplot(2,1,1);
  if (calculateErrorBasedOnPosition)
    plot(initialPointInWorld(1), initialPointInWorld(2), 'ok')
  else
    plot(0, 0, 'ok')
  end
  hold on;
  %   subplot(2,1,2);
  %   plot(cartSpline');
  %   hold on;
end

if (calculateErrorBasedOnPosition)
  errVal = 0.0;
  nanPoints = 0;
  goodPoints = 0;
  markerPointsToTransform = 1:frameSkip:length(markerTimes);
  allTransformedPoints = zeros(2, length(markerPointsToTransform));
  for markerNumberToTransform = 1:length(markerPointsToTransform)
    i = markerPointsToTransform(markerNumberToTransform);
    currentPhasespacePoint = markerMotion(i, 1:2)';
    %   currentPhasespacePoint = interp1(markerTimes, markerMotion(:, 1:2), markerTimes(i) + phasespaceOffset)';
    
    cartPose = cartSpline(:,i);
    %   cartPose = interp1(markerTimes, cartSpline', markerTimes(i) + phasespaceOffset)';
    cartPose(3) = cartPose(3) * angleScale;
    
    transformedPoint = transformPhasespacePointIntoWorld(currentPhasespacePoint, ...
      transformationFromPhasespaceToCart, phasespaceScale, cartPose);
    allTransformedPoints(:, markerNumberToTransform) = transformedPoint;
    
    if (~isnan(norm(transformedPoint - initialPointInWorld)))
      errVal = errVal + norm(transformedPoint - initialPointInWorld);
      goodPoints = goodPoints + 1;
    else
      nanPoints = nanPoints + 1;
    end
  end
  errVal = errVal / goodPoints;
else
  errVal = 0.0;
  nanPoints = 0;
  goodPoints = 0;
  markerPointsToTransform = 1:frameSkip:length(markerTimes);
  allTransformedPoints = zeros(2, length(markerPointsToTransform));
  
  %%
  cartVelocities = diff(cartSpline');
  %   cartVelocities = cartVelocities(:, 1:2);
  markerVelocities = diff(markerMotionGlobal);
  markerVelocities = markerVelocities(:, 1:2);
  
  %   cartSpeeds = sqrt(sum(cartVelocities(:, 1:2).^2, 2));
  %   markerSpeeds = sqrt(sum(markerVelocities.^2, 2));
  
  %   markerAndCartSpeedDifference = abs(cartSpeeds - markerSpeeds);
  %   cartAngularSpeed = cartVelocities(:, 3);
  %   clf;
  %   plot(markerAndCartSpeedDifference, 'r');
  %   hold on;
  %   plot(abs(cartAngularSpeed), 'g')
  %
  %   goodInds = abs(cartAngularSpeed) > 1e-4;
  %   xVals = 1:length(cartAngularSpeed);
  %   %   xVals(~goodInds) = NaN;
  %   ratio = markerAndCartSpeedDifference ./ abs(cartAngularSpeed);
  %   ratio(~goodInds) = NaN;
  %   plot(xVals, ratio, 'c');
  %
  
  %%
  
  for markerNumberToTransform = 1:(length(markerPointsToTransform)-1)
    i = markerPointsToTransform(markerNumberToTransform);
    %     currentPhasespacePoint = markerMotion(i, 1:2)';
    
    currentPhasespaceVelocity = markerVelocities(i, 1:2)';
    
    %     cartPose = cartSpline(:,i);
    
    cartVelocity = cartVelocities(i, :);
    
    velocityOfPhasespaceFrameFromCart = [cartVelocity(1:2)'; 0] + ...
      cross(-[transformationFromPhasespaceToCart(1:2); 0], [0; 0; cartVelocity(3)]);
    velocityOfPhasespaceFrameFromCart(3) = [];
    theta = transformationFromPhasespaceToCart(3);
    %     theta = -transformationFromPhasespaceToCart(3);
    velocityOfPhasespaceFrameFromCart = ...
      [cos(theta) -sin(theta);
      sin(theta) cos(theta)] ...
      * velocityOfPhasespaceFrameFromCart;
    
    
    
    velocityDiff = currentPhasespaceVelocity + velocityOfPhasespaceFrameFromCart;
    
    
    %     transformedPoint = transformPhasespacePointIntoWorld(currentPhasespacePoint, ...
    %       transformationFromPhasespaceToCart, phasespaceScale, cartPose);
    allTransformedPoints(:, markerNumberToTransform) = velocityDiff; %transformedPoint;
    
    if (~isnan(norm(velocityDiff)))
      errVal = errVal + norm(velocityDiff);
      goodPoints = goodPoints + 1;
    else
      nanPoints = nanPoints + 1;
    end
  end
  errVal = errVal / goodPoints;
end




if (findAlignmentVerbose > 2)
  subplot(2,1,1);
  plot(allTransformedPoints(1,:), allTransformedPoints(2,:), '*r');
  title(sprintf('error = %g, phasespaceOffset = %g', errVal, phasespaceOffset));
  axis equal;
  
  subplot(2,1,2);
  plot(markerPointsToTransform, allTransformedPoints, '--');
  %   hold on;
  %   plot(markerPointsToTransform, repmat(initialPointInWorld, [1, length(markerPointsToTransform)]), 'k--');
  drawnow;
end

if (findAlignmentVerbose > 1)
  fprintf('calibration scaleAndTransformation guess = [%s], errVal = %g\n', ...
    sprintf('%g, ', scaleAndTransformation), errVal);
end

end


% function [angle] = minimizedAngle(angle)
% while angle < -pi
%   angle = angle + 2*pi;
% end
% while angle >= pi
%   angle = angle - 2*pi;
% end
% end



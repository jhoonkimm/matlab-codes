function [mobileCartMotion] = integrateMobileCartMotion(mobileCartRecord, varargin)

metersPerEncoderCount = 2.1494e-04;
% updateRate = 100;
forceRecalculateLeftEncoderFromRightEncoder = 1;
for i = 1:2:length(varargin)
  if (strcmp('metersPerEncoderCount', varargin{i}))
    metersPerEncoderCount = varargin{i + 1};
  end
  if (strcmp('forceRecalculateLeftEncoderFromRightEncoder', varargin{i}))
    forceRecalculateLeftEncoderFromRightEncoder = varargin{i + 1};
  end
  %   if (strcmp('updateRate', varargin{i}))
  %     updateRate = varargin{i + 1};
  %   end
end


%% This code allowed us to calculate the wheel base distance in encoder
% ticks, when looking at a trial with both encoders working (and some
% rotation). We use this wheel base to generate encoders values when one of
% the encoders is bad, using the other encoder and the gyro.
inspectWheelBaseDistanceCalculation = 0;
if (inspectWheelBaseDistanceCalculation)
  cla;
  plot(mobileCartRecord.leftEncoder, 'r')
  hold on
  plot(mobileCartRecord.rightEncoder, 'g')
  wheelBaseInEncoderTicks = 2.240897293288014e+003; %minimizingWheelBaseInEncoderTicks; %2000;
  simulatedRightEncoder = -mobileCartRecord.leftEncoder + wheelBaseInEncoderTicks * mobileCartRecord.gyroIntegratedRadians;
  plot(simulatedRightEncoder, 'k')
  norm(-mobileCartRecord.leftEncoder + wheelBaseInEncoderTicks * mobileCartRecord.gyroIntegratedRadians - ...
    mobileCartRecord.rightEncoder)
  
  errFunction = @(wheelBaseInEncoderTicks) ...
    norm(-mobileCartRecord.leftEncoder + wheelBaseInEncoderTicks * mobileCartRecord.gyroIntegratedRadians - ...
    mobileCartRecord.rightEncoder);
  [minimizingWheelBaseInEncoderTicks errVal] = fminunc(errFunction, 2000)
end


%% fix angle wraparaounds
% angles = mobileCartRecord.gyroIntegratedRadians(1:end-1); % - mobileCartRecord.gyroIntegratedRadians(1);
angles = mobileCartRecord.gyroIntegratedRadians(1:end); % - mobileCartRecord.gyroIntegratedRadians(1);

% this removes rotations of greater than 3 radians in a single tick
% (wraparound)... introduces ever so slight an error each time it happens.
angleDiffs = diff(angles);
angleDiffs(abs(angleDiffs) > 3) = 0;
anglesFixed = cumsum(angleDiffs);
angles = [0; anglesFixed];

%% check for and fix dead encoder:
wheelBaseInEncoderTicks = 2.240897293288014e+003; %minimizingWheelBaseInEncoderTicks; %2000;
if (sum(abs(diff(mobileCartRecord.leftEncoder))) == 0)
  fprintf('integrateMobileCartMotion: warning, left encoder appears to be dead! Calculating it based on right encoder and gyro\n');
  %   mobileCartRecord.leftEncoder = mobileCartRecord.rightEncoder;
  mobileCartRecord.leftEncoder = -mobileCartRecord.rightEncoder - wheelBaseInEncoderTicks *angles;
end
if (sum(abs(diff(mobileCartRecord.rightEncoder))) == 0)
  fprintf('integrateMobileCartMotion: warning, right encoder appears to be dead! Calculating it based on left encoder and gyro\n');
  %   mobileCartRecord.rightEncoder = mobileCartRecord.leftEncoder;
  mobileCartRecord.rightEncoder = -mobileCartRecord.leftEncoder + wheelBaseInEncoderTicks * angles;
end

if (forceRecalculateLeftEncoderFromRightEncoder)
  %   fprintf('integrateMobileCartMotion: warning, the code does not trust the left encoder at all, we are recalculating it based on the right encoder and gyro\n');
  %   mobileCartRecord.leftEncoder = mobileCartRecord.rightEncoder;
  mobileCartRecord.leftEncoder = -mobileCartRecord.rightEncoder - wheelBaseInEncoderTicks *angles;
end


angles = angles(1:end-1);


% attempt to check the error... unfortunately, wraparound screws us
% again:
% endAngle = angles(end) + mobileCartRecord.gyroIntegratedRadians(1);
% while (endAngle > pi)
%   endAngle = endAngle - 1 * pi;
% end
% while (endAngle < -pi)
%   endAngle = endAngle + 1 * pi;
% end
% fprintf('fixed angles diff = %g\n', endAngle - mobileCartRecord.gyroIntegratedRadians(end));

displacementDirections = [cos(angles) sin(angles)];

% angleFromEncoders = (processedPhasespaceTrial.cartSensors.leftEncoder + processedPhasespaceTrial.cartSensors.rightEncoder);
% angleFromEncoders = angleFromEncoders * (max(rads) / max(abs(angleFromEncoders)));
% angleFromEncoders(end) = [];

forwardEncoderMotion = (mobileCartRecord.leftEncoder - mobileCartRecord.rightEncoder) / 2;
displacementMagnitudes = -diff(forwardEncoderMotion);
displacementsPacked = repmat(displacementMagnitudes, [1 2]) .* displacementDirections;
displacements = [displacementsPacked(:,1) displacementsPacked(:,2)];
phasespaceCartLocation = cumsum(displacements);

mobileCartMotion.positions = phasespaceCartLocation * metersPerEncoderCount;
mobileCartMotion.angles = angles;

%% calculate sample times from stellaris clock
stellarisClockDiff = diff(mobileCartRecord.timer);
stellarisClockDiffRolloverInds = find(stellarisClockDiff > 0);
stellarisClockDiff(stellarisClockDiffRolloverInds) = stellarisClockDiff(stellarisClockDiffRolloverInds) - 8e6;
cartTimes = [0; cumsum(-stellarisClockDiff / 8e6)];
cartTimes(end) = [];
% cartTimes(end) = [];


% mobileCartMotion.updateRate = updateRate;
% mobileCartMotion.times = (0:(length(angles)-1)) * 1/updateRate;

mobileCartMotion.times = cartTimes;

%% ok, this is to check if perhaps I am a moron. Previous tests have been conclusive, but one more can't hurt:
% fixedTimes = 0:(1/100):cartTimes(end);
% fixedPositions = spline(cartTimes, mobileCartMotion.positions', fixedTimes)';
% fixedAngles = spline(cartTimes, mobileCartMotion.angles', fixedTimes)';
% mobileCartMotion.times = fixedTimes;
% mobileCartMotion.positions = fixedPositions;
% mobileCartMotion.angles = fixedAngles;


end


function [segmentedPeriods] = calculatePhasePeriodsFromRightLegPower(rightLegCOMPower, varargin)

%%

verbose = 0;

for i = 1:2:length(varargin)
  option = varargin{i};
  value = varargin{i+1};
  switch (option)
    case 'verbose'
      verbose = value;
  end
end

%%

numTicks = length(rightLegCOMPower);

rightPushoffPeriod = [0 0];
rightPreloadPeriod = [0 0];
rightReboundPeriod = [0 0];
rightCollisionPeriod = [0 0];
% times = linspace(0, timeDuration, numTicks);



for i = numTicks : -1 : 2
  previousPower = rightLegCOMPower(i - 1);
  thisPower = rightLegCOMPower(i);
  
  powerSwitchSignsThisTick = (thisPower * previousPower) < 0;
  
  if (thisPower >= 17 && rightPushoffPeriod(2) == 0)
    rightPushoffPeriod(2) = i;
    continue
  end
  
  if (powerSwitchSignsThisTick && ...
      rightPushoffPeriod(1) == 0 && ...
      rightPushoffPeriod(2) ~= 0)
    
    rightPushoffPeriod(1) = i;
    rightPreloadPeriod(2) = i-1;
    continue
  end
  
  if (powerSwitchSignsThisTick && ...
      rightPreloadPeriod(1) == 0 && ...
      rightPreloadPeriod(2) ~= 0)
    
    rightPreloadPeriod(1) = i;
    rightReboundPeriod(2) = i-1;
    continue
  end
  
  if (powerSwitchSignsThisTick && ...
      rightReboundPeriod(1) == 0 && ...
      rightReboundPeriod(2) ~= 0)
    
    rightReboundPeriod(1) = i;
    rightCollisionPeriod(2) = i-1;
    rightCollisionPeriod(1) = 1;
    continue
  end
  
end

segmentedPeriods.collision = rightCollisionPeriod;
segmentedPeriods.rebound = rightReboundPeriod;
segmentedPeriods.preload = rightPreloadPeriod;
segmentedPeriods.pushoff = rightPushoffPeriod;

if (verbose > 0)
  %   figure();
  plot(rightLegCOMPower);
  
  minVal = min(rightLegCOMPower);
  maxVal = max(rightLegCOMPower);
  
  for i = 1:2
    line(segmentedPeriods.collision(i) * ones(2, 1), ...
      [minVal maxVal], 'color', 'r');
    line(segmentedPeriods.rebound(i) * ones(2, 1), ...
      [minVal maxVal], 'color', 'g');
    line(segmentedPeriods.preload(i) * ones(2, 1), ...
      [minVal maxVal], 'color', 'b');
    line(segmentedPeriods.pushoff(i) * ones(2, 1), ...
      [minVal maxVal], 'color', 'c');
  end
end



% collision = trapz(timeR(RHS_period(1):RHS_period(2)), signalToSegment(RHS_period(1):RHS_period(2)));
% rebound = trapz(timeR(RRB_period(1):RRB_period(2)), signalToSegment(RRB_period(1):RRB_period(2)));
% preload = trapz(timeR(RPL_period(1):RPL_period(2)), signalToSegment(RPL_period(1):RPL_period(2)));
% pushoff = trapz(timeR(RPO_period(1):RPO_period(2)), signalToSegment(RPO_period(1):RPO_period(2)));


end
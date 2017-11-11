function [collisionWork reboundWork preloadWork pushoffWork] = calculatePhaseWorkFromPower(...
  power, phasePeriods, stepDuration)


numTicks = length(power);
times = linspace(0, stepDuration, numTicks);

collisionWork = trapz(times(phasePeriods.collision(1):phasePeriods.collision(2)), ...
  power(phasePeriods.collision(1):phasePeriods.collision(2)));

reboundWork = trapz(times(phasePeriods.rebound(1):phasePeriods.rebound(2)), ...
  power(phasePeriods.rebound(1):phasePeriods.rebound(2)));

preloadWork = trapz(times(phasePeriods.preload(1):phasePeriods.preload(2)), ...
  power(phasePeriods.preload(1):phasePeriods.preload(2)));

pushoffWork = trapz(times(phasePeriods.pushoff(1):phasePeriods.pushoff(2)), ...
  power(phasePeriods.pushoff(1):phasePeriods.pushoff(2)));


end


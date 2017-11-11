function [goodIndeces] = removeBadStridesFromRunningOnTreadmill(grf, HS, TO, varargin)

% removeBadStridesFromRunningOnTreadmill identify good vs bad treadmill strides
%  For analyzing ground reactions forces (GRFs) of a split-belt
%  treadill. Input right and left GRFs each as long arrays of entire trial.
%  And inputheelstrike and toeoff events for each foot, which can be found using
%  findEventIndicesForWalkingOnTreadmill.m.
%  Function returns good stride numbers by default.
%
%  Options: (note: all thresholds have default values set in code, bu can be overwritten with these options)
%    slopeThreshold: threshold used to determine if slope of loading/unloading
%       GRFs is close enough to median slope (0.3-0.35 seems to work fairly well)
%    badStrideTimeThreshold: threshold used to determine if stride time is
%       close enough to median stride time
%    badStepTimeThreshold: threshold used to determine if step time is
%       close enough to median steptime
%    outputValues: return good strides ('goodStrideNumbers') or bad strides
%       ('badStrideNumbers)

%  adapted from walking version

%% Thresholds
slopeThreshold = 0.35; % arbitrary threshold value can be tweaked, but 0.3-0.35 seems to work fairly well
badStrideTimeThreshold = 0.25; % percent diff from mean (e.g., 0.1 = 10%)
badStepTimeThreshold = 0.25;
outputValues = 'goodStrideNumbers'; % or 'badStrideNumbers'

%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'slopeThreshold'
      slopeThreshold = val;
    case 'badStrideTimeThreshold'
      badStrideTimeThreshold = val;
    case 'badStepTimeThreshold'
      badStepTimeThreshold = val;
    case 'outputValues'
      outputValues = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end



%% Initialize variables
badOnes = [];
goodOnes = [];

%% Constants
N = 1000; % this value can be any large number

%% Check number of left and right steps
numSteps = length(HS);

%% Zero Fz swing forces (only if necesary)
findZeros = find(grf(:,3) == 0);
if length(findZeros) / length(grf) < 0.3
  zeroThreshold = mean(grf(:,3)) * 0.02;
  grf(:, 3) = grf(:, 3).*(grf(:, 3) > zeroThreshold);
end

%% Mean Step & Stride times (really calculated in indices, but doesnt matter here b/c just need these for statistics)
stepTime = TO - HS;
avgSteptime = mean(stepTime);
strideTime = diff(HS);
avgStridetime = mean(strideTime);

avgStrideTime = avgStridetime;

%% Assign forces to individual steps and interpolate all steps to one stance phase
for i = 1 : length(HS)
  GRFzr(:, i) = interpGaitCycle(grf(HS(i) : TO(i), 3), N);
end

%% Assign forces to strides and interpolate to one stride
for i = 1 : length(HS) - 1
  goodOnes = cat(2, goodOnes, i); % initially assume all strides are good
  zForcesInStride(:, i) = interpGaitCycle(grf(HS(i) : HS(i+1)-1, 3));
  
  % Calculate swing time
  swingPhaseTime(i) = length(find(zForcesInStride(:,i) == 0));
end


%% Check for Bad Strides
% Median Fz loading curve
loadingSlopeIndices = [1:50]';
medianLoadingCurveGRFzr = median(GRFzr(loadingSlopeIndices,:),2);
medianLoadingSlopeGRFzr = polyfit(loadingSlopeIndices, medianLoadingCurveGRFzr,1);
medianLoadingSlopeGRFzr = medianLoadingSlopeGRFzr(1);

% Median Fz unloading curve
unloadingSlopeIndicesGRFzr = [length(GRFzr) - 49 : length(GRFzr)]';
medianUnloadingCurveGRFzr = median(GRFzr(unloadingSlopeIndicesGRFzr,:), 2);
medianUnloadingSlopeGRFzr = polyfit(unloadingSlopeIndicesGRFzr,medianUnloadingCurveGRFzr, 1);
medianUnloadingSlopeGRFzr = medianUnloadingSlopeGRFzr(1);

% Trial-by-trial loading & unloading curves
loadingCurveGRFzr = GRFzr(loadingSlopeIndices,:);
unloadingCurveGRFzr = GRFzr(unloadingSlopeIndicesGRFzr,:);

% Right
for i = 1 : (numSteps - 1)
  clear loadingSlopeTempL loadingSlopeTempR unloadingSlopeTempL unloadingSlopeTempR
  
  % Linear fits to early loading and late unloading
  %   loadingSlopeTempL = polyfit(loadingSlopeIndices,loadingCurveLGRFzl(:,i),1); loadingSlopeLGRFzl(i,1) = loadingSlopeTempL(1);
  loadingSlopeTemp = polyfit(loadingSlopeIndices, loadingCurveGRFzr(:,i),1);
  loadingSlopeGRFzr(i,1) = loadingSlopeTemp(1);
  %   unloadingSlopeTempL = polyfit(unloadingSlopeIndicesLGRFzl,unloadingCurveLGRFzl(:,i),1); unloadingSlopeLGRFzl(i,1) = unloadingSlopeTempL(1);
  unloadingSlopeTemp = polyfit(unloadingSlopeIndicesGRFzr, unloadingCurveGRFzr(:,i),1);
  unloadingSlopeGRFzr(i,1) = unloadingSlopeTemp(1);
  
  % For right strides, includes right side loading/unloading, left side loading and left side unloading from previous (i-1) step
  %   if i==1 | ... % remove first stride b/c we are missing preceding left step
  %       diff([loadingSlopeLGRFzl(i) medianLoadingSlopeLGRFzl]) > slopeThreshold*medianLoadingSlopeLGRFzl | ... % remove if loading rate is too slow
  %       diff([loadingSlopeGRFzr(i) medianLoadingSlopeGRFzr]) > slopeThreshold*medianLoadingSlopeGRFzr | ...
  %       diff([unloadingSlopeLGRFzl(i-1) medianUnloadingSlopeLGRFzl]) < slopeThreshold*medianUnloadingSlopeLGRFzl | ... % remove if unloading rate is too slow
  %       diff([unloadingSlopeGRFzr(i) medianUnloadingSlopeGRFzr]) < slopeThreshold*medianUnloadingSlopeGRFzr | ... % note: more negative slope means steeper/faster rate of unloading
  %       abs((Rstridetime(i)-avgStrideTime)) > badStrideTimeThreshold*mean(Rstridetime) | ... % stride time too different
  %       abs((steptime(i)-avgSteptime)) > badStepTimeThreshold*mean(steptime) | ... % Right step time too different
  %       abs((Lsteptime(i-1)-avgLsteptime)) > badStepTimeThreshold*mean(Lsteptime) | ... % Previous Left step time too different
  %       abs((swingPhaseTimeR(i)-median(swingPhaseTimeR))) > badStrideTimeThreshold*median(swingPhaseTimeR) | ... % right swing phase time too different
  %       abs((LswingPhaseTimeR(i)-median(LswingPhaseTimeR))) > badStrideTimeThreshold*median(LswingPhaseTimeR) % left swing phase time too different
  %
  
  if i == 1 || ... % remove first stride b/c we are missing preceding left step
      diff([loadingSlopeGRFzr(i) medianLoadingSlopeGRFzr]) > slopeThreshold * medianLoadingSlopeGRFzr || ...
      diff([unloadingSlopeGRFzr(i) medianUnloadingSlopeGRFzr]) < slopeThreshold * medianUnloadingSlopeGRFzr || ... % note: more negative slope means steeper/faster rate of unloading
      abs((strideTime(i) - avgStrideTime)) > badStrideTimeThreshold * mean(strideTime) || ... % stride time too different
      abs((stepTime(i) - avgSteptime)) > badStepTimeThreshold * mean(stepTime) || ... % Right step time too different
      abs((swingPhaseTime(i) - median(swingPhaseTime))) > badStrideTimeThreshold * median(swingPhaseTime) % right swing phase time too different
    
    goodOnes = setdiff(goodOnes, i); % if true, then remove this stride from final averaging
    badOnes = cat(2,badOnes,i);
  end
  
end



%% Plot results
% figure(100066); clf; plot([RGRFzR(:,Rbadones(2:end)) LGRFzR(:,Rbadones(2:end))]); title('bad strides') % not plotting 1 b/c we dont have preceding LTO
% figure(100100); clf; plot([RGRFzR(:,Rgoodones) LGRFzR(:,Rgoodones)]); title('good strides')
%
% figure(100067); clf; plot([LGRFzL(:,Lbadones) RGRFzL(:,Lbadones)]); title('bad strides') % not plotting 1 b/c we dont have preceding LTO
% figure(100101); clf; plot([LGRFzL(:,Lgoodones) RGRFzL(:,Lgoodones)]); title('good strides')



%% Return values

goodIndeces = goodOnes;

% if strcmp(outputValues,'goodStrideNumbers')
%   out1=goodOnes; out2=Lgoodones;
% elseif strcmp(outputValues,'badStrideNumbers')
%   out1=badOnes; out2=Lbadones;
% else
%   out1=[]; out2=[];
%   error('invalid outputValues')
% end


function [out1, out2] = removeBadStridesFromWalkingOnTreadmill(Rgrf,Lgrf,RHS,LHS,RTO,LTO,varargin)

% REMOVEBADSTRIDESFROMWALKINGONTREADMILL identify good vs bad treadmill strides
%  For analyzing ground reactions forces (GRFs) of a split-belt
%  treadill. Input right and left GRFs each as long arrays of entire trial.
%  And inputheelstrike and toeoff events for each foot, which can be found using
%  findEventIndicesForWalkingOnTreadmill.m.
%  Function returns good stride numbers by default.
%
%  removeBadStridesFromWalkingOnTreadmill(Rgrf,Lgrf,RHS,LHS,RTO,LTO,varargin)
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

%  Karl Zelik
%  10/30/09

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
Rbadones = []; Lbadones = []; % bad strides
Rgoodones = []; Lgoodones = []; % good strides

%% Constants
N=1000; % this value can be any large number

%% Check number of left and right steps
numRsteps = length(RHS);
numLsteps = length(LHS);
if numRsteps ~= numLsteps
    error('\nError: Step Mismatch: # of left steps does not agree with # of right steps \n');
end


%% Zero Fz swing forces (only if necesary)
findZerosR = find(Rgrf(:,3)==0);
if length(findZerosR)/length(Rgrf) < 0.3
    zeroThresholdR = mean(Rgrf(:,3))*0.02;
    Rgrf(:,3) = Rgrf(:,3).*(Rgrf(:,3)>zeroThresholdR);
end
findZerosL = find(Lgrf(:,3)==0);
if length(findZerosL)/length(Lgrf) < 0.3
    zeroThresholdL = mean(Lgrf(:,3))*0.02;
    Lgrf(:,3) = Lgrf(:,3).*(Lgrf(:,3)>zeroThresholdL);
end

%% Mean Step & Stride times (really calculated in indices, but doesnt matter here b/c just need these for statistics)
Rsteptime = (RTO-RHS);
avgRsteptime = mean(Rsteptime);
Lsteptime = (LTO-LHS);
avgLsteptime = mean(Lsteptime);
Rstridetime = (diff(RHS));
avgRstridetime = mean(Rstridetime); 
Lstridetime = (diff(LHS));
avgLstridetime = mean(Lstridetime);
avgstridetime = avgRstridetime; % RL vs LR should not effect stride time

%% Assign forces to individual steps and interpolate all steps to one stance phase
% RIGHT STEPS
for i=1:length(RHS) 
    RGRFzr(:,i) = interpGaitCycle(Rgrf(RHS(i):RTO(i),3), N);
end
% LEFT STEPS
for i=1:length(LHS)
    LGRFzl(:,i) = interpGaitCycle(Lgrf(LHS(i):LTO(i),3), N);
end


%% Assign forces to strides and interpolate to one stride
% RIGHT HS to HS
for i=1:length(RHS)-1
    Rgoodones = cat(2,Rgoodones,i); % initially assume all strides are good    
    RGRFzR(:,i) = interpGaitCycle(Rgrf(RHS(i):RHS(i+1)-1,3)); 
    LGRFzR(:,i) = interpGaitCycle(Lgrf(RHS(i):RHS(i+1)-1,3));
    % Calculate swing time
    RswingPhaseTimeR(i) = length(find(RGRFzR(:,i) == 0));
    LswingPhaseTimeR(i) = length(find(LGRFzR(:,i) == 0));   
end
% LEFT HS to HS
for i=1:length(LHS)-1
    Lgoodones = cat(2,Lgoodones,i); % initially assume all strides are good    
    LGRFzL(:,i) = interpGaitCycle(Lgrf(LHS(i):LHS(i+1)-1,3));    
    RGRFzL(:,i) = interpGaitCycle(Rgrf(LHS(i):LHS(i+1)-1,3));
    % Calculate swing time
    LswingPhaseTimeL(i) = length(find(LGRFzL(:,i) == 0));
    RswingPhaseTimeL(i) = length(find(RGRFzL(:,i) == 0));
end


%% Check for Bad Strides
% Median Fz loading curve
loadingSlopeIndices = [1:50]';
medianLoadingCurveLGRFzl = median(LGRFzl(loadingSlopeIndices,:),2);
medianLoadingCurveRGRFzr = median(RGRFzr(loadingSlopeIndices,:),2);
medianLoadingSlopeLGRFzl = polyfit(loadingSlopeIndices,medianLoadingCurveLGRFzl,1); medianLoadingSlopeLGRFzl = medianLoadingSlopeLGRFzl(1);
medianLoadingSlopeRGRFzr = polyfit(loadingSlopeIndices,medianLoadingCurveRGRFzr,1); medianLoadingSlopeRGRFzr = medianLoadingSlopeRGRFzr(1);

% Median Fz unloading curve
unloadingSlopeIndicesLGRFzl = [length(LGRFzl)-49:length(LGRFzl)]';
unloadingSlopeIndicesRGRFzr = [length(RGRFzr)-49:length(RGRFzr)]';
medianUnloadingCurveLGRFzl = median(LGRFzl(unloadingSlopeIndicesLGRFzl,:),2);
medianUnloadingCurveRGRFzr = median(RGRFzr(unloadingSlopeIndicesRGRFzr,:),2);
medianUnloadingSlopeLGRFzl = polyfit(unloadingSlopeIndicesLGRFzl,medianUnloadingCurveLGRFzl,1); medianUnloadingSlopeLGRFzl = medianUnloadingSlopeLGRFzl(1);
medianUnloadingSlopeRGRFzr = polyfit(unloadingSlopeIndicesRGRFzr,medianUnloadingCurveRGRFzr,1); medianUnloadingSlopeRGRFzr = medianUnloadingSlopeRGRFzr(1);

% Trial-by-trial loading & unloading curves
loadingCurveLGRFzl = LGRFzl(loadingSlopeIndices,:);
loadingCurveRGRFzr = RGRFzr(loadingSlopeIndices,:);
unloadingCurveLGRFzl = LGRFzl(unloadingSlopeIndicesLGRFzl,:);
unloadingCurveRGRFzr = RGRFzr(unloadingSlopeIndicesRGRFzr,:);


% Right
for i=1:numRsteps-1
    clear loadingSlopeTempL loadingSlopeTempR unloadingSlopeTempL unloadingSlopeTempR
    
    % Linear fits to early loading and late unloading
    loadingSlopeTempL = polyfit(loadingSlopeIndices,loadingCurveLGRFzl(:,i),1); loadingSlopeLGRFzl(i,1) = loadingSlopeTempL(1);
    loadingSlopeTempR = polyfit(loadingSlopeIndices,loadingCurveRGRFzr(:,i),1); loadingSlopeRGRFzr(i,1) = loadingSlopeTempR(1);
    unloadingSlopeTempL = polyfit(unloadingSlopeIndicesLGRFzl,unloadingCurveLGRFzl(:,i),1); unloadingSlopeLGRFzl(i,1) = unloadingSlopeTempL(1);
    unloadingSlopeTempR = polyfit(unloadingSlopeIndicesRGRFzr,unloadingCurveRGRFzr(:,i),1); unloadingSlopeRGRFzr(i,1) = unloadingSlopeTempR(1);
    
    % For right strides, includes right side loading/unloading, left side loading and left side unloading from previous (i-1) step
    if i==1 | ... % remove first stride b/c we are missing preceding left step
            diff([loadingSlopeLGRFzl(i) medianLoadingSlopeLGRFzl]) > slopeThreshold*medianLoadingSlopeLGRFzl | ... % remove if loading rate is too slow
            diff([loadingSlopeRGRFzr(i) medianLoadingSlopeRGRFzr]) > slopeThreshold*medianLoadingSlopeRGRFzr | ...
            diff([unloadingSlopeLGRFzl(i-1) medianUnloadingSlopeLGRFzl]) < slopeThreshold*medianUnloadingSlopeLGRFzl | ... % remove if unloading rate is too slow
            diff([unloadingSlopeRGRFzr(i) medianUnloadingSlopeRGRFzr]) < slopeThreshold*medianUnloadingSlopeRGRFzr | ... % note: more negative slope means steeper/faster rate of unloading
            abs((Rstridetime(i)-avgstridetime)) > badStrideTimeThreshold*mean(Rstridetime) | ... % stride time too different
            abs((Rsteptime(i)-avgRsteptime)) > badStepTimeThreshold*mean(Rsteptime) | ... % Right step time too different
            abs((Lsteptime(i-1)-avgLsteptime)) > badStepTimeThreshold*mean(Lsteptime) | ... % Previous Left step time too different
            abs((RswingPhaseTimeR(i)-median(RswingPhaseTimeR))) > badStrideTimeThreshold*median(RswingPhaseTimeR) | ... % right swing phase time too different
            abs((LswingPhaseTimeR(i)-median(LswingPhaseTimeR))) > badStrideTimeThreshold*median(LswingPhaseTimeR) % left swing phase time too different
        
        Rgoodones = setdiff(Rgoodones, i); % if true, then remove this stride from final averaging
        Rbadones = cat(2,Rbadones,i);
    end
    
end

% Left
for i=1:numLsteps-1    
    % For left strides, includes left side loading/unloading, right side unloading and right side unloading from subsequent (i+1) step
    if i==numLsteps-1 | ... % get ride of last stride b/c do not have susequent right stride data
            diff([loadingSlopeLGRFzl(i) medianLoadingSlopeLGRFzl]) > slopeThreshold*medianLoadingSlopeLGRFzl | ... % remove if loading rate is too slow
            diff([loadingSlopeRGRFzr(i+1) medianLoadingSlopeRGRFzr]) > slopeThreshold*medianLoadingSlopeRGRFzr | ...
            diff([unloadingSlopeLGRFzl(i) medianUnloadingSlopeLGRFzl]) < slopeThreshold*medianUnloadingSlopeLGRFzl | ... % remove if unloading rate is too slow
            diff([unloadingSlopeRGRFzr(i) medianUnloadingSlopeRGRFzr]) < slopeThreshold*medianUnloadingSlopeRGRFzr | ... % note: more negative slope means steeper/faster rate of unloading
            abs((Lstridetime(i)-avgstridetime)) > badStrideTimeThreshold*mean(Lstridetime) | ... % stride time too different
            abs((Rsteptime(i)-avgRsteptime)) > badStepTimeThreshold*mean(Rsteptime) | ... % Right step time too different
            abs((Lsteptime(i)-avgLsteptime)) > badStepTimeThreshold*mean(Lsteptime) | ...% Left step time too different
            abs((RswingPhaseTimeL(i)-median(RswingPhaseTimeL))) > badStrideTimeThreshold*median(RswingPhaseTimeL) | ... % right swing phase time too different
            abs((LswingPhaseTimeL(i)-median(LswingPhaseTimeL))) > badStrideTimeThreshold*median(LswingPhaseTimeL) % left swing phase time too different
        
        Lgoodones = setdiff(Lgoodones, i); % if true, then remove this stride from final averaging
        Lbadones = cat(2,Lbadones,i);
    end
    
end


%% Plot results
% figure(100066); clf; plot([RGRFzR(:,Rbadones(2:end)) LGRFzR(:,Rbadones(2:end))]); title('bad strides') % not plotting 1 b/c we dont have preceding LTO 
% figure(100100); clf; plot([RGRFzR(:,Rgoodones) LGRFzR(:,Rgoodones)]); title('good strides')
% 
% figure(100067); clf; plot([LGRFzL(:,Lbadones) RGRFzL(:,Lbadones)]); title('bad strides') % not plotting 1 b/c we dont have preceding LTO 
% figure(100101); clf; plot([LGRFzL(:,Lgoodones) RGRFzL(:,Lgoodones)]); title('good strides')




%% Return values
if strcmp(outputValues,'goodStrideNumbers')
    out1=Rgoodones; out2=Lgoodones;
elseif strcmp(outputValues,'badStrideNumbers')
    out1=Rbadones; out2=Lbadones;    
else
    out1=[]; out2=[];
    error('invalid outputValues')
end


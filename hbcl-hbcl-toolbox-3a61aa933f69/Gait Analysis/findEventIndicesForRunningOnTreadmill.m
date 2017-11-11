function [HS, TO] = findEventIndicesForRunningOnTreadmill(grf, varargin)

% findEventIndicesForRunningOnTreadmill find HS and TO events
%  
%  [HS, TO] = findEventIndicesForRunningOnTreadmill(grf)
%  Assumes sampling rate of 120 Hz.
%
%  [HS, TO] = findEventIndicesForRunningOnTreadmill(grf, varargin)
%  Uses sampling rate from VARARGIN input.
%
%  Summary of how it works:
%  Looks at vertical force (Fz). Uses a threshold and slope approach
%  to find loading and unloading regions of Fz. Searches backwards down the
%  loading curve (within a fixed window) to find the nearest zero or local minimum --
%  heelstrike. Then searches forward down the unloading curve (within a fixed
%  window) to find nearest zero or local minimum -- toeoff. Fz threshold is
%  incrementally adjusted if method is unsuccessful on first pass.

% modified from findEvent..WalkingOnTreadmill.m

if nargin>=3
  samplingfreq = varargin{1}; % Hz
else
  samplingfreq = 120; % Hz
end


showPlots = 0; % display plots (debug mode)
figNum = 25000; % plot figure number
complete = 0; % flag to indicate when function is complete

%% Quasi-constants
% (1)
zcoeff0 = 0.6;
zcoeff=zcoeff0; % initial coefficient used to define Fz threshold --> mean(Fz)*coeff = Fz threshold used in algorithm

% (2)
windowSize = round(0.8*samplingfreq); % define search window size in which to look for events (HS, TO)
n = 7; % hack constant used to search for local minima

% (3)
wn = 0.2;
originalGRF = grf;

while ~complete
  
  %% Initialize variables each time thru loop
  HS = [];
  TO = [];
  nearHS = [];
  nearTO = [];
  
  % Define/update Fz thresholds
  zThreshold = zcoeff * mean(grf(:,3));
  
  %% (1) Find where data crosses Fz threshold -- these data points are near events (HS, TO) of interest
  for i = 2 : length(grf)-1
    % heelstrikes
    if (grf(i,3) < zThreshold && grf(i+1,3) > zThreshold)
      nearHS = [nearHS; i+1];
    end
    % toe-offs
    if (grf(i,3) > zThreshold && grf(i+1, 3) < zThreshold)
      nearTO = [nearTO; i-1];
    end
  end
  
  %% (2) Search in window to find precise heelstrike and toeoff events
  for j = 2 : max(find(nearHS + windowSize < length(grf)))
    for i = nearHS(j) : -1 : nearHS(j) - windowSize
      % find closest zero, local minimum or end of window:
      if grf(i,3) == 0 | grf(i,3) < grf(i-n, 3) | i == nearHS(j) - windowSize
        HS = [HS; i];
        break;
      end
    end
  end
  
  for j = 2 : max(find(nearTO + windowSize < length(grf)))
    for i = nearTO(j) : 1 : nearTO(j) + windowSize
      if grf(i,3) == 0 | grf(i,3) < grf(i+n, 3) | i == nearTO(j) + windowSize
        TO = [TO; i];
        break;
      end
    end
  end
  
  %% (3) Remove events if they are not unique -- can ocur when search (2) hits end of window
  HS = unique(HS);
  TO = unique(TO);
  
  %% (4) Trim strides so that we always start with a RHS and end with LTO
  % remove extra toe-off if before first heelstrike
  if TO(1) < HS(1)
    TO = TO(2:end);
  end
  
  % remove extra heelstrike at end if toe-off is not in region
  if HS(end) > TO(end)
    HS = HS(1:end-1);
  end
  
  %% (5) Plot
  %   if showPlots
  %     figure(figNum); clf;
  %     subplot(211); plot(Rgrf(:,3)); hold on; plot(RHS, Rgrf(RHS,3),'ro','MarkerSize',6); hold on; plot(RTO, Rgrf(RTO,3), 'cs','MarkerSize',6); title('Rgrf')
  %     hold on; plot(nearRHS, Rgrf(nearRHS,3),'r*','MarkerSize',3); hold on; plot(nearRTO, Rgrf(nearRTO,3),'c*','MarkerSize',3)
  %     subplot(212); plot(Lgrf(:,3),'k'); hold on; plot(LHS, Lgrf(LHS,3),'ro','MarkerSize',6); hold on; plot(LTO, Lgrf(LTO,3), 'cs','MarkerSize',6); title('Lgrf')
  %     hold on; plot(nearLHS, Lgrf(nearLHS,3),'r*','MarkerSize',3); hold on; plot(nearLTO, Lgrf(nearLTO,3),'c*','MarkerSize',3)
  %   end
  
  
  %% (6) Check Progress
  if zcoeff <= 0.1 % If z threshold gets too low, try filterng and repeating adaptive threshold searching
    if wn >= 0.1 % try filtering
      wn = wn - 0.02;
      zcoeff = zcoeff0;
      [B,A] = butter(3, wn);
      grf = filtfilt(B, A, originalGRF);
      %       Lgrf = filtfilt(B,A,originalLgrf);
    else % give up if threshold gets too low and stiff have not found walking events
      complete = 1;
      HS=[];
      TO=[];
      fprintf('\nWarning: Unable to find walking events\n')
    end
    
  elseif (length(HS) == length(TO)) % Check that number of events found are consistent
    complete = 1;
  else
    zcoeff = zcoeff - 0.02; % incrementally lower Fz threshold and try again
  end
  
end


% close(figNum)


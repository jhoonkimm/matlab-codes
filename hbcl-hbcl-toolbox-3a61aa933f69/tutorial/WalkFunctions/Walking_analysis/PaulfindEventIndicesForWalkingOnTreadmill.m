function [RHS,LHS,RTO,LTO] = PaulfindEventIndicesForWalkingOnTreadmill(Rgrf, Lgrf, varargin)

% FINDEVENTINDICESFORWALKINGONTREADMILL find HS and TO events
%  Find heelstrike (HS) and toeoff (TO) events from walking on
%  a split-belt force treadmill. Input right and left ground
%  reaction forces in Nx3 array and sampling rate (optional). 
%  Outputs indices for right HS, left HS, right TO and left TO. 
%  If method is unsuccessful it will print a warning to the 
%  command line and output an empty-set.
%
%  [RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(RGRF, LGRF)
%  Assumes sampling rate of 120 Hz.
%
%  [RHS,LHS,RTO,LTO] = findEventIndicesForWalkingOnTreadmill(RGRF, LGRF, VARARGIN)
%  Uses sampling rate from VARARGIN input.
% 
%  Summary of how it works:
%  Looks at vertical force (Fz). Uses a threshold and slope approach
%  to find loading and unloading regions of Fz. Searches backwards down the 
%  loading curve (within a fixed window) to find the nearest zero or local minimum -- 
%  heelstrike. Then searches forward down the unloading curve (within a fixed
%  window) to find nearest zero or local minimum -- toeoff. Fz threshold is 
%  incrementally adjusted if method is unsuccessful on first pass.

%  Karl Zelik
%  updated 12/17/09

if nargin>=3
    samplingfreq = varargin{1}; % Hz
else
    samplingfreq = 120; % Hz
end


showPlots = 1; % display plots (debug mode)
figNum = 25000; % plot figure number
complete = 0; % flag to indictae when function is complete

%% Quasi-constants
% (1) 
zcoeff0 = 0.2; zcoeff=zcoeff0; % initial coefficient used to define Fz threshold --> mean(Fz)*coeff = Fz threshold used in algorithm
 
% (2) 
windowSize = round(0.8*samplingfreq); % define search window size in which to look for events (HS, TO)
n=7; % hack constant used to search for local minima

% (3)
wn=0.2;
originalRgrf = Rgrf; originalLgrf = Lgrf;

while ~complete

    %% Initialize variables each time thru loop
    RHS = []; LHS = []; RTO = []; LTO = [];
    nearRHS = []; nearLHS = []; nearRTO = []; nearLTO = [];

     % Define/update Fz thresholds
    zRthreshold = zcoeff*mean(Rgrf(:,3));
    zLthreshold = zcoeff*mean(Lgrf(:,3));

    %% (1) Find where data crosses Fz threshold -- these data points are near events (HS, TO) of interest
    for i=2:length(Rgrf)-1
        % Right heelstrikes
        if (Rgrf(i,3)<zRthreshold && Rgrf(i+1,3)>zRthreshold)
            nearRHS = [nearRHS; i+1];
        end
        % Left heelstrikes
        if (Lgrf(i,3)<zLthreshold && Lgrf(i+1,3)>zLthreshold)
            nearLHS = [nearLHS; i+1];
        end

        % Right toe-offs
        if (Rgrf(i,3)>zRthreshold && Rgrf(i+1,3)<zRthreshold)
            nearRTO = [nearRTO; i-1];
        end
        % Left toe-offs
        if (Lgrf(i,3)>zLthreshold && Lgrf(i+1,3)<zLthreshold)
            nearLTO = [nearLTO; i-1];
        end
    end



    %% (2) Search in window to find precise heelstrike and toeoff events
    for j=3:max(find(nearRHS+windowSize<length(Rgrf))) % length(nearRHS)-1
        for i = nearRHS(j):-1:nearRHS(j)-windowSize
            if Rgrf(i,3)==0 & Rgrf(i-n,3)==0%| i==nearRHS(j)-windowSize % find closest zero, local minimum or end of window
                RHS = [RHS; i];
                break;
            end
        end
    end

    for j=3:max(find(nearLHS+windowSize<length(Rgrf))) % length(nearLHS)-1
        for i = nearLHS(j):-1:nearLHS(j)-windowSize
            if Lgrf(i,3)==0 & Lgrf(i-n,3)==0% | i==nearLHS(j)-windowSize
                LHS = [LHS; i];
                break;
            end
        end
    end

    for j=3:max(find(nearRTO+windowSize<length(Rgrf))) % length(nearRTO)-1
        for i = nearRTO(j):1:nearRTO(j)+windowSize
            if Rgrf(i,3)==0 &Rgrf(i+n,3)==0% | i==nearRTO(j)+windowSize
                RTO = [RTO; i];
                break;
            end
        end
    end

    for j=3:max(find(nearLTO+windowSize<length(Rgrf))) % length(nearLTO)-1
        for i = nearLTO(j):1:nearLTO(j)+windowSize
            if Lgrf(i,3)==0 &Lgrf(i+n,3)==0 %| i==nearLTO(j)+windowSize
                LTO = [LTO; i];
                break;
            end
        end
    end


    %% (3) Remove events if they are not unique -- can ocur when search (2) hits end of window
    RHS = unique(RHS); RTO = unique(RTO);
    LHS = unique(LHS); LTO = unique(LTO);


    %% (4) Trim strides so that we always start with a RHS and end with LTO
    % remove extra toe-off if before first heelstrike
    if RTO(1) < RHS(1)
        RTO = RTO(2:end);
    end
    if LTO(1) < LHS(1)
        LTO = LTO(2:end);
    end

    % remove extra heelstrike at end if toe-off is not in region
    if RHS(end) > RTO(end)
        RHS = RHS(1:end-1);
    end
    if LHS(end) > LTO(end)
        LHS = LHS(1:end-1);
    end

    % arbitrarily always start with a right heelstrike, end with left toe off
    if LHS(1) < RHS(1)
        LHS = LHS(2:end);
        LTO = LTO(2:end);
    end

    if RTO(end) > LTO(end)
        RTO = RTO(1:end-1);
        RHS = RHS(1:end-1);
    end


    %% (5) Plot
    if showPlots
        figure(figNum); clf;
        subplot(211); plot(Rgrf(:,3)); hold on; plot(RHS, Rgrf(RHS,3),'ro','MarkerSize',6); hold on; plot(RTO, Rgrf(RTO,3), 'cs','MarkerSize',6); title('Rgrf')
        hold on; plot(nearRHS, Rgrf(nearRHS,3),'r*','MarkerSize',3); hold on; plot(nearRTO, Rgrf(nearRTO,3),'c*','MarkerSize',3)
        subplot(212); plot(Lgrf(:,3),'k'); hold on; plot(LHS, Lgrf(LHS,3),'ro','MarkerSize',6); hold on; plot(LTO, Lgrf(LTO,3), 'cs','MarkerSize',6); title('Lgrf')
        hold on; plot(nearLHS, Lgrf(nearLHS,3),'r*','MarkerSize',3); hold on; plot(nearLTO, Lgrf(nearLTO,3),'c*','MarkerSize',3)
    end
      
    
    %% (6) Check Progress
    if zcoeff<=0.1 % If z threshold gets too low, try filterng and repeating adaptie threshold searching
        if wn >= 0.1 % try filtering
            wn = wn-0.02;
            zcoeff=zcoeff0;
            [B,A] = butter(3,wn);
            Rgrf = filtfilt(B,A,originalRgrf); Lgrf = filtfilt(B,A,originalLgrf);
        else % give up if threshold gets too low and stiff have not found walking events
            complete = 1;
            RHS=[]; LHS=[]; RTO=[]; LTO=[];
            fprintf('\nWarning: Unable to find walking events\n')
        end
        
    elseif (length(RHS)==length(LHS) && length(RTO)==length(LTO) && length(RHS)==length(RTO)) % Check that number of events found are consistent
        complete = 1;
    else
        zcoeff = zcoeff - 0.001; % incrementally lower Fz threshold and try again  
    end
    
end


% close(figNum)


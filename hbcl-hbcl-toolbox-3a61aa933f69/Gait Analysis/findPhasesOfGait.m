function [indices] = findPhasesOfGait(var,varargin)

% FINDPHASESOFGAIT find indices for gait phases based on CoM work rate
%  Find start and end indices for 4 phases of gait from individual limbs method
%  center of mass (CoM) work rate plot.
% 
%  Input VAR must contain COM work rate data for the full stride (not just stance 
%  phase). The right limb work rate must be plotted from right heelstrike 
%  to right heelstrike (and vice versa for the left side).
% 
%  Output contains indices for heelstrike (HS) and the end of each of the 4 CoM 
%  phases of gait -- collision (CO), rebound (RB), preload (PL), pushoff (PO). 
%  If this method is unable to find phases of gait from original work rate data, 
%  it use filtering to try to locate phase transitions. If iterative filtering 
%  is still unsuccessful then the function will output an empty-set 
%  and print a failure warning in the command window.
%
%  Set VARARGIN==1 will enable debug mode which will show plots of phase 
%  transitions and filtering (if necessary)
%  
%  Examples:
%   array_of_HS_COend_RBend_PLend_POend = findPhasesOfGait(Right_Leg_Individual_Limbs_Method_CoM_Work_Rate );
%   findPhasesOfGait(Right_Leg_Individual_Limbs_Method_CoM_Work_Rate, 1); % show debug plots

%  Karl Zelik
%  updated 10/29/09


complete = 0; % flag to indicate function is finished running
wn = 0.25; % initial filter cutoff freq
var0=var; % save an original copy of input variable

if nargin>=2 & 1 % chnge this to (& 0) to manually enforce debug
    debg=varargin{1};
else
    debg=0; % if debg==1, plot will update each time thru while loop & show adaptive filtering and phase identification
end

while ~complete

    foundPhases = 1; % initally assume phases will be found

    % initialize arrays that contains indices that have zero crossings
    pos_crossing = []; % zero crossing with positive slope
    neg_crossing = []; % zero crossing with negative slope
    
    zplus = 0 + 0.001*max(var); % create band around zero to ignore small noise during swing phase
    zminus = 0 - 0.001*max(var);

    % find zero crossings
    i0=round(0.075*length(var));
    for i=i0:length(var)-1 % note: starts at sample i0 to skip any quick positive blips right at heelstrike collision
        % crosses zero with positive slope
        if (var(i)<zminus) && (var(i+1)>zminus)
            pos_crossing = cat(2,pos_crossing,i);
        end
        % crosses zero with negative slope
        if (var(i)>zplus) && (var(i+1)<zplus)
            neg_crossing = cat(2,neg_crossing,i);
        end
    end
    

    % if algorithm finds at least 2 positive & 2 negative slope zero crossings
    if length(pos_crossing)>=2 && length(neg_crossing)>=2 
        indices = [1 pos_crossing(1) neg_crossing(1) pos_crossing(2) neg_crossing(2)];

        % Perform a few simple checks to determine if phases were identified properly
        % (1) check proper ordering of positive and negative slope crossings
        if indices(2)>indices(3) | indices(3)>indices(4) | indices(4)>indices(5) 
            foundPhases=0;
        end

        % (2) for stride, end of pushoff (indice 5) should be arond 60% of gait cycle
        if indices(5)>(0.75*length(var)) | indices(5)<(0.55*length(var))
            foundPhases=0;
        end
        
        % (3) only the end of pushoff should occur after 60% of gait cycle
        if indices(4) > 0.6*length(var) 
            foundPhases=0;
        end
        
        % (4) no phases should be too short in duration
        if min(diff(indices)) < 0.03*length(var)
            foundPhases=0;
        end
        
    else
        foundPhases = 0;
    end


    % if phases are found, then its done. otherwise, try filtering the work rate plot and retrying to find phases of gait
    if foundPhases
        complete=1;
    elseif wn<=0.02 % if cut-off freq gets too low, concede that this method will not find phases
        indices = [];
        fprintf('\nWarning: Phases Not Found\n'); % Output warning to command line
        complete=1;
    else % if it does not find phases initially, try filtering; each time thru loop reduce cut-off freq
        indices = [];
        wn = wn - 0.01;
        [b,a] = butter(3,wn);
        var = filtfilt(b,a,var0);
    end
    
    
    
    %% plot for debugging
    if debg==1
        if foundPhases
            figure(727374); clf; plot(var,'LineWidth',3); ax=axis; hold on; plot([ax(1) ax(2)],[0 0],'k','LineWidth',1)
            hold on; plot([indices(2) indices(2)],[ax(3) ax(4)],'k');
            hold on; plot([indices(3) indices(3)],[ax(3) ax(4)],'k');
            hold on; plot([indices(4) indices(4)],[ax(3) ax(4)],'k');
            hold on; plot([indices(5) indices(5)],[ax(3) ax(4)],'k');
        else
            figure(727374); clf; plot(var,'r','LineWidth',3);
            hold on; plot(var0,'b','LineWidth',3);
            ax=axis; hold on; plot([ax(1) ax(2)],[0 0],'k','LineWidth',1);
        end
    end
    
    
end
    

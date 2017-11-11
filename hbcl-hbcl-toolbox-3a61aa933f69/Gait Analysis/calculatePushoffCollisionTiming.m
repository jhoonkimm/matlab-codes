function [percentpPO1beforeNpercentCO2 pPO1 nCO2] = calculatePushoffCollisionTiming(comWorkRate1, comWorkRate2, varargin)

% CALCULATEPUSHOFFCOLLISIONTIMING based on percentage of pushoff work performed before "half-collision"
%  "Half-collision" is defiend as point in gait cycle when 50% of collision work has been performed.
%  Function returns 3 outputs: (1) percent of pushoff before 50% collision
%  (2) total positive push-off work**, (3) total negative collision work**
%
%  ** Function needs to be updated b/c currently integrates over indices rather than time, so work outputs are not in real units (as of 6/22/10)
%
%  Variable (optional) input arguments:
%    'N': instead of finding pushoff work before 50% (half) collision, find before N% of collision (nominally, N=50)
%
%  Examples:
%   


N = 50; % define percentage of collision

debg=0;

%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'N'
      N = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


%% 
comPhases1 = findPhasesOfGait(comWorkRate1);

startComPushoff = comPhases1(4);
endComCollision = comPhases1(5);

% find end of collision if it ends after push-off ends
comWorkRate2zeros = find(comWorkRate2(1:length(comWorkRate2)-1)<0 & comWorkRate2(2:length(comWorkRate2))>0); % find comWorkRate2 zero crossings
comWorkRate2zeros = comWorkRate2zeros(find(comWorkRate2zeros>comPhases1(5) & comWorkRate2zeros<comPhases1(5)+round(0.1*length(comWorkRate2))));
if ~isempty(comWorkRate2zeros)
    endComCollision = comWorkRate2zeros(1);
end

% find start of collision if it begins before push-off
comWorkRate2zeros = find(comWorkRate2(1:length(comWorkRate2)-1)>=0 & comWorkRate2(2:length(comWorkRate2))<0); % find comWorkRate2 zero crossings
comWorkRate2zeros = comWorkRate2zeros(find(comWorkRate2zeros<comPhases1(4) & comWorkRate2zeros>comPhases1(4)-round(0.1*length(comWorkRate2))));
if ~isempty(comWorkRate2zeros)
    startComPushoff = comWorkRate2zeros(end);
end


po2co = [startComPushoff:endComCollision]'; % array indices from starts of pushoff to end of collision

if debg
    figure(7890); clf; plot([comWorkRate1 comWorkRate2]); ax = axis;
    hold on; plot([po2co(1) po2co(1)],[ax(3) ax(4)],'k');
    hold on; plot([po2co(end) po2co(end)],[ax(3) ax(4)],'k');
end

[pPO1 nPO1] = findPosNegWork(comWorkRate1(po2co));
[pCO2 nCO2] = findPosNegWork(comWorkRate2(po2co));

%% Find percentage of push-off completed before N % of collsion
% initialize
iN = po2co(1);
nCO2_fraction = 0;

while (nCO2_fraction < N/100)
    iN=iN+1;
    [partialpCO2 partialnCO2] = findPosNegWork(comWorkRate2(po2co(1):iN));
    nCO2_fraction = partialnCO2/nCO2;
end

if debg, hold on; plot([iN iN],[ax(3) ax(4)],'r'); end

[partialpPO1 partialnPO1] = findPosNegWork(comWorkRate1(po2co(1):iN));
pPO1beforeNpercentCO2 = partialpPO1;
percentpPO1beforeNpercentCO2 = partialpPO1/pPO1;

NpercentOfnCO2 = N/100*nCO2;

return



% For CESR data
if strcmp(footStrikeOrder,'IP')
    calculatePushoffCollisionTiming(IntaILMpowerStride, ProsILMpowerStride), title('IP foot strike order');
else
    calculatePushoffCollisionTiming(ProsILMpowerStride, IntaILMpowerStride), title('PI foot strike order');
end




function [ComWorkRate, v] = ...
  calculateComWorkRateForSteadyStateGaitRunning( ...
  GRF, avgGaitSpeed, timeDuration, varargin)

% calculateComWorkRateForSteadyStateGaitRunning
%  Returns center of mass (COM) work rate and  COM velocity. Input arguments
%  are GRFs, average walking speed
%  and time duration of stride(s).
%
%  GRFs should be in Nx3 matrix with columns of X,Y,Z and should
%  be data for an integer number of strides (e.g., 1 or 10,
%  but NOT 2.25 strides).
%  +Z is defined as UP.
%  +Y is defined as FORWARD.
%  +X is defined as to the RIGHT.
%
%  Variable input arguments:
%    'mass': used a fixed subject mass that you input (nominally, mass is calculated from vertical GRFs)
%    'subtractVelocitySlope': if set to 1, it removes velcoity slope by assuming no net acceleration across stride (nominally, this is set to 0/false)
%


%% Default variables & constants
mass = [];
g = 9.8;
time = linspace(0,timeDuration,length(GRF));
vgait = avgGaitSpeed;
subtractVelocitySlope = 0;

%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'mass'
      mass = val;
    case 'subtractVelocitySlope'
      subtractVelocitySlope = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end

%% Calculate mass from average vertical GRF over 1 stride
if isempty(mass)
  mass = mean([GRF(:,3)])/g;
end

%% COM acceleration
a = (GRF)/mass - g*[zeros(length(GRF),2) ones(length(GRF),1)];

%% Remove velocity slope. Assumes no net acceleration from start to end of stride
if subtractVelocitySlope
  a =  a - ones(length(a),1)*mean(a);
end

%% COM velocity = integral of accelerations plus constant of integration
v = cumtrapz(time,a);

%% Remove velocity slope. Assumes velocity and start and end of stride is equal - OLD WAY
if subtractVelocitySlope
  vSlope = (v(end,:)-v(1,:))/(length(v)-1);
  vNew = v - [vSlope(1)*[0:1:length(v)-1]' vSlope(2)*[0:1:length(v)-1]' vSlope(3)*[0:1:length(v)-1]'];
  v = vNew;
end

v = v - ones(length(v),1)*mean(v) + ones(length(v),1)*[0 vgait 0]; % assume average velocty across stride(s)

% %% Remove velocity mean/bias
% if subtractVelocityMean
%     v = v - [mean(v(:,1))*ones(length(v),1) mean(v(:,2))*ones(length(v),1)-vgait mean(v(:,3))*ones(length(v),1)];
% end


%% Change in kinetic energy
deltaKE = 0.5*mass*(v(end)^2 - v(1)^2);

%% COM work rate
ComWorkRate = dot(v, GRF, 2);
% LComWorkRate = dot(v,LGRF,2);





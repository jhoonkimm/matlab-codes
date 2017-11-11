function [pos_work neg_work] = findPosNegWork(power, varargin)

% FINDPOSNEGWORK integrate to find positive and negative areas under a curve
%  Calculates positive and negative area under a curve.
% 
%  Same as FINDPOSNEGAREAUNDERCURVE
% 
%  Examples:
%    time = linspace(0,2*pi,100); y = sin(time); figure; plot(time,y); 
%    [pos neg] = findPosNegWork(y,time); % integrates with respect to time
%    [pos neg] = findPosNegWork(y); % integrates over 1:1:length(y)

%  Karl Zelik
%  updated 10/29/09

if nargin==1
    time = 1:length(power);
elseif nargin==2
    time = varargin{1};
else
    fprintf('Error: Incorrect number of input arguments');
    return;
end

if (length(time) == 1)
  % interpret the time input as the stride period
  time = linspace(0, time, length(power));
end

pos_rect = [power.*(power > 0)];
pos_work = trapz(time, pos_rect); 
neg_rect = [power.*(power < 0)];
neg_work = trapz(time, neg_rect);

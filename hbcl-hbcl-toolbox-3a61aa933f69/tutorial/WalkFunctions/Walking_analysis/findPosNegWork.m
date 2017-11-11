function [pos_work neg_work] = findPosNegWork(varin, varargin)

% FINDPOSNEGWORK integrate to find positive and negative areas under a curve
%  Calculates positive and negative area under a curve.
% 
%  Same as FINDPOSNEGAREAUNDERCURVE
% 
%  Examples:
%    time = linspace(0,2*pi,100); y = sin(t); figure; plot(t,y); 
%    [pos neg] = findPosNegWork(y,t); % integrates with respect to t
%    [pos neg] = findPosNegWork(y); % integrates over 1:1:length(y)

%  Paul
%  updated 08/09/10

if nargin==1
    time=1:length(varin);
elseif nargin==2
    time=linspace(0,varargin{1},length(varin));
else
    fprintf('Error: Incorrect number of input arguments');
    return;
end

pos_rect = [varin.*(varin > 0)];
pos_work = trapz(time, pos_rect); 
neg_rect = [varin.*(varin < 0)];
neg_work = trapz(time, neg_rect);

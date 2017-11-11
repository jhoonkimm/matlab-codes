function [F1raw, F2raw, M1raw, M2raw, varargout] = loadANCdata(datafile,varargin)

% LOADANCDATA read in analog data from .anc file and convert to forces, moments
%  Load in analog data. Use force plate calibration matrix to convert to forces
%  (N) and moments (Nm). If LOADANCDATA.m is not located in the working
%  directory, then DATAFILE must include full path.
%
%  VARARGIN Options:
%   'samplefreq': analog sampling frequency of force plate (default 1200 Hz)
%   'adcResolution': ADC bit resolution (default 12 bit)
%   'adcVoltageRange': ADC input voltage range (default 10V)
%   'applyForcePlaCalMatrix': binary option to convert from analog voltage to 
%               forces in N and moments in Nm (default == 1, apply cal matrix) 
%   'filterForces': binary option to filter forces (default == 1, filter)
%   'filterFreq': low-pass cutoff frequency for filtering forces (default 25 Hz)
%   'forcePlaCalMatrix1': 6x6 force plate calibration matrix for force plate 1
%               (default uses Bertec default cal matrix -- diagnonal matrix)
%   'forcePlaCalMatrix1': 6x6 force plate calibration matrix for force plate 1
%               (default uses Bertec default cal matrix -- diagnonal matrix)
%   'forcePlaCalMatrix2': 6x6 force plate calibration matrix for force plate 2
%               (default uses Bertec default cal matrix -- diagnonal matrix)
%
%  IMPORTANT NOTE: BE CAREFUL INPUTTING YOUR OWN FORCE PLATE CAL MATRIX
%   In this code, the calibration matrix is applied in the following way:
%   F1M1raw = [F1volts M1volts]*forcePlaCalMatrix1;
% 
%   This may vary from other software packages. Check force plate technical 
%   documentation to determine correct orientation your calibration matrix,
%   as it may or may not need to be transposed.
% 
%  Examples:
%   datafile = '/Users/kzelik/ ... /hx_shod_walking_on_both_belts_1.anc';
%   [F1raw, F2raw, M1raw, M2raw, varargout] = loadANCdata(datafile); % load in using all default settings
%
%   [F1raw, F2raw, M1raw, M2raw, varargout] = loadANCdata(datafile,'adcResoltuion;,14); % update ADC resolution to match DAQ system used

%  Karl Zelik
%  updated 10/29/09


%% Default variables
samplefreq = 1200; % HNL Motion Analysis motion capture system sampling rate (note: Vicon system may be different)
applyForcePlaCalMatrix = 1;
filterForces = 1; 
filterFreq = 25;
forcePlaCalMatrix1 = diag([500 500 1000 600 300 300]); forcePlaCalMatrix2 = forcePlaCalMatrix1; % Bertec default calibration matrix
adcResolution = 12; % 12 bit ADC
adcVoltageRange = 10; % 10V input range


%% Optional input arguments
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'samplefreq'
        samplefreq = val;
    case 'adcResolution'
        adcResolution = val;
    case 'adcVoltageRange'
        adcVoltageRange = val;
    case 'applyForcePlaCalMatrix'
        applyForcePlaCalMatrix = val;
    case 'filterForces'
        filterForces = val;
    case 'filterFreq'
        filterFreq = val;
    case 'forcePlaCalMatrix1'
        forcePlaCalMatrix1 = val;
    case 'forcePlaCalMatrix2'
        forcePlaCalMatrix2 = val;
    otherwise
      error('\nError: incorrect varargin\n')
  end
end


%% Load data
alldata = textread(datafile,'%s','delimiter','\n','whitespace','');
valuesstart = min(findcellstr(alldata,'Range')) + 1;
data = dlmread(datafile,'\t',valuesstart,0);
F1analog = data(:,2:4); F2analog = data(:,8:10);
M1analog = data(:,5:7); M2analog = data(:,11:13);
% Extra analog data
if size(data,2) > 13
    varargout{1} = data(:,14:size(data,2)); % any other analog data
end


%% Convert analog bits to voltages
% HNL force treadmill: convert bits to volts (10V/4095bits) ; -5V-5V gives 10V range; 12bit ADC: 2^12 = 4096 bits (4095 bit steps)
bits2volts = (2^adcResolution-1)/adcVoltageRange;
F1volts = F1analog./bits2volts; 
F2volts = F2analog./bits2volts;
M1volts = M1analog./bits2volts;
M2volts = M2analog./bits2volts;


%% Convert voltages to forces using forcePlaCal matrices
if applyForcePlaCalMatrix 
    F1M1raw = [F1volts M1volts]*forcePlaCalMatrix1; F1raw = F1M1raw(:,1:3); M1raw = F1M1raw(:,4:6);
    F2M2raw = [F2volts M2volts]*forcePlaCalMatrix2; F2raw = F2M2raw(:,1:3); M2raw = F2M2raw(:,4:6);  
else
    F1raw = F1analog; F2raw = F2analog; % just use analog voltages
end


%% Filter forces
if filterForces
    % Filter forces
    [w,q] = butter(3,filterFreq/(samplefreq/2)); % wn of 1 --> 1/2 sample freq.  Forces sampled at 1200 Hz 
    F1raw_unfilt = F1raw; F2raw_unfilt = F2raw;
    F1raw = filtfilt(w,q,F1raw); F2raw = filtfilt(w,q,F2raw);
end


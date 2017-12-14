function [metOutput] = oxyconProcessing(metFilename,varargin)
%  
% [metOutput] = oxyconProcessing(metFilename)
%  Allow for processing of Oxycon (Kuo) metabolic data into metabolic rates (metOutput). Outputs in kJ s-1
%  metFilename must be .txt file.
%
%  Flags:
%
%   brockway
%      Default is to use RER and its caloric equivalent (Péronnet, 1991),
%      with the default value of 0. Brockway will use equation: E = 16.58kJ/L *VO2
%      + 4.51 kJ/L *VCO2. (Brockway, 1987)
%   smooth
%      Default is 1. 3-point moving average smoothing filter.
%   massSpecific
%      Default is 1. Mass-specific metabolic rates given for RER caloric
%      eq.
%
%  RER file 'EnergyConversionTable_1991.txt' is needed for caloric
%  equivalent calculations
%
%  Created by Jay Kim (11/09/2017)

%% Default input arguments
brockway = 0;
filt = 1;
massSpecific = 1;

%% Optional input arguments

for i = 1:2:length(varargin)
  switch varargin{i}
    case 'brockway'
    brockway = varargin{i + 1};
    case 'smooth'
    filt = varargin{i + 1};
    case 'massSpecific'
    massSpecific = varargin{i + 1};
  end
end

%% Load textfile

f = fopen(metFilename);
strData = textscan(f,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t');
fclose(f);

time = vertcat(strData{1,1}{4:end});
VO2kg = str2double([{strData{1,8}{4:end}}']);
VO2 = str2double([{strData{1,7}{4:end}}']);
VCO2 = str2double([{strData{1,9}{4:end}}']);
RER = str2double([{strData{1,10}{4:end}}']);

if filt
  VO2kg = smooth(VO2kg,3,'moving');
  VO2 = smooth(VO2,3,'moving');
  VCO2 = smooth(VCO2,3,'moving');
end



%% Processing textfile string
t = [];
for i = 1:length(time) %converting string time to value time
    tempTime = str2double(regexp(time(i,:),':','split'));
    t = [t;tempTime*[60;1]];
end

%% Metabolic analysis

%RER Caloric equivalent 
if brockway == 0
%retrieving energy conversion table
  f = fopen(fullfile(cd,'EnergyConversionTable_1991.txt'));
  textscan(f, '%s %s', 1,'Delimiter', '\t');
  enConv = textscan(f, '%f %f', 61,'Delimiter', '\t');
  enConv = [enConv{1} enConv{2}];
  fclose(f);

  cal = [];
   for i = 1:length(RER)
          r = (round(RER(i)*100))/100;
          if r < 0.7050 %lower caloric end
              r = 0.7036;
          elseif r > 0.995 %higher caloric end
              r = 0.996;
          end
      [rw,~] = find(enConv==r);
      if massSpecific == 1;
        cal = [cal; VO2kg(i)*enConv(rw,2)*(1/60)]; %gives W kg-1
      elseif massSpecific == 0;
        cal = [cal; VO2(i)*enConv(rw,2)*(1/60)]; %gives W kg-1
      end
   end
   
elseif brockway == 1
     
   cal = ((VCO2*4.51)+(VO2*16.58))/60; %gives W
   
end

  metOutput = cal;



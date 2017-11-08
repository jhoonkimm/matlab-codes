function [metOutput] = oxyconProcessing(metFilename,varargin)
% ------------------------------------------------------ %
%  [metOutput] = oxyconProcessing(metFilename)
%
%  Created by Jay Kim (Nov. 8th 2017)
%  Allow for processing of Oxycon (Kuo) metabolic data into metabolic rates and
%  cost of transport (metOutput). Outputs in kJ s-1
%
%  Need RERFile to proceed. metFilename must be .txt file.
% ------------------------------------------------------ %

%% Default input arguments
COT = 0;

%% Optional input arguments

for i = 1:2:length(varargin)
  if (strcmp('COT', varargin{i}))
    argCOT = varargin{i + 1};
  end
end

%% Load textfile

f = fopen(metFilename);
strData = textscan(f,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t');
fclose(f);

time = vertcat(strData{1,1}{4:end});
VO2kg = smooth(str2double([{strData{1,8}{4:end}}']),3,'moving');
RER = str2double([{strData{1,10}{4:end}}']);

%% Processing textfile string
t = [];
for i = 1:length(time) %converting string time to value time
    tempTime = str2double(regexp(time(i,:),':','split'));
    t = [t;tempTime*[60;1]];
end


%% Metabolic analysis

%retrieving energy conversion table
f = fopen(fullfile(cd,'EnergyConversionTable_1991.txt'));
textscan(f, '%s %s', 1,'Delimiter', '\t');
enConv = textscan(f, '%f %f', 61,'Delimiter', '\t');
enConv = [enConv{1} enConv{2}];
fclose(f);

cal = [];
 for i = 1:length(RER)
        r = (round(RER(i)*2,2))/2;
        if r < 0.7050 %lower caloric end
            r = 0.7036;
        elseif r > 0.995 %higher caloric end
            r = 0.996;
        end
    [rw,~] = find(enConv==r);
    cal = [cal; VO2kg(i)*enConv(rw,2)*(1/60)]; %gives kJ/s
 end
    
metOutput = cal;



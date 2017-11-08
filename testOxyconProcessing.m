clear all; close all; clc;

metFilename = {'C:\Users\BertramLab\Documents\MATLAB\Jay\SBD_001\10146402_resting.txt'
               'C:\Users\BertramLab\Documents\MATLAB\Jay\SBD_001\10146402_pws.txt'
               'C:\Users\BertramLab\Documents\MATLAB\Jay\SBD_001\10146402_150baseline.txt'
               'C:\Users\BertramLab\Documents\MATLAB\Jay\SBD_001\10146402_75baseline.txt'
               'C:\Users\BertramLab\Documents\MATLAB\Jay\SBD_001\10146402_split.txt'};
tot = [];
for i = 1:length(metFilename)
    metOutput = oxyconProcessing(metFilename{i},'COT',0);
    tot = [tot;metOutput];
end



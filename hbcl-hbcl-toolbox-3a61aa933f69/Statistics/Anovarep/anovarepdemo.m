%ANOVAREPDEMO
%This is a demo on the use of ANOVAREP function
%Example
%Run anovarepdemo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Anovarep: compute the Anova for repeated measures and
% Holm-Sidak test for multiple comparisons if Anova is positive. 
% http://www.mathworks.com/matlabcentral/fileexchange/18746

clear all
load anrepdemo
Table{13,5}=[];   
Table(1,1:5)={'Soldier' 'At begin of training' 'At seize' '12 hours after seize' '48 hours after seize'};
Table(2:13,1)={'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
Table(2:13,2:5)=num2cell(x);
H=figure('Position',[360 80 644 213]);
statdisptable(Table, 'Data', 'Testosterone levels (ng/dL)', 'Biol Psychiatry, 2000; 47:89-1901 - press a key to continue',[-1 -1 0 -1 2 4],H);
pause
set(H,'Position',[600 50 644 213])
anovarep(x)

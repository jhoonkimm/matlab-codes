function [typestr,folderstr,namestr,matlabvarstr]=buildV3Dexportstrs(types,folders,names,matlabvars)
% builds export strings from cell array columns of signal names

typestr=strBuilder(types,'/SIGNAL_TYPES=');
folderstr=strBuilder(folders,'/SIGNAL_FOLDER=');
namestr=strBuilder(names,'/SIGNAL_NAMES=');
matlabvarstr=strBuilder(matlabvars,'/OUTPUT_NAMES=');

function str=strBuilder(sigs, command)
sigs=strcat(sigs,'+');
prestr=strcat(sigs{:});
str=[command prestr(1:end-1) '\r\n'];

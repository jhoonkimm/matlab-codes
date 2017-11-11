function [out1,out2,out3] = blah(in);
% [genscale, scales, offsets] = petergetscales(c3dobj)
% accepts an input c3dserver object (formed by command "test = c3dserver")
% returns: 
%   general Analog scaling factor, 
%   individual channel scaling factors, 
%   individual channel offsets (in # of levels)


nIndex = in.GetParameterIndex('Analog','Gen_Scale');
out1 = in.GetParameterValue(nIndex,0);

nIndex = in.GetParameterIndex('Analog','Scale');
nLength = in.GetParameterLength(nIndex);
for I = 0:nLength-1
    out2{I+1} = in.GetParameterValue(nIndex,I);
end

nIndex = in.GetParameterIndex('Analog','Offset');
nLength = in.GetParameterLength(nIndex);
for I = 0:nLength-1
    out3{I+1} = in.GetParameterValue(nIndex,0);
end


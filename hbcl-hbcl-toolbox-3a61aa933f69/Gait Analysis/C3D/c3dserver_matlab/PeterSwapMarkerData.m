clear all

[filenames, pathnames] = uigetfile('*.c3d','Pick a C3D file to analyze','MultiSelect','on');
if ~iscell(filenames); filenames = {filenames}; end

namestoswap = input('Enter Cell Array, 2 Point Names to Switch: ');

itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package

for ii = 1:length(filenames)
    
openc3d(itf, 0, [pathnames filenames{ii}])    % applies the correct open options to the file you choose    

framestoswap = input([filenames{ii} ': Enter Frames in which to Swap Data: ']);
if isempty(framestoswap)
    framestoswap = [itf.GetVideoFrame(0):itf.GetVideoFrame(1)];
end

PointLabelsIndex = itf.GetParameterIndex('Point','Labels');  % typ. 5

for jj = 0:itf.GetParameterLength(PointLabelsIndex)
    labels{jj+1} = itf.GetParameterValue(PointLabelsIndex,jj);
end

swap{1} = find(strcmp(namestoswap{1},labels))-1;
swap{2} = find(strcmp(namestoswap{2},labels))-1;

for jj = framestoswap
    for comp = 0:2
        tmp = itf.GetPointData(swap{1},comp,jj,'0');
        itf.SetPointData(swap{1},comp,jj,itf.GetPointData(swap{2},comp,jj,'0'));
        itf.SetPointData(swap{2},comp,jj,tmp);
    end
end


% itf.SaveFile([pathnames filenames{ii}], 1)
% closec3d(itf)
end



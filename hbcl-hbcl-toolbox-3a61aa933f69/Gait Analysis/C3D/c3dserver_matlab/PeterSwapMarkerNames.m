clear all

[filenames, pathnames] = uigetfile('*.c3d','Pick a C3D file to analyze','MultiSelect','on');
if ~iscell(filenames); filenames = {filenames}; end

namestoswap = input('Enter Cell Array, 2 Point Names to Switch ');

itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package

for ii = 1:length(filenames)
    
openc3d(itf, 0, [pathnames filenames{ii}])    % applies the correct open options to the file you choose    

PointLabelsIndex = itf.GetParameterIndex('Point','Labels');  % typ. 5

for jj = 0:itf.GetParameterLength(PointLabelsIndex)
    labels{jj+1} = itf.GetParameterValue(PointLabelsIndex,jj);
end

swap{1} = find(strcmp(namestoswap{1},labels))-1;
swap{2} = find(strcmp(namestoswap{2},labels))-1;

itf.SetParameterValue(PointLabelsIndex,swap{1},'blah')
itf.SetParameterValue(PointLabelsIndex,swap{2},namestoswap{1})
itf.SetParameterValue(PointLabelsIndex,swap{1},namestoswap{2})

itf.SaveFile([pathnames filenames{ii}], 1)
closec3d(itf)
end


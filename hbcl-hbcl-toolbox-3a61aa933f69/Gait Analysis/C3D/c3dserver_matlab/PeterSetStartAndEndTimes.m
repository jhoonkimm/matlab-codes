clear all

[filenames, pathnames] = uigetfile('*.c3d','Pick a C3D file to Set Start and End Times','MultiSelect','off');
if ~iscell(filenames); filenames = {filenames}; end

itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package

for ii = 1:length(filenames)
    
openc3d(itf, 0, [pathnames filenames{ii}])    % applies the correct open options to the file you choose    
lastframe = itf.GetVideoFrame(1);
videorate = itf.GetVideoFrameRate();

startgood = 1; startbad = lastframe;
startgood = input('Enter First Good Frame [1] ');
startbad = input('Enter First BAD Frame [end+1] ');


eventnum(1) = itf.AddEvent('GOOD','0',startgood)
eventnum(2) = itf.AddEvent('GOOD','1',startbad)

if any(eventnum==-1)
    fprintf('Events Not All Created')
end

% itf.SaveFile([pathnames filenames{ii}], 1)
% closec3d(itf)
end



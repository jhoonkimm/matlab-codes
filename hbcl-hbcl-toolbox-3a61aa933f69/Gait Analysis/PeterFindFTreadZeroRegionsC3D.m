% PETERFINDFTREADZEROREGIONSC3D
%
% This script reads C3D files (selected with a gui) and, from FZ, finds a
% time when each foot is off the ground. It assigns indices for the
% relevant frames to the new parameter "Zeros" in a copy of the C3D file
% (marked "_modified"). This operation should be done for ALL force
% treadmill data before processing in Visual3D, because otherwise Visual3D
% assumes that the first 10 frames of data (on Force channels) are
% Zero, and uses these to offset the rest of the data. For force treadmill
% data, this Visual3D default usually gives good results for one foot but
% not for the other, since it was on the treadmill during the assumed
% "zero" time.  
%
% The script works based on FZ data. 
% 1) It reads the channel, then smoooths the hell out of it (a simple
% moving average filter run 700 times over the data set), so that each
% step's force profile is reduced to a simple sinusoid-like form. 
% 2) It picks out the local maxima (These are actually local minima in
% force, but the analog signal comes in Negative, so we look for maxima to
% find smallest magnitude). These maxima are typically in the middle of the
% time when the foot is actually off the force plate. 
% 3) It computes the mean of the real signal within the window specified by
% "windowinds" around each local max, and weeds out any results that are
% not within 1 SD of the MEDIAN, and any that are too close to the
% beginning or end. This is to get rid of any "freaky" steps.
% 4) It selects the one zero region whose mean is closest to the MEAN of the
% remaining means.  
% 5) It reports a video-frame region around this point as the "zero
% region" in the new parameter "Zeros". The window goes from windowinds(1)
% to windowinds(2) frames relative to the center point chosen. This is
% determined by the interaction of the raw force with the smoothing
% algorithm, and may need to be tuned for a given data set. 
% 6) It graphs the original Force, and highlights the Zero Region in Red.
% Visual inspection is important.
%
% Parameters to be tuned are: 
% fthresh: A force threshold that the local maxima need to exceed to not be
% excluded. This may need to be tuned, especially if step frequency gets
% high.
% windowinds: Relative Indices for the window to report, w.r.t. the chosen
% local maximum. Can be anything, but visual inspection should confirm that
% regions thus chosen really are zero. This is originally specified to
% yield 20 frames (achieved by windowinds(2)-windowinds(1) = 19, with one
% for the nominal frame.) 
% smoothreps: Repetitions of the smoothing filter. 
%
% Note that this function could be generalized to use the analog-to-video
% ratio, instead of the assumed "10", but that has not been done. 

% Peter Adamczyk 2009-11-23

clear ; close all; 

fthresh = 250;  % force threshold to determine which among the "local maxima" of FZ are to be included in the "good" set.
windowinds = [-5 14];
smoothreps = 700;

analogwindowinds = 10*windowinds + [-9 10];

filenames = {};
pathnames = {};
for i = 1:100
    fnms = {}; pnms = {};
    [fnms, pnms{1}] = uigetfile('*.c3d','Pick C3D files to analyze','MultiSelect','on');
    if ~iscell(fnms); if ~fnms; break; end; fnms = {fnms}; end
    pnms = repmat(pnms,length(fnms),1);
    filenames = cat(1,filenames,fnms{:});
    pathnames = cat(1,pathnames,pnms{:});
end
if ~filenames{end}; error('File Selection Cancelled'); end

for qq = 1:length(filenames)  %% do this for many files. 
    
pathname = pathnames{qq};
filename = [pathname filenames{qq}]; 
subjname = filenames{qq}(1:2);

alreadymodified = strcmp(filename(end-11:end),'modified.c3d');

if alreadymodified  % if this has file has gone thru the processing loop of MATLAB-->Visual3D-->(now MATLAB)
    modifiedpathname = pathname;
    modifiedfilename = filename;
else  % if this file has not gone through the loop, treat it like the first time. 
    modifiedpathname = [pathname ' ModifiedForZeros\'];
    modifiedfilename = [modifiedpathname filenames{qq}(1:end-4) '_modified.c3d'];
    goodcopy = mkdir(pathname,' ModifiedForZeros');
    goodcopy = copyfile( filename, modifiedfilename, 'f');
    if ~goodcopy
        error(['File copy failed:: \n    ' filename]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start C3DServer stuff

itf = actxserver('C3DServer.C3D');  % makes "itf" a COM object for the "c3dserver" package
openc3d(itf, 0, [modifiedfilename])    % applies the correct open options to the file you choose


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Force, Moment, COP
%% First get a "zero" region for each treadmill to remove offsets. 

%% IFF there is already a "ZEROS" parameter in the file, skip this part and
%% use the values it already has:
zerosparmind = itf.GetParameterIndex('Force_Platform','ZEROS');
zeroparmind = itf.GetParameterIndex('Force_Platform','ZERO');
if zerosparmind ~= -1  %% there is already a "ZEROS" parameter... (i.e., this file has not been modified yet)
    cont = input('This File already has a Zeros parameter. Overwrite?  ');
else
    cont = 1;
end

if cont
    zeroregion(1) = itf.GetParameterValue(zeroparmind, 0);
    zeroregion(2) = itf.GetParameterValue(zeroparmind, 1);
    itf.SetParameterValue(zeroparmind,0,0); %%% set the MINimum zeroing frame to 0, so that no zero adjustment is made.
    itf.SetParameterValue(zeroparmind,1,0); %%% set the MAXimum zeroing frame to 0, so that no zero adjustment is made.
    % get UNCORRECTED values for fz1 andfz2
    fz1 = double(cell2mat(itf.GetForceData(2,0,-1,-1)));
    fz2 = double(cell2mat(itf.GetForceData(2,1,-1,-1)));

%% find the zone that FP1 should treat as "zero"
    [mag, index] = lmax(smooth(fz1,smoothreps), 3);  %% fz1 comes in with negative sign, so use LMAX
    goods = find((index<size(fz1,1)-(analogwindowinds(2)+10)) & (index>-analogwindowinds(1) + 10) & (abs(mag) < fthresh)); % 32767 is 2^15-1.  This is the upper limit of the parameter FORCE_PLATFORM::ZERO because of the 16-bit signed integer number format
    index = index(goods);
    %     mag = mag(goods);
    mag = arrayfun(@(x) mean(fz1(x+analogwindowinds(1):x+analogwindowinds(2))), index);
    magmn = median(mag);
    magsd = std(mag);
    goods = find( (abs(mag - magmn) < magsd) );
    magmn = mean(mag(goods));
    [mag,idd] = min( abs(mag-magmn) );
    index = floor(index(idd)/10);  % makes this an even multiple of 10, so it can easily be written in terms of Video Frames
    zeroindsfp1{qq} = [index+windowinds(1), index+windowinds(2)];
    index = -9+10*zeroindsfp1{qq}(1):10*zeroindsfp1{qq}(2) ;
    offset1 = mean(fz1( index ));
    figure(1000*qq+1); plot(fz1); hold on; plot(index, fz1(index), 'r','MarkerSize',12); text(mean(xlim),mean(ylim),num2str(median(index))); title({'zGRF for "ZEROS"',filenames{qq}},'Interpreter','none')
    if isempty(zeroindsfp1{qq}); zeroindsfp1{qq} = [0 0]; end    
    display(['FP1 zero indices: ' int2str(zeroindsfp1{qq})]);

%% find the zone that FP2 should treat as "zero"
    [mag, index] = lmax(smooth(fz2,smoothreps), 3);  %% fz2 comes in with negative sign, so use LMAX
    goods = find((index<size(fz2,1)-(analogwindowinds(2)+10)) & (index>-analogwindowinds(1) + 10) & (abs(mag) < fthresh)); % 32767 is 2^15-1.  This is the upper limit of the parameter FORCE_PLATFORM::ZERO because of the 16-bit signed integer number format
    index = index(goods);
    %     mag = mag(goods);
    mag = arrayfun(@(x) mean(fz2(x+analogwindowinds(1):x+analogwindowinds(2))), index);
    magmn = median(mag);
    magsd = std(mag);
    goods = find( (abs(mag - magmn) < magsd) );
    magmn = mean(mag(goods));
    [mag,idd] = min( abs(mag-magmn) );
    index = floor(index(idd)/10);  % makes this an even multiple of 10, so it can easily be written in terms of Video Frames
    zeroindsfp2{qq} = [index+windowinds(1), index+windowinds(2)];
    index = -9+10*zeroindsfp2{qq}(1):10*zeroindsfp2{qq}(2) ;
    offset2 = mean(fz2( index ));
    figure(1000*qq+2); plot(fz2); hold on; plot(index, fz2(index), 'r','MarkerSize',12); text(mean(xlim),mean(ylim),num2str(median(index))); title({'zGRF for "ZEROS"',filenames{qq}},'Interpreter','none')
    if isempty(zeroindsfp2{qq}); zeroindsfp2{qq} = [0 0]; end    
    display(['FP2 zero indices: ' int2str(zeroindsfp2{qq})]);

    %% Now Save the Modified file with the new Offset regions:
        modsave = itf.AddParameter('ZEROS','Array of Baseline Indices for Multiple Force Platforms', 'FORCE_PLATFORM', '0', int16(2), int16(2), int16([2 2]), int16([zeroindsfp1{qq} zeroindsfp2{qq}]));
        modsave = itf.SaveFile('',-1);
        if ~modsave; error('Modified C3D File Not Saved \n'); else; fprintf('Modified C3D File Saved Successfully!\n'); end

else  % Else there _is_ a "ZEROS" parameter, indicating that this file has been processed before.
    %     zeroregion(1) = itf.GetParameterValue(zeroparmind, 0);
    %     zeroregion(2) = itf.GetParameterValue(zeroparmind, 1);
    %     zeroindsfp1{qq}(1) = itf.GetParameterValue(zerosparmind,0);  % the File stores "video frames" for Visual3D, but C3DServer uses Analog Sample Number. The equivalent for Video [5:10] is Analog [41:100].
    %     zeroindsfp1{qq}(2) = itf.GetParameterValue(zerosparmind,1);
    %     zeroindsfp2{qq}(1) = itf.GetParameterValue(zerosparmind,2);
    %     zeroindsfp2{qq}(2) = itf.GetParameterValue(zerosparmind,3);
    continue
end
drawnow; 

% %%%NOW Get the Data
% 
% %FP1
% %% tell the c3d file that this zone (201 points?) is
% %% the correct "zeroregion"
% itf.SetParameterValue(zeroparmind,0,zeroindsfp1{qq}(1));
% itf.SetParameterValue(zeroparmind,1,zeroindsfp1{qq}(2));
% 
% % %% Get the data
% % fx1 = -double(cell2mat(itf.GetForceData(0,0,-1,-1)));
% % fy1 = double(cell2mat(itf.GetForceData(1,0,-1,-1)));
% % fz1 = -double(cell2mat(itf.GetForceData(2,0,-1,-1)));
% % cx1 = double(cell2mat(itf.GetMomentData(0,0,-1,-1)));
% % cy1 = double(cell2mat(itf.GetMomentData(1,0,-1,-1)));
% % mz1 = -double(cell2mat(itf.GetMomentData(2,0,-1,-1)));
% 
% [fx1 fy1 fz1 cx1 cy1 mz1] = getFfromAnalog(itf,1);
% 
% %FP2
% %% tell the c3d file that this zone (201 points?) is
% %% the correct "zeroregion"
% itf.SetParameterValue(zeroparmind,0,zeroindsfp2{qq}(1));
% itf.SetParameterValue(zeroparmind,1,zeroindsfp2{qq}(2));
% 
% % %% Get the data
% % fx2 = -double(cell2mat(itf.GetForceData(0,1,-1,-1)));
% % fy2 = double(cell2mat(itf.GetForceData(1,1,-1,-1)));
% % fz2 = -double(cell2mat(itf.GetForceData(2,1,-1,-1)));
% % cx2 = double(cell2mat(itf.GetMomentData(0,1,-1,-1)));
% % cy2 = double(cell2mat(itf.GetMomentData(1,1,-1,-1)));
% % mz2 = -double(cell2mat(itf.GetMomentData(2,1,-1,-1)));
% 
% [fx2 fy2 fz2 cx2 cy2 mz2] = getFfromAnalog(itf,2);

% % Reset the "zeroregion" to its initial values:
% itf.SetParameterValue(zeroparmind,0,zeroregion(1)); 
% itf.SetParameterValue(zeroparmind,1,zeroregion(2)); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

closec3d(itf)
end

function status = getdatafromv3d_perturbation()
%
%
%


%% prompt for data file / signaltable and open pipeline file for editing
% [f,p]=uigetfile('*.cmo','Select workspace');
% p='C:\Users\xyf\Desktop\Xiao-Yu\Obese Data\';
% f='OBGM06_125_0.cmo';
% p=uigetdir();
% p=[p '\'];
workingDirectory = pwd;
p='D:/xyfu/Local Perturbation/Sarah 26Jul2016/';
cd(p)
fullPath = [pwd filesep];
cd(workingDirectory)
% p='C:\Users\Xiao-Yu\Dropbox\Perturbation Xiao-Yu\Xiao-Yu 10Nov12\';
% p='C:\Users\xyf\Dropbox\Perturbation Xiao-Yu\Xiao-Yu 10Nov12\';
modelFile=[fullPath 'Sarah_ModelBuilding.cmo'];
files=dir([p '*_intWithPert.c3d']);
filenames={files.name};


% [sigf,sigp]=uigetfile([p '*.xlsx'],'Select signals table');
sigp=[p 'pipeline\'];
sigf='Signals.xlsx';

fid1=fopen([p 'pipeline\processingPipeline.v3s']);
processing=fread(fid1);
fclose(fid1);

fid2=fopen([p 'pipeline\Angular_Momentum_XY_lite.v3s']);
angMom=fread(fid2);
fclose(fid2);

[inertiaPipeline,inertiaType,inertiaFolder,inertiaNamestr,inertiaMatlabstr]=generateV3DInertiaPipeline([p 'pipeline\InertiaMass.v3s'],{'RFT','RSK','RTH','RAR','RFA','RPV','RTA','LFT','LSK','LTH','LAR','LFA'});

fid3=fopen(inertiaPipeline);
inertiaPipe=fread(fid3);
fclose(fid3);

[num,txt]=xlsread([sigp,sigf]);
for k=1:8%5:7%length(filenames)
  pipelinefile = [p filenames{k}(1:end-4) '_exp.v3s'];
  pipelinefileFull = [fullPath filenames{k}(1:end-4) '_exp.v3s'];
  fid=fopen(pipelinefile,'w');
  if fid > 0
    
    %% Open model builder cmo
    fprintf(fid,'File_Open\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n',modelFile);
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    %% Open motion file
    fprintf(fid,'File_Open\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n',[fullPath filenames{k}]);
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');
    
    %% Assign calibration file
    fprintf(fid,'Assign_Model_File\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',modelFile);
    fprintf(fid,'/MOTION_FILE_NAMES=%s\r\n', [fullPath filenames{k}]);
    fprintf(fid,'! /REMOVE_EXISTING_ASSIGNMENTS=FALSE\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    %% write processing
    fwrite(fid,processing);
    
    %% write angular momentum
    fwrite(fid,angMom);
    
     %% write inertia
    fwrite(fid,inertiaPipe);
    
    %% Save CMO
    fprintf(fid,'File_Save_As\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n', [fullPath filenames{k}(1:end-4) '_workspace.cmo']);
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');
    
    %% write export
    fprintf(fid,'Export_Data_To_Matfile\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n',[fullPath filenames{k}(1:end-4) '_raw_data' '.mat']);

    [typestr,folderstr,namestr,matlabvarstr]=buildV3Dexportstrs([txt(2:end,1); inertiaType],...
                                                                [txt(2:end,2); inertiaFolder],...
                                                                [txt(2:end,3); inertiaNamestr],...
                                                                [txt(2:end,4); inertiaMatlabstr]);
    fprintf(fid,typestr);
    fprintf(fid,namestr);
    fprintf(fid,folderstr);
    fprintf(fid,matlabvarstr);
    fprintf(fid,'! /PARAMETER_NAMES=\r\n');
    fprintf(fid,'! /PARAMETER_GROUPS=\r\n');    
    fprintf(fid,'! /OUTPUT_PARAMETER_NAMES=\r\n');    
    fprintf(fid,'/USE_NAN_FOR_DATANOTFOUND=TRUE\r\n');    
    fprintf(fid,';\r\n');  
    fprintf(fid,'\r\n');    

    %% write footer
    %exit
    fprintf(fid,'Exit_Workspace\r\n');
    fprintf(fid,';\r\n');
    %% close v3s and execute

    status = fclose(fid);
  else
    status = -1;
    disp(['Error opening pipeline file ' pipelinefile]);
  end
    dos(['"C:\Program Files (x86)\Visual3D v6\Visual3D.exe" /s "' pipelinefileFull '"']) % '&'])
%     dos(['"C:\Program Files (x86)\Visual3D v4_96\Visual3D.exe" /s "' pipelinefileFull '"']) % '&'])
%     dos(['"C:\Program Files (x86)\Visual3D v4\Visual3D.exe" /s "' pipelinefile '"']) % '&'])
end

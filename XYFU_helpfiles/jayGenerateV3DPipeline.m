function status = jayGenerateV3DPipeline(path,folder,varargin)
%   jayGenerateV3DPipeline(strDir)  
%   This function demonstrates how to generate and
%   save a Visual3D pipeline script to a user specified directory.  This
%   particular example will:
%   1.  Create a new workspace (File|New).
%   2.  Create a hybrid model (Model|Create|Hybrid Model from C3DFile).
%   3.  Apply a model template (Model|Apply Model Template).
%   4.  Insert a motion file (File|Open).
%   5.  Assign the model to the motion file (Model|Assign Model to Motion
%       files).
%   6.  Calculate the position of the left shank with respect to the left
%       thigh (Model|Compute Model Based Data).
%   7.  Export a LINK_MODEL_BASED signal as a MATLAB .mat file.
%   8.  Save the workspace (File|Save Workspace).
%
%
%   Copyright (C) 2008 C-Motion, Inc.  All rights reserved.
%   For support contact support@c-motion.com.  
%%  Flags:
%
%   trial
%       Specify which trial in the folder you are looking to make pipeline for.
%       i.e. 'pws','fast','slow','split','post'
%   allFiles
%       Make pipeline that loads all files in folder. Default is 0
%
%  Created by Jay Kim (12/17/2017)

%% Default input arguments
trial = [];
allFiles = 0;

%% Optional input arguments
for i = 1:2:length(varargin)
  switch varargin{i}
    case 'trial'
    trial = varargin{i + 1};
    case 'allFiles'
    allFiles = varargin{i + 1};
  end
end
%% Loading filenames for the function

dirSBD = dir(fullfile(path,folder));
dirSBDnames = {dirSBD.name};

idx = indexProcessing(path,folder,'integrated',1);

strFNPipeline = sprintf('%s_%s.v3s',folder,upper(trial));
strFNModelC3D = dirSBDnames{idx.static};
strFNTemplate = sprintf('%s.mdh',folder);
%% If making pipeline for single trial
if allFiles==0
    
strFNPipeline = sprintf('%s_%s.v3s',folder,upper(trial));
% Error messages
if isempty(trial)
    error('No trial has been selected. Please add a trial to the function')
end

if ~isstr(trial)
    error('Specified trial value must be string value!')
end
trialCell = {'pws','fast','slow','split','post'};
if mean(strcmp(trial,trialCell))==0
    error('Incorrect trial string! Try with pws, fast, slow, split or post.')
end
switch trial
    case 'pws'
    strFNMotionC3D = dirSBDnames{idx.c3d.pws};
    case 'fast'
    strFNMotionC3D = dirSBDnames{idx.c3d.fast};
    case 'slow'
    strFNMotionC3D = dirSBDnames{idx.c3d.slow};
    case 'split'
    strFNMotionC3D = dirSBDnames{idx.c3d.split};
    case 'post'
    strFNMotionC3D = dirSBDnames{idx.c3d.post};
end
strFNMAT = sprintf('%s_%s.mat',trial,folder);
strFNWorkspace = sprintf('%s_%s.cmo',trial,folder);

fid = fopen(fullfile(path,folder,strFNPipeline),'w');

if fid > 0
    % Create a new workspace (File|New).
    fprintf(fid,'File_New\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Create a hybrid model (Model|Create|Hybrid Model from C3DFile).
    fprintf(fid,'Create_Hybrid_Model\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'! /RANGE=ALL_FRAMES\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Apply a model template (Model|Apply Model Template).
    fprintf(fid,'Apply_Model_Template\r\n');
    fprintf(fid,'/MODEL_TEMPLATE=%s\r\n',fullfile(path,folder,strFNTemplate));
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Build model
    fprintf(fid,'Build_Model\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'! /REBUILD_ALL_MODELS=FALSE\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Insert a motion file (File|Open).
    fprintf(fid,'Open_File\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n',fullfile(path,folder,strFNMotionC3D));
    fprintf(fid,'! /CALIBRATION_FILE=\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Assign the model to the motion file (Model|Assign Model to Motion
    % files).
    fprintf(fid,'Assign_Model_File\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'/MOTION_FILE_NAMES=%s\r\n',fullfile(path,folder,strFNMotionC3D));
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    computeId = fopen('E:\RESEARCH\MDSC 508\SimplifiedProcessingPipeline_2.txt');
    while ~feof(computeId)
        tline = fgetl(computeId);
        fprintf(fid,'%s\r\n',tline);
    end
    fclose(computeId);

    fprintf(fid,'Exit_Workspace\r\n');
    fprintf(fid,';\r\n');

    status = fclose(fid);
else
    status = -1;
    fprintf('Error opening pipeline file %s',strFNPipeline);
end

end

%% All Files
if allFiles == 1
    trialCell = {'pws','fast','slow','split','post'};
    strFNPipeline = sprintf('%s_%s.v3s',folder,'ALL');
    fid = fopen(fullfile(path,folder,strFNPipeline),'w');
    
if fid > 0
    % Create a new workspace (File|New).
    fprintf(fid,'File_New\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Create a hybrid model (Model|Create|Hybrid Model from C3DFile).
    fprintf(fid,'Create_Hybrid_Model\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'! /RANGE=ALL_FRAMES\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Apply a model template (Model|Apply Model Template).
    fprintf(fid,'Apply_Model_Template\r\n');
    fprintf(fid,'/MODEL_TEMPLATE=%s\r\n',fullfile(path,folder,strFNTemplate));
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Build model
    fprintf(fid,'Build_Model\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'! /REBUILD_ALL_MODELS=FALSE\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');
    
    for k = 1:length(trialCell)
    strFNMAT = sprintf('%s_%s.mat',trialCell{k},folder);
        switch trialCell{k}
            case 'pws'
            strFNMotionC3D = dirSBDnames{idx.c3d.pws};
            case 'fast'
            strFNMotionC3D = dirSBDnames{idx.c3d.fast};
            case 'slow'
            strFNMotionC3D = dirSBDnames{idx.c3d.slow};
            case 'split'
            strFNMotionC3D = dirSBDnames{idx.c3d.split};
            case 'post'
            strFNMotionC3D = dirSBDnames{idx.c3d.post};
        end
 
    % Insert a motion file (File|Open).
    fprintf(fid,'Open_File\r\n');
    fprintf(fid,'/FILE_NAME=%s\r\n',fullfile(path,folder,strFNMotionC3D));
    fprintf(fid,'! /CALIBRATION_FILE=\r\n');
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');

    % Assign the model to the motion file (Model|Assign Model to Motion
    % files).
    fprintf(fid,'Assign_Model_File\r\n');
    fprintf(fid,'/CALIBRATION_FILE=%s\r\n',fullfile(path,folder,strFNModelC3D));
    fprintf(fid,'/MOTION_FILE_NAMES=%s\r\n',fullfile(path,folder,strFNMotionC3D));
    fprintf(fid,';\r\n');
    fprintf(fid,'\r\n');
    
    end
    

    status = fclose(fid);
    
else
    status = -1;
    fprintf('Error opening pipeline file %s',strFNPipeline);

end
end
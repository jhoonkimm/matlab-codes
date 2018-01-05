clear all; close all; clc;

% Kinetic Analysis for Jay's study
% January 5th 2018
%% ---------------------------------------------------------------------------------------------------------------------------------------%
system = 'WIN';

if strcmp(system,'OS')
    path  = '\\Mac\Home\Google Drive\3-Research\2016\Trials';
    path2 = '\\Mac\Home\Google Drive\3-Research';
    path3 = '\\Mac\Home\Google Drive\3-Research\PaperComplianceStudy\FigureData';
elseif strcmp(system,'WIN')
    path  = 'E:\RESEARCH\MDSC508';
end

exclude = {'pws','fast','slow'};
excludeFolder = {'SBD_004','SBD_005','SBD_006','SBD_007','SBD_008','SBD_009','SBD_010','SBD_011','SBD_012'};

folder = {};
pathDir = dir(path);
for i = 1:length(pathDir)
    if pathDir(i).isdir
        if regexp(lower(pathDir(i).name),'sbd')
            folder{i} = pathDir(i).name;
        end
    end
end

%% --------------------------------------------------------------Executing V3D matfiles -------------------------------------------------------------------%
for  m = find(~cellfun(@isempty,folder))
    if ~strcmp(excludeFolder,folder{m})
        folderDir = dir(fullfile(path,folder{m}));
        matfiles = {};
        for i = 1:length(folderDir)
            if regexp(lower(folderDir(i).name),'processed')
                matfiles{i} = folderDir(i).name;
            end
        end
        dirSBDnames = {folderDir.name};
        idx = indexProcessing(path,folder{m});
        for  i = find(~cellfun(@isempty,matfiles))    
            fprintf('Reading matfile %s...\r\n',matfiles{i})
            tempName = matfiles{i};
            tempStruct = load(char(fullfile(path,folder{m},tempName)));
            
                    tempStruct.l_ank_power  = tempStruct.l_ank_power{1};
                    tempStruct.r_ank_power  = tempStruct.r_ank_power{1};
                    tempStruct.l_kne_power  = tempStruct.l_kne_power{1};
                    tempStruct.r_kne_power  = tempStruct.r_kne_power{1};
                    tempStruct.l_hip_power  = tempStruct.l_hip_power{1};
                    tempStruct.r_hip_power  = tempStruct.r_hip_power{1};
                    tempStruct.l_ank_angle  = tempStruct.l_ank_angle{1};
                    tempStruct.r_ank_angle  = tempStruct.r_ank_angle{1};
                    tempStruct.l_kne_angle  = tempStruct.l_kne_angle{1};
                    tempStruct.r_kne_angle  = tempStruct.r_kne_angle{1};
                    tempStruct.l_hip_angle  = tempStruct.l_hip_angle{1};
                    tempStruct.r_hip_angle  = tempStruct.r_hip_angle{1};
                    tempStruct.l_ank_moment = tempStruct.l_ank_moment{1};
                    tempStruct.r_ank_moment = tempStruct.r_ank_moment{1};
                    tempStruct.l_kne_moment = tempStruct.l_kne_moment{1};
                    tempStruct.r_kne_moment = tempStruct.r_kne_moment{1};
                    tempStruct.l_hip_moment = tempStruct.l_hip_moment{1};
                    tempStruct.r_hip_moment = tempStruct.r_hip_moment{1};
                    tempStruct.l_ank_force  = tempStruct.l_ank_force{1};
                    tempStruct.r_ank_force  = tempStruct.r_ank_force{1};
                    tempStruct.l_kne_force  = tempStruct.l_kne_force{1};
                    tempStruct.r_kne_force  = tempStruct.r_kne_force{1};
                    tempStruct.l_hip_force  = tempStruct.l_hip_force{1};
                    tempStruct.r_hip_force  = tempStruct.r_hip_force{1};
                    tempStruct.l_heel       = tempStruct.l_heel{1};
                    tempStruct.r_heel       = tempStruct.r_heel{1};
                    tempStruct.l_cop        = tempStruct.l_cop{1};
                    tempStruct.r_cop        = tempStruct.r_cop{1};
                                    
            kinetics.(tempName(11:17)).(tempName(19:end-4)) = tempStruct;
            
            if strcmp(lower(tempName(19:end-4)),'split')
                kinetics.(tempName(11:17)).(tempName(19:end-4)).treadmill = loadForcesFromHBCLBertecTreadmillMatFile(fullfile(path,folder{m},dirSBDnames{idx.mat.split}));
            elseif strcmp(lower(tempName(19:end-4)),'pws')
                kinetics.(tempName(11:17)).(tempName(19:end-4)).treadmill = loadForcesFromHBCLBertecTreadmillMatFile(fullfile(path,folder{m},dirSBDnames{idx.mat.pws}));
            elseif strcmp(lower(tempName(19:end-4)),'fast')
                kinetics.(tempName(11:17)).(tempName(19:end-4)).treadmill = loadForcesFromHBCLBertecTreadmillMatFile(fullfile(path,folder{m},dirSBDnames{idx.mat.fast}));
            elseif strcmp(lower(tempName(19:end-4)),'slow')
                kinetics.(tempName(11:17)).(tempName(19:end-4)).treadmill = loadForcesFromHBCLBertecTreadmillMatFile(fullfile(path,folder{m},dirSBDnames{idx.mat.slow}));
            elseif strcmp(lower(tempName(19:end-4)),'post')
                kinetics.(tempName(11:17)).(tempName(19:end-4)).treadmill = loadForcesFromHBCLBertecTreadmillMatFile(fullfile(path,folder{m},dirSBDnames{idx.mat.post}));
            end
        end
    end
end
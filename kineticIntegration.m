clear all; close all; clc;

% Kinetic Analysis for Jay's study
% November 3rd 2017
%% ---------------------------------------------------------------------------------------------------------------------------------------%
system = 'WIN';

if strcmp(system,'OS')
    path  = '\\Mac\Home\Google Drive\3-Research\2016\Trials';
    path2 = '\\Mac\Home\Google Drive\3-Research';
    path3 = '\\Mac\Home\Google Drive\3-Research\PaperComplianceStudy\FigureData';
elseif strcmp(system,'WIN')
    path  = 'E:\RESEARCH\MDSC508';
end

all = 1;
exclude = {'pws','fast','slow'};
excludeFolder = {'SBD_003','SBD_005','SBD_006','SBD_007','SBD_008','SBD_009','SBD_010','SBD_011','SBD_012'};

if all
folder = {};
pathDir = dir(path);
for i = 1:length(pathDir)
    if pathDir(i).isdir
        if regexp(lower(pathDir(i).name),'sbd')
            folder{i} = pathDir(i).name;
        end
    end
end
elseif ~all
    folder = {'SBD_009'};
end

%% -------------------------------------------------------------------- Integrating -------------------------------------------------------------------%

for  m = find(~cellfun(@isempty,folder))
    fprintf('Analyzing %s...\r\n',folder{m})
end
    fprintf('\n');
    fprintf('--------------------Starting integration----------------------------\n');
    fprintf('\n');

for m = 1:length(folder)
    integrationC3dAndTreadmill(path,folder{m},'exclude',exclude)
end

    fprintf('\n');
    fprintf('--------------------Ending integration----------------------------\n');
    fprintf('\n');
%% --------------------------------------------------------------Generating V3D pipelines -------------------------------------------------------------------%

for  m = find(~cellfun(@isempty,folder))
    if ~strcmp(excludeFolder,folder{m})
        fprintf('Generating pipelines for %s...\r\n',folder{m})
        V3DpipelineGenerator(path,folder{m})
    end
end

%% --------------------------------------------------------------Executing V3D pipelines -------------------------------------------------------------------%

for  m = find(~cellfun(@isempty,folder))
    if ~strcmp(excludeFolder,folder{m})
        folderDir = dir(fullfile(path,folder{m}));
        pipeline = {};
        for i = 1:length(folderDir)
            if regexp(lower(folderDir(i).name),'.v3s')
                pipeline{i} = folderDir(i).name;
            end
        end
        for  i = find(~cellfun(@isempty,pipeline))    
            fprintf('Executing pipeline %s...\r\n',pipeline{i})
            status = runV3DPipeline(char(fullfile(path,folder{m},pipeline(i))));
        end
    end
end


clc; close all; clear all;

% metabolic for Jay's study
% December 1st 2017
%---------------------------------------------------------------------------------------------------------------------------------------%
system = 'WIN';

if strcmp(system,'OS')
    path  = '\\Mac\Home\Google Drive\3-Research\2016\Trials';
    path2 = '\\Mac\Home\Google Drive\3-Research';
    path3 = '\\Mac\Home\Google Drive\3-Research\PaperComplianceStudy\FigureData';
elseif strcmp(system,'WIN')
    path  = 'E:\RESEARCH\MDSC 508';
end
%---------------------------------------------------------------------------------------------------------------------------------------%

folder = 'SBD_003';

%---------------------------------------------------------------------------------------------------------------------------------------%

dirSBD = dir(fullfile(path,folder));
dirSBDnames = {dirSBD.name};
idxDirSBD = find(not(cellfun('isempty',regexp(dirSBDnames,'.xlsx'))));

idxCPET = []; idxRest = [];
for i = idxDirSBD
    if regexp(lower(dirSBDnames{i}),'cpet')
    idxCPET =    [idxCPET,   i];
    elseif regexp(lower(dirSBDnames{i}),'ree')
    idxRest =    [idxRest,   i];
    end
end

metRest = cosmedProcessing(fullfile(path,folder,dirSBDnames{idxRest}));
calRest = caloricEquivalent(metRest.vo2kg(metRest.time>180),metRest.rer(metRest.time>180));

for k = 1:length(idxCPET)
    cpet{k,1} = dirSBDnames{idxCPET(k)};
end
sortedCPET = sort(cpet);

metPWS  = cosmedProcessing(fullfile(path,folder,sortedCPET{1}));
met150  = cosmedProcessing(fullfile(path,folder,sortedCPET{2}));
met75   = cosmedProcessing(fullfile(path,folder,sortedCPET{3}));
metSplt = cosmedProcessing(fullfile(path,folder,sortedCPET{4}),'timeConstraint',[1*60,16*60]);
metPost = cosmedProcessing(fullfile(path,folder,sortedCPET{4}),'timeConstraint',[16*60,26*60]);

idxSplitE = ((metSplt.time>240) + (metSplt.time<360))==2;
idxPostE  = ((metPost.time>1140) + (metPost.time<1200))==2;


calPWS = caloricEquivalent(smooth(metPWS.vo2kg(metPWS.time>180),5,'moving'),metPWS.rer(metPWS.time>180));
cal150 = caloricEquivalent(smooth(met150.vo2kg(met150.time>180),5,'moving'),met150.rer(met150.time>180));
cal75  = caloricEquivalent(smooth(met75.vo2kg(met75.time>180),5,'moving'),met75.rer(met75.time>180));
calSpltE  = caloricEquivalent(smooth(metSplt.vo2kg(idxSplitE),5,'moving'),metSplt.rer(idxSplitE));
calSpltL  = caloricEquivalent(smooth(metSplt.vo2kg(metSplt.time>840),5,'moving'),metSplt.rer(metSplt.time>840));
calPostE  = caloricEquivalent(smooth(metPost.vo2kg(idxPostE),5,'moving'),metPost.rer(idxPostE));
calPostL  = caloricEquivalent(smooth(metPost.vo2kg(metPost.time>1440),5,'moving'),metPost.rer(metPost.time>1440));

[mean(calPWS)
mean(cal150 )
mean(cal75  )
mean(calSpltE)
mean(calSpltL )
mean(calPostE )
mean(calPostL)]





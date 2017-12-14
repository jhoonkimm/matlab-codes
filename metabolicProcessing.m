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

pathDir = dir(path);

folder = {};
for i = 1:length(pathDir)
    if pathDir(i).isdir
        if regexp(lower(pathDir(i).name),'sbd')
            folder{i} = pathDir(i).name;
        end
    end
end

%---------------------------------------------------------------------------------------------------------------------------------------%
%% Metabolic Processing

met = struct;
metcell = [];
spltcell = [];

for m = find(~cellfun(@isempty,folder))
    clear dirSBDnames
    dirSBD = dir(fullfile(path,folder{m}));
    dirSBDnames = {dirSBD.name};
    idxDirSBD = find(not(cellfun('isempty',regexp(dirSBDnames,'.xlsx'))));

idxCPET = []; idxRest = [];
for i = idxDirSBD
    if regexp(lower(dirSBDnames{i}),'cpet')
        idxCPET =    [idxCPET,   i];
    elseif regexp(lower(dirSBDnames{i}),'ree') %cosmed #2
        idxRest =    [idxRest,   i];
    elseif regexp(lower(dirSBDnames{i}),'rest') %cosmed #1
        idxRest =    [idxRest,   i];
    end
end

if regexp(lower(dirSBDnames{idxRest}),'rest')
    metRest = cosmedProcessing(fullfile(path,folder{m},dirSBDnames{idxRest}),'cosmedOne',true);
elseif regexp(lower(dirSBDnames{idxRest}),'ree')
    metRest = cosmedProcessing(fullfile(path,folder{m},dirSBDnames{idxRest}));
end
    calRest = caloricEquivalent(metRest.vo2kg(metRest.time>180),metRest.rer(metRest.time>180));

for k = 1:length(idxCPET)
    cpet{k,1} = dirSBDnames{idxCPET(k)};
end
sortedCPET = sort(cpet);

if regexp(lower(dirSBDnames{idxRest}),'rest')
    metPWS  = cosmedProcessing(fullfile(path,folder{m},sortedCPET{1}),'cosmedOne',true);
    met150  = cosmedProcessing(fullfile(path,folder{m},sortedCPET{2}),'cosmedOne',true);
    met75   = cosmedProcessing(fullfile(path,folder{m},sortedCPET{3}),'cosmedOne',true);
    if strcmp(folder{m},'SBD_009')
        metSplt = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[1*60,18*60],'cosmedOne',true);
        metPost = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[18*60,28*60],'cosmedOne',true);
    else
        metSplt = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[1*60,16*60],'cosmedOne',true);
        metPost = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[16*60,26*60],'cosmedOne',true);
    end
elseif regexp(lower(dirSBDnames{idxRest}),'ree')
    metPWS  = cosmedProcessing(fullfile(path,folder{m},sortedCPET{1}));
    met150  = cosmedProcessing(fullfile(path,folder{m},sortedCPET{2}));
    met75   = cosmedProcessing(fullfile(path,folder{m},sortedCPET{3}));
    if strcmp(folder{m},'SBD_009')
        metSplt = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[1*60,18*60]);
        metPost = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[18*60,28*60]);
    else
        metSplt = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[1*60,16*60]);
        metPost = cosmedProcessing(fullfile(path,folder{m},sortedCPET{4}),'timeConstraint',[16*60,26*60]);
    end
end



idxSplitE = ((metSplt.time>240) + (metSplt.time<360))==2;
idxPostE  = ((metPost.time>1140) + (metPost.time<1260))==2;
idxSplitL = ((metSplt.time>540) + (metSplt.time<660))==2;
idxPostL  = ((metPost.time>1440) + (metPost.time<1560))==2;

calPWS = caloricEquivalent(smooth(metPWS.vo2kg(metPWS.time>180),5,'moving'),metPWS.rer(metPWS.time>180));
cal150 = caloricEquivalent(smooth(met150.vo2kg(met150.time>180),5,'moving'),met150.rer(met150.time>180));
cal75  = caloricEquivalent(smooth(met75.vo2kg(met75.time>180),5,'moving'),met75.rer(met75.time>180));
calSpltE  = caloricEquivalent(smooth(metSplt.vo2kg(idxSplitE),5,'moving'),metSplt.rer(idxSplitE));
calSpltL  = caloricEquivalent(smooth(metSplt.vo2kg(idxSplitL),5,'moving'),metSplt.rer(idxSplitL));
calPostE  = caloricEquivalent(smooth(metPost.vo2kg(idxPostE),5,'moving'),metPost.rer(idxPostE));
calPostL  = caloricEquivalent(smooth(metPost.vo2kg(idxPostL),5,'moving'),metPost.rer(idxPostL));

spltMet = [];
per = [4,5,6,7,8,9,10,11,12,13,14,15,16]*60;
for k = 1:length(per)-1
    idxPer  = ((metSplt.time>per(k)) + (metSplt.time<per(k+1)))==2;
    spltMet = [spltMet;mean(caloricEquivalent(smooth(metSplt.vo2kg(idxPer),5,'moving'),metSplt.rer(idxPer)))-mean(calRest)];
end

met.(folder{m}) = [ mean(calPWS)-mean(calRest)
                    mean(cal150 )-mean(calRest)
                    mean(cal75  )-mean(calRest)
                    mean(calSpltE)-mean(calRest)
                    mean(calSpltL )-mean(calRest)
                    mean(calPostE )-mean(calRest)
                    mean(calPostL)-mean(calRest)];

metcell = [metcell, met.(folder{m})];
spltcell = [spltcell, spltMet];
end

%%
metcelli = metcell';

sizeRaw = size(metcelli);  repAn = [];
    for k = 1:sizeRaw(2)
        for i = 1:sizeRaw(1)
            repAn = [repAn;metcelli(i,k) k i];
        end
    end
p = rmaov1(repAn);

[pvalues, results] = hbclstats(metcelli);

%%

plot(mean(spltcell,2))
hold on
errorbar(1:length(mean(spltcell,2)),mean(spltcell,2),std(spltcell'))
hold off



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
    path  = 'E:\RESEARCH\MDSC508';
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



idxSplitE = ((metSplt.time>4*60) + (metSplt.time<6*60))==2;
idxPostE  = ((metPost.time>19*60) + (metPost.time<21*60))==2;
idxSplitL = ((metSplt.time>14*60) + (metSplt.time<16*60))==2;
idxPostL  = ((metPost.time>24*60) + (metPost.time<26*60))==2;

calPWS = caloricEquivalent(smooth(metPWS.vo2kg(metPWS.time>180),5,'moving'),metPWS.rer(metPWS.time>180));
cal150 = caloricEquivalent(smooth(met150.vo2kg(met150.time>180),5,'moving'),met150.rer(met150.time>180));
cal75  = caloricEquivalent(smooth(met75.vo2kg(met75.time>180),5,'moving'),met75.rer(met75.time>180));
calSpltE  = caloricEquivalent(smooth(metSplt.vo2kg(idxSplitE),5,'moving'),metSplt.rer(idxSplitE));
calSpltL  = caloricEquivalent(smooth(metSplt.vo2kg(idxSplitL),5,'moving'),metSplt.rer(idxSplitL));
calPostE  = caloricEquivalent(smooth(metPost.vo2kg(idxPostE),5,'moving'),metPost.rer(idxPostE));
calPostL  = caloricEquivalent(smooth(metPost.vo2kg(idxPostL),5,'moving'),metPost.rer(idxPostL));

met.(folder{m}) = [ mean(calPWS)-mean(calRest)
                    mean(cal150 )-mean(calRest)
                    mean(cal75  )-mean(calRest)
                    mean(calSpltE)-mean(calRest)
                    mean(calSpltL )-mean(calRest)
                    mean(calPostE )-mean(calRest)
                    mean(calPostL)-mean(calRest)];

metcell = [metcell, met.(folder{m})];

end

%% Statistics
metcelli = metcell';

sizeRaw = size(metcelli);  repAn = [];
    for k = 1:sizeRaw(2)
        for i = 1:sizeRaw(1)
            repAn = [repAn;metcelli(i,k) k i];
        end
    end
p = rmaov1(repAn);

[pvalues, results] = hbclstats(metcelli);
%% Labeling Statistic Results
resultcell = {};
for m = 1:size(results,1)
    for k = 1:size(results,2)
        if k<5
            if results(m,k)==1
                resultcell{m,k}='PWS';
            elseif results(m,k)==2
                resultcell{m,k}='150';
            elseif results(m,k)==3
                resultcell{m,k}='75';
            elseif results(m,k)==4
                resultcell{m,k}='Early adaptation';
            elseif results(m,k)==5
                resultcell{m,k}='Late adaptation';
            elseif results(m,k)==6
                resultcell{m,k}='Early postadaptation';
            elseif results(m,k)==7
                resultcell{m,k}='Late postadaptation';
            else
                resultcell{m,k}=results(m,k);
            end
        else
            resultcell{m,k}=results(m,k);
        end
    end
end
%% Bar graphs
%     Xlabels = {'slow baseline','fast baseline','early adaptation','late adaptation','early postadaptation','late postadaptation'};
    Xlabels = {'NB','FB','SB','EA','LA','EP','LP'};
    x = metcelli;
    stdx = std(metcelli);
    for i = 1:length(mean(x))
        sem = stdx(i)/sqrt(length(mean(x)));
        t = tinv(0.975,length(mean(x))-1);
        CI(i) = (t*sem);
    end
    
    figure('Units', 'pixels', ...
           'Position', [500 500 650 450]);
    hold on;

     hBar = bar(mean(x));
     hE = errorbar(1:7,mean(x),CI);
     
     set(hBar,                              ...
        'FaceColor',       'black',         ...
        'LineWidth',       0.5,             ...
        'BarWidth',        0.5,             ...
        'FaceAlpha',       0.35)
     set(hE                            ,    ...
        'LineWidth'       , 1           ,   ...
        'Marker'          , 'o'         ,   ...
        'MarkerSize'      , 6           ,   ...
        'MarkerEdgeColor' , [.2 .2 .2]  ,   ...
        'MarkerFaceColor' , [.7 .7 .7]  ,   ...
        'LineStyle'       , 'none'      ,   ...
        'Color'           , 'black');
     set(gca,...
        'XTick',        1:7,            ...
        'XTickLabels',  Xlabels,        ...
        'YLim',         [1.5,5],        ...
        'FontName',     'Times',        ...
        'FontSize',     10,             ...
        'Box',          'off',          ...
        'TickDir',      'out',          ...
        'TickLength',   [.02 .02],      ...
        'YMinorTick',   'on',           ...
        'YGrid',        'on',           ...
        'XColor',       [0 0 0],     ...
        'YColor',       [0 0 0],     ...
        'LineWidth',    1               );
    
    hBarXLabel = xlabel('Conditions',                       ...
                        'FontSize',     12,                 ...
                        'Position',     [3.5 1.2 -1]      );
                    
    hBarYLabel = ylabel('Net metabolic power (W kg^{-1})',  ...
                        'FontSize',     12,                 ...
                        'Position',     [-0.75 3.25 -1]      );
                    
                    hold off

clear all; clc; close all;

path  = 'E:\RESEARCH\MDSC 508';
pathDir = dir(path);

folder = {};
for i = 1:length(pathDir)
    if pathDir(i).isdir
        if regexp(lower(pathDir(i).name),'sbd')
            folder{i} = pathDir(i).name;
        end
    end
end

age = [];
for m = find(~cellfun(@isempty,folder))
    output = subjectInfoLoading(fullfile(path,folder{m},sprintf('%s.xlsx',folder{m})));
    age = [age;output.age];
end


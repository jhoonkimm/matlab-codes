function [output] = subjectInfoLoading(xlsFilename,varargin)
%  
% [output] = subjectInfoLoading(xlsFilename,varargin)
%  Allow for processing of subject files.
%  xlsFilename must be .txt file.
%
%  Flags:
%
%
%  Created by Jay Kim (12/14/2017)

%% Default input arguments


%% Optional input arguments

% for i = 1:2:length(varargin)
%   switch varargin{i}
%     case 'brockway'
%     brockway = varargin{i + 1};
%     case 'smooth'
%     filt = varargin{i + 1};
%     case 'massSpecific'
%     massSpecific = varargin{i + 1};
%   end
% end

%%
output = struct;

[~,~,raw] = xlsread(xlsFilename);

output.date         = raw{1,2};
output.time         = raw{1,4};
output.id           = raw{2,2};
output.fname        = raw{3,2};
output.lname        = raw{4,2};
output.sex          = raw{5,2};
output.dob          = raw{6,2};
output.height       = raw{7,2};
output.weight       = raw{8,2};
output.leg          = raw{9,2};
output.domLeg       = raw{10,2};
output.nonDomLeg    = raw{10,4};
output.speedPWS     = raw{12,2};
output.speed150     = raw{13,2};
output.speed075     = raw{14,2};
output.notes        = raw{16,2};
output.rightheel    = raw{1,8};
output.leftheel     = raw{2,8};
output.raw          = raw;

age = datevec(datenum(output.date,'yyyy-mm-dd')-datenum(output.dob,'yyyy-mm-dd'));
output.age = age(1,1:3);



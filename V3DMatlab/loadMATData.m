function status = loadMATData(strFNData)
%loadMATData(strFNData)   This function demonstrates how to access the
%   fields of a user specified MAT file, which could be produced by the
%   Visual3D export to MATLAB pipeline commands.
%
%   strFNData:  a MAT file containing data in which to load.
%
%   Example:  status = loadMATData('C:\data\example.mat');
%
%   Copyright (C) 2008 C-Motion, Inc.  All rights reserved.
%   For support contact support@c-motion.com.  

disp(['loadMATData() -- Loading data from ' strFNData '...']);

status = 0;

% Load the MAT file into memory.
matData = load(strFNData);

% Get a list of the field names in the file.
strFieldNames = fieldnames(matData);

% For each field, display its name and corresponding data (if appropriate).
if numel(strFieldNames) > 0
    disp('Field Names:');
    disp(strFieldNames);

    for fn = 1:numel(strFieldNames)
        disp(strFieldNames{fn});
        % Get and display the data associated with each field.
        tmp = getfield(matData,strFieldNames{fn});
       
        if numel(tmp) > 0
            data = tmp{1};
            
            % Display the data according to its type.
            if ischar(data)
                disp(data);
                disp(' ');
            else
                disp(data);
            end
        else
            disp(['Field ' strFieldNames{fn} ' data missing.']);
        end
    end
else
    status = -1;
    disp(['MAT file ' strFNData ' contains no fields.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

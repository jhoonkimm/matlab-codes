function status = runV3DPipeline(strFNPipeline)
%runV3DPipeline(strFNPipeline)   This function demonstrates how to run
%   Visual3D and execute a user specified pipeline script.
%
%   strFNPipeline:  the filename of a Visual3D pipeline to execute.  If
%       strFNPipeline is not a full path then the location of the pipeline
%       file must be in the path environment variable.
%
%   Example:  status = runV3DPipeline('C:\data\example.v3s');
%
%   Copyright (C) 2008 C-Motion, Inc.  All rights reserved.
%   For support contact support@c-motion.com.  

disp(['runV3DPipeline() -- Running Visual3D with pipeline ' strFNPipeline '...']);

% Important:  Visual3D.exe must be in the user's path in order for the
% MATLAB 'dos' command to work properly.  Alternatively, the full path to
% the application can be used instead
% (e.g., '"C:\Program Files\Visual3D v3\Visual3D.exe" /s' strFNPipeline).
% [status,result] = dos(['Visual3D.exe /s "' strFNPipeline]);
[status,result] = dos(['"C:\Program Files\Visual3D v6 x64\Visual3D.exe" /s ' strFNPipeline '&'],'-echo');

if status == 0
    fprintf('Visual3D executed properly')
end

% To edit the path on Windows XP:
% 1.  Click Start menu, then right click My Computer and select Properties.
% 2.  Click the Advanced tab and then the Environment Variables button.
% 3.  The path variable is in the System variables section.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

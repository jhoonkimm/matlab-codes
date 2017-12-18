function [inertiaPipelinePath,typestr,folderstr,resultname,varname,status] = generateV3DInertiaPipeline(inertiaPipelinePath,segments)

typestr={};
folderstr={};
resultname={};
varname={};
% opt_argin = varargin;
fid=fopen(inertiaPipelinePath,'w');

if fid > 0
  while length(segments) >= 1,
    val = segments{1};
    segments = segments(2:end);
    switch lower(val)
      case {'rft','rsk','rth','rar','rfa','lft','lsk','lth','lar','lfa','rpv','rta'}
        fprintf(fid,'Evaluate_Expression\r\n');
        fprintf(fid,'/EXPRESSION=MODEL::SEGMENT::%s::IXX\r\n',upper(val));
        fprintf(fid,'/RESULT_NAME=%s_Ixx\r\n',upper(val));
        fprintf(fid,'/RESULT_TYPE=DERIVED\r\n');
        fprintf(fid,'/RESULT_FOLDER=Mass_Inertia\r\n');
        fprintf(fid,';\r\n');
        fprintf(fid,'\r\n');
        resultname = vertcat(resultname,sprintf('%s_Ixx',upper(val)));
        varname = vertcat(varname,sprintf('%s_Ixx',lower(val)));
        
        fprintf(fid,'Evaluate_Expression\r\n');
        fprintf(fid,'/EXPRESSION=MODEL::SEGMENT::%s::IXX\r\n',upper(val));
        fprintf(fid,'/RESULT_NAME=%s_Iyy\r\n',upper(val));
        fprintf(fid,'/RESULT_TYPE=DERIVED\r\n');
        fprintf(fid,'/RESULT_FOLDER=Mass_Inertia\r\n');
        fprintf(fid,';\r\n');
        fprintf(fid,'\r\n');
        resultname = vertcat(resultname,sprintf('%s_Iyy',upper(val)));
        varname = vertcat(varname,sprintf('%s_Iyy',lower(val)));
        
        fprintf(fid,'Evaluate_Expression\r\n');
        fprintf(fid,'/EXPRESSION=MODEL::SEGMENT::%s::IZZ\r\n',upper(val));
        fprintf(fid,'/RESULT_NAME=%s_Izz\r\n',upper(val));
        fprintf(fid,'/RESULT_TYPE=DERIVED\r\n');
        fprintf(fid,'/RESULT_FOLDER=Mass_Inertia\r\n');
        fprintf(fid,';\r\n');
        fprintf(fid,'\r\n');
        resultname = vertcat(resultname,sprintf('%s_Izz',upper(val)));
        varname = vertcat(varname,sprintf('%s_Izz',lower(val)));

        fprintf(fid,'Evaluate_Expression\r\n');
        fprintf(fid,'/EXPRESSION=MODEL::SEGMENT::%s::MASS\r\n',upper(val));
        fprintf(fid,'/RESULT_NAME=%s_Mass\r\n',upper(val));
        fprintf(fid,'/RESULT_TYPE=DERIVED\r\n');
        fprintf(fid,'/RESULT_FOLDER=Mass_Inertia\r\n');
        fprintf(fid,';\r\n');
        fprintf(fid,'\r\n');
        resultname = vertcat(resultname,sprintf('%s_Mass',upper(val)));
        varname = vertcat(varname,sprintf('%s_Mass',lower(val)));
        
        % output rotation matrix
        fprintf(fid,'Compute_Model_Based_Data\r\n');
        fprintf(fid,'/RESULT_NAME=%s_Rotmat\r\n',upper(val));
        fprintf(fid,'/FUNCTION=JOINT_ROTATION\r\n');
        fprintf(fid,'/SEGMENT=%s\r\n',upper(val));
        fprintf(fid,'! /REFERENCE_SEGMENT=LAB\r\n');
        fprintf(fid,'! /RESOLUTION_COORDINATE_SYSTEM=LAB\r\n');
        fprintf(fid,'! /USE_CARDAN_SEQUENCE=FALSE\r\n');
        fprintf(fid,'! /NORMALIZATION=FALSE\r\n');
        fprintf(fid,'! /NORMALIZATION_METHOD=\r\n');
        fprintf(fid,'! /NORMALIZATION_METRIC=\r\n');
        fprintf(fid,'! /NEGATEX=FALSE\r\n');
        fprintf(fid,'! /NEGATEY=FALSE\r\n');
        fprintf(fid,'! /NEGATEZ=FALSE\r\n');
        fprintf(fid,'! /AXIS1=X\r\n');
        fprintf(fid,'! /AXIS2=Y\r\n');
        fprintf(fid,'! /AXIS3=Z\r\n');
        fprintf(fid,'! /TREADMILL_DATA=FALSE\r\n');
        fprintf(fid,'/TREADMILL_DIRECTION=UNIT_VECTOR(0, 1, 0)\r\n');
        fprintf(fid,'! /TREADMILL_SPEEED=0.0\r\n');
        fprintf(fid,';\r\n');
        fprintf(fid,'\r\n');
        resultname = vertcat(resultname,sprintf('%s_Rotmat',upper(val)));
        varname = vertcat(varname,sprintf('%s_Rotmat',lower(val)));
        
        typestr=[vertcat(typestr,repmat({'DERIVED'},4,1)); 'LINK_MODEL_BASED'];
        folderstr=[vertcat(folderstr,repmat({'Mass_Inertia'},4,1)); 'Original'];
        
      otherwise
        warning('\nWarning: unrecognized segment ''%s''. Segment will be ignored.\n',val)
    end
  end
  status = fclose(fid);
else
  status = -1;
  disp(['Error opening pipeline file ' pipelinefile]);
end
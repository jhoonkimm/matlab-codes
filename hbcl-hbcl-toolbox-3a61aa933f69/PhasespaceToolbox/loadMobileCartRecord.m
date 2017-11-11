function [mobileCartRecord] = loadMobileCartRecord(mobileCartRecordFilename, varargin)

verbose = 0;
forceReprocess = 0;
for i=1:length(varargin)
  if (strcmp(varargin{i}, 'verbose'))
    verbose = varargin{i+1};
  end
  if (strcmp(varargin{i}, 'forceReprocess'))
    forceReprocess = varargin{i+1};
  end
end

%%
outputFilename = createMatOutputFilename(mobileCartRecordFilename);

if (forceReprocess || ~exist(outputFilename, 'file'))
  %%
  cartSensorData = importdata(mobileCartRecordFilename);
  mobileCartRecord.packetNumber = cartSensorData(:,1);
  mobileCartRecord.timer = cartSensorData(:,2);
  mobileCartRecord.leftEncoder = cartSensorData(:,3);
  mobileCartRecord.rightEncoder = cartSensorData(:,4);
  mobileCartRecord.gyroIntegratedRadians = (cartSensorData(:,5) / 10) * (pi / 180);
  
  save(outputFilename, 'mobileCartRecord');
else
  load(outputFilename, 'mobileCartRecord');
end


if (verbose > 0)
  figure;
  subplot(2,1,1);
  plot(mobileCartRecord.packetNumber, mobileCartRecord.leftEncoder, 'r');
  hold on;
  plot(mobileCartRecord.packetNumber, mobileCartRecord.rightEncoder, 'g');
  legend('leftEncoderCounts', 'rightEncoderCounts');
  subplot(2,1,2);
  plot(mobileCartRecord.packetNumber, mobileCartRecord.gyroIntegratedRadians);
  legend('gyroIntegratedRadians');
end

end


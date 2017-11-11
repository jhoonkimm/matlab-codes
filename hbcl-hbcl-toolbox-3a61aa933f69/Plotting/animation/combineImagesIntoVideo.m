function [outputVideo] = combineImagesIntoVideo(imageFilenames, videoOutputFilename, frameRate, varargin)
%%

videoProfile = [];

for i = 1 : 2 : length(varargin)
  option = varargin{i};
  value = varargin{i + 1};
  switch(option)
    case 'videoProfile'
    videoProfile = value;
  end
end


%%

if (isempty(videoProfile))
  outputVideo = VideoWriter(videoOutputFilename);
else
  outputVideo = VideoWriter(videoOutputFilename, videoProfile);
end
outputVideo.FrameRate = frameRate;
open(outputVideo);

for i = 1 : length(imageFilenames)
  img = imread(imageFilenames{i});
  writeVideo(outputVideo, img);
end

close(outputVideo);

end


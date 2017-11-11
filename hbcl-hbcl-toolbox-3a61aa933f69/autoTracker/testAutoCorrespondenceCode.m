
%% test our code for figuring out, from a cloud of points,
% which of those points correspond to a given rigid body

clear all;
close all;
format compact;

minSquaredErrorLimit = 0.2;

pointsInTemplate = [[0 0 0]', ...
  [1 0 0]', ...
  [0 2 0]', ...
  [0 0 3]'];

%%
templateBody = RigidBody;
for i = 1:size(pointsInTemplate, 2)
  templateBody = templateBody.addTrackerPositionInBodyFrame(pointsInTemplate(:, i));
end

% templateBody.position = [12 6 -1]';
% templateBody.orientation = [pi/8 -pi 3*pi/2]';

templateBody.position = [0 0 1]';
templateBody.orientation = [0 0 0]';

transformedTemplatePoints = templateBody.getTrackerPositionsInWorld();

templateBody.position = [0 0 0 ]';
templateBody.orientation = [0 0 0]';

% add a bad marker on the template
templateBody = templateBody.addTrackerPositionInBodyFrame([-3 -3 -3]');

%%
markerNoiseStandardDeviation = 0.05; %3;
cloudBody = RigidBody;

cloudBody = cloudBody.addTrackerPositionInBodyFrame([-2 2 1]');
cloudBody = cloudBody.addTrackerPositionInBodyFrame([1 0 2]');

noisyTransformedPoints = transformedTemplatePoints + ...
  randn(size(transformedTemplatePoints)) .* markerNoiseStandardDeviation;
% fail to add one of the template markers in the cloud body;
for i = 1:(size(transformedTemplatePoints, 2) - 1)
  cloudBody = cloudBody.addTrackerPositionInWorldFrame(noisyTransformedPoints(:, i));
end

% cloudBody = cloudBody.addTrackerPositionInBodyFrame([0 2 1]');
% cloudBody = cloudBody.addTrackerPositionInBodyFrame([1 0 3]');

% for i = 1:10
%   templateBody = templateBody.addTrackerPositionInBodyFrame([10 10 10]');
%   templateBody = templateBody.addTrackerPositionInBodyFrame([1 2 3]');
%   templateBody = templateBody.addTrackerPositionInBodyFrame([4 -6 5]');
% end

%%
figure();
subplot(2,2,[1]);
% templateBody.plotTheseWorldPoints(pointsInTemplate, {'ko'});
templateBody.plotTrackersInWorldStyled({'ko'});
cloudBody.plotTrackersInWorldStyled({'go'});
axis equal;

%%
[correspondeces] = templateBody.findCorrespondecesOfThisBodyToGivenPoints(cloudBody.trackersInBodyFrame);


%% now get a few template bodies, and add them all to the cloud, and figure them out again:

pointsInTemplates = {[[0 0 0]', ...
  [1 0 0]', ...
  [0 2 0]', ...
  [0 0 3]'], ...
  ...
  [[0 0 0]', ...
  [1 0 0]', ...
  [0 2 0]', ...
  [0 0 3]'] * 1.5, ...
  ...
  [[0 0 -1]', ...
  [1 0 0]', ...
  [0 2 0]', ...
  [0 0 -3]']};

cloudOffsets = {[0 0 0]', ...
  [0 0.3 0.1]', ...
  [-0.6 0 0]'};
cloudOrientation = {rand(3,1), ...
  rand(3,1), ...
  rand(3,1)};


templateBodies = {};
for i =  1:length(pointsInTemplates)
  templateBodies{i} = RigidBody;
  templateBodies{i}.trackersInBodyFrame = pointsInTemplates{i};
end

cloudBody = RigidBody;

for i = 1:length(templateBodies)
  templateBodies{i}.position = cloudOffsets{i};
  templateBodies{i}.orientation = cloudOrientation{i};
  
  cloudBody.trackersInBodyFrame = [cloudBody.trackersInBodyFrame templateBodies{i}.getTrackerPositionsInWorld];
  
  templateBodies{i}.position = [0 0 0]';
  templateBodies{i}.orientation = [0 0 0]';
end

figure;
cloudBody.plotTrackersInWorldStyled({'ro'})
axis equal;

%%
correspondeces = [];
numMarkersOverall = 0;
for i = 1:length(templateBodies)
  fprintf('finding correspondences for template body %g\n', i);
  
  [theseCorrepondeces] = ...
    templateBodies{i}.findCorrespondecesOfThisBodyToGivenPoints(cloudBody.trackersInBodyFrame, 'verbose', 1);
  
  theseCorrepondeces(:, 1) = theseCorrepondeces(:, 1) + numMarkersOverall;
  numMarkersOverall = numMarkersOverall + size(templateBodies{i}.trackersInBodyFrame, 2);
  
  correspondeces = [correspondeces; theseCorrepondeces];
end

correspondeces


%% now do a test with some common markers.

pointsInTemplateSet.M1 = [0 0 2]';
% pointsInTemplateSet.M2 = [-0.1 0 1.9]';
% pointsInTemplateSet.M3 = [0.1 0 1.9]';
% pointsInTemplateSet.M4 = [0 -0.2 1.7]';
% pointsInTemplateSet.M5 = [0 -0.1 1.6]';
% 
% rigidBodyMarkerNames.torso = {'M1', 'M2', 'M3', 'M4', 'M5'};

% pointsInTemplateSet.M6 = [-0.2 0.2 1.5]';
% pointsInTemplateSet.M7 = [0.2 0.2 1.5]';
% pointsInTemplateSet.M8 = [-0.18 -0.2 1.45]';
% pointsInTemplateSet.M9 = [0.18 -0.2 1.45]';
pointsInTemplateSet.M10LTrochanter = [0.25 0 1.0]';
% pointsInTemplateSet.M11RTrochanter = [-0.25 0 1.0]';
% 
% rigidBodyMarkerNames.pelvis = {'M6', 'M7', 'M8', 'M9', 'M10LTrochanter', 'M11RTrochanter'};

pointsInTemplateSet.M12 = [0.25 0.1 0.9]';
pointsInTemplateSet.M13 = [0.20 -0.1 0.85]';
pointsInTemplateSet.M14LKneeLateral = [0.2 0.0 0.5]';
pointsInTemplateSet.M15LKneeMedial = [0.05 0.0 0.5]';

rigidBodyMarkerNames.leftThigh = {'M10LTrochanter', 'M12', 'M13', 'M14LKneeLateral', 'M15LKneeMedial'};

% pointsInTemplateSet.M16 = [-0.25 0.1 0.9]';
% pointsInTemplateSet.M17 = [-0.20 -0.1 0.85]';
% pointsInTemplateSet.M18RKneeLateral = [-0.2 0.0 0.5]';
% pointsInTemplateSet.M19RKneeMedial = [-0.05 0.0 0.5]';
% 
% rigidBodyMarkerNames.rightThigh = {'M11RTrochanter', 'M16', 'M17', 'M18RKneeLateral', 'M19RKneeMedial'};

pointNames = fieldnames(pointsInTemplateSet);
pointCloud = zeros(3, length(pointNames));
for i = 1:length(pointNames)
  pointCloud(:, i) = pointsInTemplateSet.(pointNames{i});
end

noise = 0.001;
allTemplatePoints = pointCloud;
pointCloud = pointCloud + randn(size(pointCloud)) .* noise;


rigidBodyNames = fieldnames(rigidBodyMarkerNames);
templateBodies = [];
for i = 1:length(rigidBodyNames)
  rigidBodyName = rigidBodyNames{i};
  theseRigidBodyMarkerNames = rigidBodyMarkerNames.(rigidBodyName);
  templateBodies.(rigidBodyName) = RigidBody;
  
  for markerNumber = 1:length(theseRigidBodyMarkerNames)
    markerName = theseRigidBodyMarkerNames{markerNumber};
    templateBodies.(rigidBodyName) = ...
      templateBodies.(rigidBodyName).addTrackerPositionInBodyFrame(pointsInTemplateSet.(markerName));
  end
end


%
% for now, let's assume we already know a correspondence:
leftTrochanterIndex = find(strcmp(fieldnames(pointsInTemplateSet), 'M10LTrochanter'));

correspondeces = [leftTrochanterIndex leftTrochanterIndex 
];

rigidBodyName = 'leftThigh'; %rigidBodyNames{i};
fprintf('finding correspondences for template body %s\n', rigidBodyName);

thisTemplate = templateBodies.(rigidBodyName);

theseConstrainedCorrespondences = [];
theseMarkerNames = rigidBodyMarkerNames.(rigidBodyName);

for i = 1:length(theseMarkerNames)
  thisMarkerIndex = find(strcmp(fieldnames(pointsInTemplateSet), theseMarkerNames{i}));
  knownCorrespondenceIndex = find(correspondeces(:, 1) == thisMarkerIndex);
  if (~isempty(knownCorrespondenceIndex))
    theseConstrainedCorrespondences = [theseConstrainedCorrespondences; ...
      [i correspondeces(knownCorrespondenceIndex, 2)]];
  end
end

%%
[theseCorrepondeces] = ...
  thisTemplate.findCorrespondecesOfThisBodyToGivenPoints(pointCloud, ...
  'verbose', 1, ...
  'constrainedCorrespondences', theseConstrainedCorrespondences);

%%
% numMarkersOverall = 0;
% for i = 1:length(rigidBodyNames)
%   rigidBodyName = rigidBodyNames{i};
%   thisTemplate = templateBodies.(rigidBodyName);
%   
%   fprintf('finding correspondences for template body %s\n', rigidBodyName);
%   
%   [theseCorrepondeces] = ...
%     templateBodies{i}.findCorrespondecesOfThisBodyToGivenPoints(cloudBody.trackersInBodyFrame, ...
%     'verbose', 1, ...
%     'constrainedCorrespondences', []);
%   
%   theseCorrepondeces(:, 1) = theseCorrepondeces(:, 1) + numMarkersOverall;
%   numMarkersOverall = numMarkersOverall + size(templateBodies{i}.trackersInBodyFrame, 2);
%   
%   correspondeces = [correspondeces; theseCorrepondeces];
% end


%%
clear all;
%%
% clear MarkeredRigidBodyChain

pointsInTemplateSet.M1 = [0 0 2]';
pointsInTemplateSet.M2 = [-0.1 0 1.9]';
pointsInTemplateSet.M3 = [0.1 0 1.9]';
pointsInTemplateSet.M4 = [0 -0.2 1.7]';
pointsInTemplateSet.M5 = [0 -0.2 1.6]';

rigidBodyMarkerNames.torso = {'M1', 'M2', 'M3', 'M4', 'M5'};

pointsInTemplateSet.M6 = [-0.2 0.2 1.3]';
pointsInTemplateSet.M7 = [0.2 0.2 1.3]';
pointsInTemplateSet.M8 = [-0.18 -0.2 1.25]';
pointsInTemplateSet.M9 = [0.18 -0.2 1.25]';
pointsInTemplateSet.M10LTrochanter = [0.25 0 1.0]';
pointsInTemplateSet.M11RTrochanter = [-0.25 0 1.0]';

rigidBodyMarkerNames.pelvis = {'M6', 'M7', 'M8', 'M9', 'M10LTrochanter', 'M11RTrochanter'};
% rigidBodyMarkerNames.pelvis = {'M6', 'M10LTrochanter', 'M11RTrochanter'};

pointsInTemplateSet.M12 = [0.25 0.1 0.85]';
pointsInTemplateSet.M13 = [0.20 -0.1 0.7]';
pointsInTemplateSet.M14LKneeLateral = [0.2 0.0 0.5]';
pointsInTemplateSet.M15LKneeMedial = [0.05 0.0 0.5]';

rigidBodyMarkerNames.leftThigh = {'M10LTrochanter', 'M12', 'M13', 'M14LKneeLateral', 'M15LKneeMedial'};
% rigidBodyMarkerNames.leftThigh = {'M10LTrochanter', 'M12', 'M14LKneeLateral', 'M15LKneeMedial'};

pointsInTemplateSet.M16 = [-0.25 0.1 0.85]';
pointsInTemplateSet.M17 = [-0.20 -0.1 0.7]';
pointsInTemplateSet.M18RKneeLateral = [-0.2 0.0 0.5]';
pointsInTemplateSet.M19RKneeMedial = [-0.05 0.0 0.5]';

rigidBodyMarkerNames.rightThigh = {'M11RTrochanter', 'M16', 'M17', 'M18RKneeLateral', 'M19RKneeMedial'};

rigidBodyChain = MarkeredRigidBodyChain(pointsInTemplateSet, rigidBodyMarkerNames);

% generate the fake cloud data:
noise = 0.003;
names = rigidBodyChain.getMarkerNames();
for i = 1:length(names)
  pointsInCloudSet.([names{i} 'Cloud']) = pointsInTemplateSet.(names{i}) + randn(1) .* noise;
end

% test closest point finding
correspondedRigidBodyChain = MarkeredRigidBodyChain(pointsInCloudSet, []);
[closestPointIndexInCloud, closestPointIndexInTemplate] = ...
  rigidBodyChain.rigidBodySegments.pelvis.findClosestPointInCloudToAnyPointInBody( ...
  correspondedRigidBodyChain.getMarkersInMatrix());


%%
rigidBody = rigidBodyChain.rigidBodySegments.pelvis;
numberOfClosestPointsToFind = size(rigidBody.trackersInBodyFrame, 2);
closeCloudPoints = zeros(3, numberOfClosestPointsToFind);
closeTemplatePoints = zeros(3, numberOfClosestPointsToFind);
possibleCorrespondingPoints = correspondedRigidBodyChain.getMarkersInMatrix();
templatePoints = rigidBody.trackersInBodyFrame;
for i = 1:numberOfClosestPointsToFind
  [closestPointIndexInCloud, closestPointIndexInTemplate] = ...
    rigidBody.findClosestPointInCloudToAnyPointInBody( ...
    possibleCorrespondingPoints);
  
  closeCloudPoints(:, i) = possibleCorrespondingPoints(:, closestPointIndexInCloud);
  closeTemplatePoints(:, i) = templatePoints(:, closestPointIndexInTemplate);
  possibleCorrespondingPoints(:, closestPointIndexInCloud) = [];
  
  plotLine = @(a, b, color, lineWidth) line([a(1) b(1)], ...
    [a(2) b(2)], [a(3) b(3)], 'color', color, 'LineWidth', lineWidth);
  
  plotLine(closeTemplatePoints(:, i), closeCloudPoints(:, i), 'k', 2);
  hold on;
end

rigidBody.plotTheseWorldPoints(templatePoints, {'ok'});
rigidBody.plotTheseWorldPoints(correspondedRigidBodyChain.getMarkersInMatrix(), {'or'});
axis equal

% test correspondence finding
minSquaredErrorLimit = 0.1;

[correspondedRigidBodyChain, correspondences] = rigidBodyChain.findCorrespondencesOfRigidBodiesToPointCloud(pointsInCloudSet, ...
  'minSquaredErrorLimit', minSquaredErrorLimit)



%
for i = 1:length(correspondences)
  
correspondences{i}
end




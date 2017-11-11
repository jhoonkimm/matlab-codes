classdef MarkeredRigidBodyChain
  %MARKEREDRIGIDBODYCHAIN Set of points that move with rigid bodies
  % conneced together in a kinematic chain
  %   Useful for representing motion capture data of a human. Has methods
  %   for finding correpondences between different trials of the same
  %   chain to automatically name markers, filling in missing data using 
  %   rigid body assumptions, etc.
  
  properties
    
    rigidBodySegments = [];
    
    allMarkers = [];
    
    markersInRigidBodies = [];
    
  end
  
  methods
    
    function [this] = MarkeredRigidBodyChain(markers, markersInRigidBodies)
      this.allMarkers = markers;
      this.markersInRigidBodies = markersInRigidBodies;
      
      if (~isempty(markersInRigidBodies))
        rigidBodyNames = fieldnames(markersInRigidBodies);
        this.rigidBodySegments = [];
        for i = 1:length(rigidBodyNames)
          rigidBodyName = rigidBodyNames{i};
          theseRigidBodyMarkerNames = markersInRigidBodies.(rigidBodyName);
          this.rigidBodySegments.(rigidBodyName) = RigidBody;
          
          for markerNumber = 1:length(theseRigidBodyMarkerNames)
            markerName = theseRigidBodyMarkerNames{markerNumber};
            this.rigidBodySegments.(rigidBodyName) = ...
              this.rigidBodySegments.(rigidBodyName).addTrackerPositionInBodyFrame(markers.(markerName));
          end
        end
      end
    end
    
    %% some accessors:
    function [markerNames] = getMarkerNames(this)
      markerNames = fieldnames(this.allMarkers);
    end
    
    function [markerMatrix] = getMarkersInMatrix(this)
      names = this.getMarkerNames();
      markerMatrix = zeros(3, length(names));
      for i = 1:length(names)
        markerMatrix(:, i) = this.allMarkers.(names{i});
      end
    end

    function [markerMatrix] = getMarkersWithNames(this, theseNames)
      %       names = this.getMarkerNames();
      markerMatrix = zeros(3, length(theseNames));
      for i = 1:length(theseNames)
        markerMatrix(:, i) = this.allMarkers.(theseNames{i});
      end
    end
    
    function [index] = getIndexOfMarkerInBody(this, markerName, body)
      namesInBody = this.markersInRigidBodies.(body);
      for index = 1:length(namesInBody)
        if (strcmp(namesInBody{index}, markerName))
          return;
        end
      end
    end
    
    
    function [index] = getIndexOfMarkerInChain(this, markerName)
      names = this.getMarkerNames();
      for index = 1:length(names)
        if (strcmp(names{index}, markerName))
          return;
        end
      end
    end
    
    
    
    %% mutators:
    
    function [this] = removeMarkersWithNames(this, theseNames)
      for i = 1:length(theseNames)
        this.allMarkers = rmfield(this.allMarkers, theseNames{i});
      end
    end
    
    function [this] = removeAllMarkersExceptTheseNames(this, theseNames)
      names = this.getMarkerNames();
      for i = 1:length(names)
        if (isempty(cell2mat(strfind(theseNames, names{i}))))
          this.allMarkers = rmfield(this.allMarkers, names{i});
        end
      end
    end
    
    function [this] = addToAllMarkers(this, offset)
      names = this.getMarkerNames();
      for i = 1:length(names)
        this.allMarkers.(names{i}) = this.allMarkers.(names{i}) + offset';
      end
    end
    

    
    %% correspondence finding... let's us find which unnamed markers in a point cloud
    % correspond to which markers in this set of rigid bodies (say this
    % MarkeredRigidBodyChain represents a standing trial with known anatomical marker locations,
    % and the point cloud is from some later walking trial of the sameish setup, but
    % where we don't know which marker is at which anatomical location).
    function [correspondedRigidBodyChain, correspondences] = ...
        findCorrespondencesOfRigidBodiesToPointCloud(this, ...
        pointCloud, varargin)
      
      verbose = 0;
      correspondences = [];
      minSquaredErrorLimit = 0.1;
      additionalPointsToIncludeInCorrespondenceSearch = 1;
      
      for i = 1:2:length(varargin)
        if (strcmp('correspondences', varargin{i}))
          correspondences = varargin{i + 1};
        end
        if (strcmp('minSquaredErrorLimit', varargin{i}))
          minSquaredErrorLimit = varargin{i + 1};
        end
        if (strcmp('verbose', varargin{i}))
          verbose = varargin{i + 1};
        end
        if (strcmp('additionalPointsToIncludeInCorrespondenceSearch', varargin{i}))
          additionalPointsToIncludeInCorrespondenceSearch = varargin{i + 1};
        end
      end
      
      if (isa(pointCloud, 'MarkeredRigidBodyChain'))
        correspondedRigidBodyChain = pointCloud;
      else
        correspondedRigidBodyChain = MarkeredRigidBodyChain(pointCloud, []);
      end
      
      rigidBodyNames = fieldnames(this.rigidBodySegments);
      rigidBodyNumbersToCorrespond = 1:length(rigidBodyNames);
      
      % right, screw it, just order processing of bodies from high z to low z
      heights = zeros(1, length(rigidBodyNames));
      for i = 1:length(rigidBodyNames)
        segment = this.rigidBodySegments.(rigidBodyNames{i});
        markers = segment.trackersInBodyFrame;
        heights(i) = mean(markers(3,:));
      end
      
      while (~isempty(rigidBodyNumbersToCorrespond))
        heightsLeftToCorrespond = heights(rigidBodyNumbersToCorrespond);
        [maxHeight, maxHeightIndex] = max(heightsLeftToCorrespond);
        rigidBodyIndexToCorrespond = rigidBodyNumbersToCorrespond(maxHeightIndex);
        rigidBodyNumbersToCorrespond(maxHeightIndex) = [];
        
        rigidBodyName = rigidBodyNames{rigidBodyIndexToCorrespond};
        fprintf('finding correspondences for template body: %s\n', rigidBodyName);
        
        thisTemplate = this.rigidBodySegments.(rigidBodyName);
        thisCorrespondedRigidBodyChain = correspondedRigidBodyChain;
        
        
        theseConstrainedCorrespondences = [];
        theseConstrainedNotCorresponded = [];
        theseConstrainedCorrespondencesNamed = {};
        theseConstrainedNotCorrespondedNamed = {};
        
        theseMarkerNames = this.markersInRigidBodies.(rigidBodyName);
        cloudMarkerNames = thisCorrespondedRigidBodyChain.getMarkerNames();
        
        for correspondenceNumber = 1:length(correspondences)
          thisCorrespondence = correspondences{correspondenceNumber};
          indexIntoMarkerInTemplateBody = find(strcmp(theseMarkerNames, thisCorrespondence{1}));
          indexIntoMarkerInCloudBody = find(strcmp(cloudMarkerNames, thisCorrespondence{2}));
          if (~isempty(indexIntoMarkerInTemplateBody) && ~isempty(indexIntoMarkerInCloudBody))
            theseConstrainedCorrespondencesNamed{length(theseConstrainedCorrespondencesNamed) + 1} = ...
              thisCorrespondence;
          end
          if (isempty(indexIntoMarkerInTemplateBody) && ~isempty(indexIntoMarkerInCloudBody))
            % this marker isn't in this template, but we already know what
            % it is (so it is associated with a marker that is not in this
            % template, argal, is not associated with a marker in this
            % template.
            theseConstrainedNotCorrespondedNamed = [theseConstrainedNotCorrespondedNamed thisCorrespondence{2}];
          end
        end
        
        numberOfClosestPointsToFind = size(thisTemplate.trackersInBodyFrame, 2) + ...
          additionalPointsToIncludeInCorrespondenceSearch;
        
        % if a correspondence is already known for any of the points in
        % this template, we are going to shift the entire cloud to match
        % that point
        if (numel(theseConstrainedCorrespondencesNamed) > 0)
          templatePoint = this.getMarkersWithNames({theseConstrainedCorrespondencesNamed{1}{1}});
          cloudPoint = thisCorrespondedRigidBodyChain.getMarkersWithNames({theseConstrainedCorrespondencesNamed{1}{2}});
          shiftForAllPointsInCloud = templatePoint - cloudPoint;
          thisCorrespondedRigidBodyChain = thisCorrespondedRigidBodyChain.addToAllMarkers(shiftForAllPointsInCloud);
        end
        
        % now remove any points from the
        % possible cloud points that have been corresponded to, but not
        % from a marker that is in this template. We know that they
        % aren't part of this template.
        
        fprintf('theseConstrainedNotCorresponded:\n');
        theseConstrainedNotCorrespondedNamed
        thisCorrespondedRigidBodyChain = ...
          thisCorrespondedRigidBodyChain.removeMarkersWithNames(theseConstrainedNotCorrespondedNamed);
        
        
        % we need to keep doing this till we have ensured that we have
        % all of the constrained points in the set
        haveAllConstrainedPoints = 0;
        while (~haveAllConstrainedPoints)
          % use this to only use the closest n points to the current
          % template body
          [closestIndeces] = ...
            thisTemplate.getClosestNPointsToPointsInThisBody( ...
            thisCorrespondedRigidBodyChain.getMarkersInMatrix(), numberOfClosestPointsToFind, ...
            'shouldPlot', 0);
          
          closeMarkersToUse = thisCorrespondedRigidBodyChain.getMarkerNames();
          closeMarkersToUse = closeMarkersToUse(closestIndeces);
          
          % have to check and make sure constrained correspondences are in
          % the subset of close points
          haveAllConstrainedPoints = 1;
          for i = 1:length(theseConstrainedCorrespondencesNamed)
            thisCorrespondence = theseConstrainedCorrespondencesNamed{i};
            if (isempty(cell2mat(strfind(thisCorrespondence, thisCorrespondence{2}))))
              haveAllConstrainedPoints = 0;
              numberOfClosestPointsToFind = numberOfClosestPointsToFind + 1;
              fprintf('closest points don''t include constrained points, adding more points for correspondence search (%g)!\n', ...
                numberOfClosestPointsToFind);
              break;
            end
          end
        end
        
        thisCorrespondedRigidBodyChain = ...
          thisCorrespondedRigidBodyChain.removeAllMarkersExceptTheseNames(closeMarkersToUse);
        
        % convert from names of correspondences to indeces
        constrainedCorrespondencesIndeces = zeros(length(theseConstrainedCorrespondencesNamed), 2);
        for i = 1:length(theseConstrainedCorrespondencesNamed)
          constrainedCorrespondencesIndeces(i, 1) = this.getIndexOfMarkerInBody(theseConstrainedCorrespondencesNamed{i}{1}, rigidBodyName);
          constrainedCorrespondencesIndeces(i, 2) = thisCorrespondedRigidBodyChain.getIndexOfMarkerInChain(theseConstrainedCorrespondencesNamed{i}{2});
        end
        
        % find new correspondence indeces
        [theseCorrepondeces, correspondenceError] = ...
          thisTemplate.findCorrespondecesOfThisBodyToGivenPoints( ...
          thisCorrespondedRigidBodyChain.getMarkersInMatrix(), ...
          'verbose', 0, ...verbose, ...
          'constrainedCorrespondences', constrainedCorrespondencesIndeces, ...
          'minSquaredErrorLimit', minSquaredErrorLimit);
        
        % and convert new correspondence indeces back to names        
        correspondedMarkerNamesInTemplate = this.markersInRigidBodies.(rigidBodyName);
        correspondedMarkerNamesInTemplate = correspondedMarkerNamesInTemplate(theseCorrepondeces(:, 1));
        
        correspondedMarkerNamesInCloud = thisCorrespondedRigidBodyChain.getMarkerNames()';
        correspondedMarkerNamesInCloud = correspondedMarkerNamesInCloud(theseCorrepondeces(:, 2));
        
        latestCorrespondences = [correspondedMarkerNamesInTemplate; correspondedMarkerNamesInCloud];

        fprintf('latest correspondences are:\n%s', ...
          sprintf('%s\n', ...
          sprintf('%s  ', latestCorrespondences{1, :}), ...
          sprintf('%s  ', latestCorrespondences{2, :})));
        
        % add the new correspondences
        for correspondenceNumber = 1:length(correspondedMarkerNamesInTemplate)
          correspondences{length(correspondences) + 1} = ...
            {correspondedMarkerNamesInTemplate{correspondenceNumber}, ...
            correspondedMarkerNamesInCloud{correspondenceNumber}};
        end
        
        
        % check for conflicting marker correspondences (one marker
        % in the template set corresponded to two different markers in the
        % cloud set, and below we check that one marker in the cloud set is only corresponded
        % to one marker in the template set).
        uniqueRenamings = [];
        for i = 1:length(correspondences)
          first = correspondences{i}{1};
          second = correspondences{i}{2};
          if (isfield(uniqueRenamings, first))
            if (~strcmp(uniqueRenamings.(first), second))
              fprintf('current known renaming is %s -> %s, we just found a correspondence, %s -> %s\n', ...
                first, uniqueRenamings.(first), first, second);
              error('got conflicting marker correspondences!');
            end
          else
            uniqueRenamings.(first) = second;
          end
        end
        
        uniqeRenamingsFirsts = fieldnames(uniqueRenamings);
        correspondences = cell(length(uniqeRenamingsFirsts), 1);
        for i = 1:length(uniqeRenamingsFirsts)
          correspondences{i} = {uniqeRenamingsFirsts{i}, uniqueRenamings.(uniqeRenamingsFirsts{i})};
        end
        
        
        uniqueRenamings = [];
        for i = 1:length(correspondences)
          second = correspondences{i}{1};
          first = correspondences{i}{2};
          if (isfield(uniqueRenamings, first))
            if (~strcmp(uniqueRenamings.(first), second))
              error('got conflicting marker correspondences!');
            end
          else
            uniqueRenamings.(first) = second;
          end
        end
      end
    end
    
    
    
     %% correspondence finding... let's us find which unnamed markers in a point cloud 
    % correspond to which markers in this set of rigid bodies (say this
    % MarkeredRigidBodyChain represents a standing trial with known anatomical marker locations, 
    % and the point cloud is from some later walking trial of the sameish setup, but
    % where we don't know which marker is at which anatomical location).
    function [correspondedRigidBodyChain, correspondences] = ...
        findCorrespondencesOfRigidBodiesToPointCloudOld(this, ...
        pointCloud, varargin)
      
      verbose = 0;
      correspondences = [];
      minSquaredErrorLimit = 0.1;
      additionalPointsToIncludeInCorrespondenceSearch = 1;
      %       knownCorrespondences = [];
      
      for i = 1:2:length(varargin)
        if (strcmp('correspondences', varargin{i}))
          correspondences = varargin{i + 1};
        end
        if (strcmp('minSquaredErrorLimit', varargin{i}))
          minSquaredErrorLimit = varargin{i + 1};
        end
        if (strcmp('verbose', varargin{i}))
          verbose = varargin{i + 1};
        end
        if (strcmp('additionalPointsToIncludeInCorrespondenceSearch', varargin{i}))
          additionalPointsToIncludeInCorrespondenceSearch = varargin{i + 1};
        end
      end
      
      if (isa(pointCloud, 'MarkeredRigidBodyChain'))
        correspondedRigidBodyChain = pointCloud;
      else
        correspondedRigidBodyChain = MarkeredRigidBodyChain(pointCloud, []);
      end
      
      rigidBodyNames = fieldnames(this.rigidBodySegments);
      rigidBodyNumbersToCorrespond = 1:length(rigidBodyNames);
      
      %       numberOfCommonMarkersInChain = this.getNumberOfCommonMarkersInChain(rigidBodyNames);
      
      % right, screw it, just order processing of bodies from high z to low z
      heights = zeros(1, length(rigidBodyNames));
      for i = 1:length(rigidBodyNames)
        segment = this.rigidBodySegments.(rigidBodyNames{i});
        markers = segment.trackersInBodyFrame;
        heights(i) = mean(markers(3,:));
      end
      
      while (~isempty(rigidBodyNumbersToCorrespond))
        
        heightsLeftToCorrespond = heights(rigidBodyNumbersToCorrespond);
        [maxHeight, maxHeightIndex] = max(heightsLeftToCorrespond);
        rigidBodyIndexToCorrespond = rigidBodyNumbersToCorrespond(maxHeightIndex);
        rigidBodyNumbersToCorrespond(maxHeightIndex) = [];
        
        rigidBodyName = rigidBodyNames{rigidBodyIndexToCorrespond};
        fprintf('finding correspondences for template body: %s\n', rigidBodyName);
        
        %         if (strcmp('RPV', rigidBodyName))
        %           rigidBodyName
        %         end
        
        thisTemplate = this.rigidBodySegments.(rigidBodyName);
        
        theseConstrainedCorrespondences = [];
        theseConstrainedNotCorresponded = [];
        theseMarkerNames = this.markersInRigidBodies.(rigidBodyName);
        cloudMarkerNames = correspondedRigidBodyChain.getMarkerNames();
        
        for correspondenceNumber = 1:length(correspondences)
          thisCorrespondence = correspondences{correspondenceNumber};
          indexIntoMarkerInTemplateBody = find(strcmp(theseMarkerNames, thisCorrespondence{1}));
          indexIntoMarkerInCloudBody = find(strcmp(cloudMarkerNames, thisCorrespondence{2}));
          if (~isempty(indexIntoMarkerInTemplateBody) && ~isempty(indexIntoMarkerInCloudBody))
            thisCorrespondence
            theseConstrainedCorrespondences = [theseConstrainedCorrespondences; ...
              [indexIntoMarkerInTemplateBody indexIntoMarkerInCloudBody]];
          end
          if (isempty(indexIntoMarkerInTemplateBody) && ~isempty(indexIntoMarkerInCloudBody))
            % this marker isn't in this template, but we already know what
            % it is (so it is associated with a marker that is not in this
            % template, argal, is not associated with a marker in this
            % template.
            
            theseConstrainedNotCorresponded = [theseConstrainedNotCorresponded; ...
              indexIntoMarkerInCloudBody];
          end
        end
        % theseConstrainedCorrespondences
        
        numberOfClosestPointsToFind = size(thisTemplate.trackersInBodyFrame, 2) + ...
          additionalPointsToIncludeInCorrespondenceSearch;
        cloudPoints = correspondedRigidBodyChain.getMarkersInMatrix();
        
        % if a correspondence is already known for any of the points in
        % this template, we are going to shift the entire cloud to match
        % that point
        if (numel(theseConstrainedCorrespondences) > 0)
          indexIntoTemplate = theseConstrainedCorrespondences(1, 1);
          indexIntoCloud = theseConstrainedCorrespondences(1, 2);
          cloudPoint = cloudPoints(:, indexIntoCloud);
          templatePoint = thisTemplate.trackersInBodyFrame(:, indexIntoTemplate);
          cloudPoints = cloudPoints - repmat(cloudPoint - templatePoint, [1 size(cloudPoints, 2)]);
        end
        
        
        if (verbose > 9 && (numel(theseConstrainedCorrespondences) > 0) && size(theseConstrainedCorrespondences, 1) == 2)
          constrainedCorrespondenceFig = figure();
          subplot(2, 2, 2);
          thisTemplate.plotTheseWorldPoints(thisTemplate.trackersInBodyFrame, {'ok'});
          
          thisTemplate.plotTheseWorldPoints(cloudPoints(:, theseConstrainedCorrespondences(1, 2)), {'sr'})
          thisTemplate.plotTheseWorldPoints(cloudPoints(:, theseConstrainedCorrespondences(2, 2)), {'sg'})
          
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondences(1, 1)), {'*r'});
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondences(2, 1)), {'*g'});
          axis equal;
          
          title('before removing non correspondences');
        end
        
        % now remove any points from the
        % possible cloud points that have been corresponded to, but not
        % from a marker that is in this template. We know that they
        % aren't part of this template.
        allCloudIndeces = 1:size(cloudPoints, 2);
        indecesFromCloudWeCanConsider = 1:size(cloudPoints, 2);
        indecesFromCloudWeCanConsider(theseConstrainedNotCorresponded) = [];
        cloudPoints = cloudPoints(:, indecesFromCloudWeCanConsider);
        
        fprintf('theseConstrainedNotCorresponded:\n');
        theseConstrainedNotCorresponded'
        
        if (numel(theseConstrainedCorrespondences) > 0)
          constrainedIndecesInCloudPoints = zeros(size(allCloudIndeces));
          [sorted, ordering] = sort(theseConstrainedCorrespondences(:, 2));
          constrainedIndecesInCloudPoints(theseConstrainedCorrespondences(:, 2)) = 1;
          %           constrainedIndecesInCloudPoints(theseConstrainedCorrespondences(:, 2)) = theseConstrainedCorrespondences(:, 2);
          
          constrainedIndecesInCloudPoints(theseConstrainedNotCorresponded) = [];
          %           goodIndeces = find(constrainedIndecesInCloudPoints);
          %           theseConstrainedCorrespondences(:, 2) = constrainedIndecesInCloudPoints(goodIndeces);
          
          newMarkerNamesSorted = find(constrainedIndecesInCloudPoints);
          theseConstrainedCorrespondences(:, 2) = newMarkerNamesSorted(ordering);
        end
        
        if (verbose > 9 && (numel(theseConstrainedCorrespondences) > 0) && size(theseConstrainedCorrespondences, 1) == 2)
          figure(constrainedCorrespondenceFig);
          subplot(2, 2, 3);
          thisTemplate.plotTheseWorldPoints(thisTemplate.trackersInBodyFrame, {'ok'});
          
          thisTemplate.plotTheseWorldPoints(cloudPoints(:, theseConstrainedCorrespondences(1, 2)), {'sr'})
          thisTemplate.plotTheseWorldPoints(cloudPoints(:, theseConstrainedCorrespondences(2, 2)), {'sg'})
          
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondences(1, 1)), {'*r'});
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondences(2, 1)), {'*g'});
          axis equal;
          
          title('after removing non correspondences');
        end
        
        % cloudPointsWeCanConsider = cloudPoints(:, indecesFromCloudWeCanConsider);
        
        haveAllConstrainedPoints = 0;
        while (~haveAllConstrainedPoints)
          % we need to keep doing this till we have ensured that we have
          % all of the constrained points in the set
          
          
          % use this to only use the closest n points to the current
          % template body
          [closestIndeces] = ...
            thisTemplate.getClosestNPointsToPointsInThisBody( ...
            cloudPoints, numberOfClosestPointsToFind, ...
            'shouldPlot', 0);
          %           drawnow;
          %           pause(0.01);
          
          cloudPointsToUse = cloudPoints(:, closestIndeces);
          
          % have to switch constrained correspondences to be for the new
          % subset of cloud points
          theseConstrainedCorrespondencesIndexedIntoSubset = ...
            theseConstrainedCorrespondences;
          haveAllConstrainedPoints = 1;
          for i = 1:size(theseConstrainedCorrespondences, 1)
            if (isempty(find(closestIndeces == theseConstrainedCorrespondences(i, 2))))
              % this means that the points that were passed in to use were
              % not enough to include a constrained point!
              haveAllConstrainedPoints = 0;
              numberOfClosestPointsToFind = numberOfClosestPointsToFind + 1;
              fprintf('closest points don''t include constrained points, adding more points for correspondence search (%g)!\n', ...
                numberOfClosestPointsToFind);
              break;
            end
            theseConstrainedCorrespondencesIndexedIntoSubset(i, 2) = ...
              find(closestIndeces == theseConstrainedCorrespondences(i, 2));
          end
        end
        
        if (verbose > 9 && (numel(theseConstrainedCorrespondences) > 0) && size(theseConstrainedCorrespondences, 1) == 2)
          figure(constrainedCorrespondenceFig);
          subplot(2, 2, 4);
          thisTemplate.plotTheseWorldPoints(thisTemplate.trackersInBodyFrame, {'ok'});
          
          thisTemplate.plotTheseWorldPoints(cloudPointsToUse(:, theseConstrainedCorrespondencesIndexedIntoSubset(1, 2)), {'sr'})
          thisTemplate.plotTheseWorldPoints(cloudPointsToUse(:, theseConstrainedCorrespondencesIndexedIntoSubset(2, 2)), {'sg'})
          
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondencesIndexedIntoSubset(1, 1)), {'*r'});
          thisTemplate.plotTheseWorldPoints( ...
            thisTemplate.trackersInBodyFrame(:, theseConstrainedCorrespondencesIndexedIntoSubset(2, 1)), {'*g'});
          axis equal;
          title('after finding closest point set');
        end
        %%
        
        [theseCorrepondeces, correspondenceError] = ...
          thisTemplate.findCorrespondecesOfThisBodyToGivenPoints( ...
          cloudPointsToUse, ...
          'verbose', 0, ...verbose, ...
          'constrainedCorrespondences', theseConstrainedCorrespondencesIndexedIntoSubset, ...
          'minSquaredErrorLimit', minSquaredErrorLimit);
        
        % index conversion from
        % cloud indeces we did consider to 
        % cloud indeces we could consider to 
        % all cloud indeces.
        theseCorrepondeces(:, 2) = closestIndeces(theseCorrepondeces(:, 2));
        corresponds = theseCorrepondeces(:, 2);
        theseCorrepondeces(:, 2) = indecesFromCloudWeCanConsider(corresponds);

        % convert theseCorrespondences (indeces) to correspondences (named markers)
        
        latestCorrespondences = cell(2, size(theseCorrepondeces, 1));
        for correspondenceNumber = 1:size(theseCorrepondeces, 1)
          thisCorrespondence = theseCorrepondeces(correspondenceNumber, :);
          thisTemplateMarkerName = theseMarkerNames{thisCorrespondence(1)};
          thisCloudMarkerName = cloudMarkerNames{thisCorrespondence(2)};
          
          %           fprintf('new correspondence: ');
          %           {thisTemplateMarkerName, thisCloudMarkerName}
          latestCorrespondences{1, correspondenceNumber} = thisTemplateMarkerName;
          latestCorrespondences{2, correspondenceNumber} = thisCloudMarkerName;
          
          correspondences{length(correspondences) + 1} = ...
            {thisTemplateMarkerName, thisCloudMarkerName};
        end
        fprintf('latest correspondences are:\n%s', ...
          sprintf('%s\n', ...
          sprintf('%s  ', latestCorrespondences{1, :}), ...
          sprintf('%s  ', latestCorrespondences{2, :})));
        
        % check for conflicting marker correspondences (one marker
        % in the template set corresponded to two different markers in the
        % cloud set, and below we check that one marker in the cloud set is only corresponded
        % to one marker in the template set).
        uniqueRenamings = [];
        for i = 1:length(correspondences)
          first = correspondences{i}{1};
          second = correspondences{i}{2};
          if (isfield(uniqueRenamings, first))
            if (~strcmp(uniqueRenamings.(first), second))
              fprintf('current known renaming is %s -> %s, we just found a correspondence, %s -> %s\n', ...
                first, uniqueRenamings.(first), first, second);
              error('got conflicting marker correspondences!');
            end
          else
            uniqueRenamings.(first) = second;
          end
        end
        
        uniqeRenamingsFirsts = fieldnames(uniqueRenamings);
        correspondences = cell(length(uniqeRenamingsFirsts), 1);
        for i = 1:length(uniqeRenamingsFirsts)
          correspondences{i} = {uniqeRenamingsFirsts{i}, uniqueRenamings.(uniqeRenamingsFirsts{i})};
        end
        
        
        uniqueRenamings = [];
        for i = 1:length(correspondences)
          second = correspondences{i}{1};
          first = correspondences{i}{2};
          if (isfield(uniqueRenamings, first))
            if (~strcmp(uniqueRenamings.(first), second))
              error('got conflicting marker correspondences!');
            end
          else
            uniqueRenamings.(first) = second;
          end
        end
        
        
      end
    end
    
    %%
    function [numbersOfCommonMarkersInChain] = getNumberOfCommonMarkersInChain(this, rigidBodyNames)
      numbersOfCommonMarkersInChain = zeros(1, length(rigidBodyNames));
      for rigidBodyNumber = 1:length(rigidBodyNames)
        markersInThisBody = this.markersInRigidBodies.(rigidBodyNames{rigidBodyNumber});
        numbersOfCommonMarkersInChain(rigidBodyNumber) = 0;
        
        for markerNumber = 1:length(markersInThisBody)
          numbersOfCommonMarkersInChain(rigidBodyNumber) = ...
            numbersOfCommonMarkersInChain(rigidBodyNumber) + ...
            this.getNumberOfSegmentsWithThisMarker(markersInThisBody{markerNumber}) - 1;
        end
      end
    end
    
    %%
    function [numberOfSegmentsWithThisMarker] = ...
        getNumberOfSegmentsWithThisMarker(this, markerName)
      rigidBodyNames = fieldnames(this.markersInRigidBodies);
      numberOfSegmentsWithThisMarker = 0;
      for rigidBodyNumber = 1:length(rigidBodyNames)
        numberOfSegmentsWithThisMarker = numberOfSegmentsWithThisMarker + ...
          sum(strcmp(markerName, this.markersInRigidBodies.(rigidBodyNames{rigidBodyNumber})));
      end
    end
    
    
  end
  
end


function [] = visualizeMobileCartMotion(mobileCartMotion, varargin)

recordName = '';
orientationPlotInterval = 480/10;
arrowLength = []; %100;
% axisHandlesToPlotOn = [];
for i = 1:2:length(varargin)
  if (strcmp('recordName', varargin{i}))
    recordName = varargin{i + 1};
  end
  if (strcmp('orientationPlotInterval', varargin{i}))
    orientationPlotInterval = varargin{i + 1};
  end
  if (strcmp('arrowLength', varargin{i}))
    arrowLength = varargin{i + 1};
  end
end


%%
plot(mobileCartMotion.positions(:,1), mobileCartMotion.positions(:,2), 'k');
hold on;
axis equal
if (isempty(arrowLength))
  motionBounds = max(max(mobileCartMotion.positions) - min(mobileCartMotion.positions));
  arrowLength = motionBounds / 20;
end
for i = 1:orientationPlotInterval:length(mobileCartMotion.angles)
  position = mobileCartMotion.positions(i, :);
  angle = mobileCartMotion.angles(i);
  endPoint = position' + arrowLength * [cos(angle) sin(angle)]';
  
  %   [position(1) endPoint(1)]
  %   [position(2) endPoint(2)]
  line([position(1) endPoint(1)], [position(2) endPoint(2)], 'color', 'b', 'LineWidth', 4);
  plot(position(1), position(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
  
  %   arrow(position,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir);
  %   arrow(position, position + endPoint', arrowLength, 'FaceColor', 'none', 'EdgeColor', 'b'); %'MarkerFaceColor', 'off');
end
axis equal
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('mobile cart motion %s', recordName));




end


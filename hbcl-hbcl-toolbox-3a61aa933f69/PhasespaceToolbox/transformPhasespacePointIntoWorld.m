function [pointInWorld] = transformPhasespacePointIntoWorld(point, ...
  transformationFromPhasespaceToCart, phasespaceScale, cartPose)

pointInCartFrame = transformPoint2D(point, transformationFromPhasespaceToCart);
cartPose(1:2) = cartPose(1:2) * phasespaceScale;
pointInWorld = transformPoint2D(pointInCartFrame, cartPose);

end
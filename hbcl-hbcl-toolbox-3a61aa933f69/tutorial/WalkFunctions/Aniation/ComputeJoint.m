function [Hip Lknee Rknee Lankle Rankle]=ComputeJoint(ARhip, ALhip, ARknee, ALknee )
%% This function calculates the joints position in sagittal plan. It assume
%% hip joint is at the origin. ARhip and ALhip are hip joint angle with zero
%% average and extension is in positive direction. ARknee and ALknee are
%% knee joint angles and extension is in positive direction. Knee is full
%% extended with zero knee angle.
%% Writed by Paul 2010/11/4

Hip=[0 0]';
Lknee=Hip+[-sin(ALhip) -cos(ALhip)]';
Rknee=Hip+[-sin(ARhip) -cos(ARhip)]';
Lankle=Lknee+[-sin(ALhip-ALknee) -cos(ALhip-ALknee)];
Rankle=Rknee+[-sin(ARhip-ARknee) -cos(ARhip-ARknee)];

end
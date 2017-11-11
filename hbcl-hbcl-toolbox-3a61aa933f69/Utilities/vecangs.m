function [ang3d, xang, yang, zang] = vecangs(v1,v2)

% VECANGS find angle between two vectors
%  Computes angle between two vecotrs in 
%  3d and in all 2d planes

% Peter Adamczyk

warning off MATLAB:divideByZero

ang3d = acos( dot(v1,v2)/(vecnorm(v1)*vecnorm(v2)) ) ; % acos( (v1*v2')/(norm(v1)*norm(v2)) ) ;
%             if isnan(ang3d);
%                 [v1;v2]
%                 (v1*v2')/(norm(v1)*norm(v2))
%             end
%             up3d = atan2(v2(3),norm(v2(1:2))) - atan2(v1(3),norm(v1(1:2))) ;
xang = acos( dot(v1([2,3]),v2([2,3]))/(vecnorm(v1([2,3]))*vecnorm(v2([2,3]))) ) ;   % atan2(v2(3),v2(2)) - atan2(v1(3),v1(2)) ;
yang = acos( dot(v1([1,3]),v2([1,3]))/(vecnorm(v1([1,3]))*vecnorm(v2([1,3]))) ) ;   % atan2(v2(3),v2(1)) - atan2(v1(3),v1(1)) ;
zang = acos( dot(v1([1,2]),v2([1,2]))/(vecnorm(v1([1,2]))*vecnorm(v2([1,2]))) ) ;   % atan2(v2(1),v2(2)) - atan2(v1(1),v1(2)) ;
end

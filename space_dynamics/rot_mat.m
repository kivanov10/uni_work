function [rot_mat]=rot_mat(theta,rot_mode)
% DCM for sequential rotiations selection function. It accepts the angle
% and the axis around which it happens. The notation is kept to the
% standard 1-2-3. In other words for a rotation of 3-1-3 the function would
% be called three times with the corresponding angles and the numbers
% 3,1,3.
% Input: theta - angle by which the vector/matrix will be rotated
%        rot_mode - around which axis will the rotation be done. It uses
%                   the numerical representation 1,2,3 rather than using 
%                   x,y,z.
% Output: rot_mat- rotation matrix for selected axis and angle               

% simple switch case for ease of use in larger scripts
switch rot_mode
    case 1
        %around the x-axis
        rot_mat = [1     0          0
                   0 cos(theta) sin(theta)
                   0 -sin(theta) cos(theta)];
    case 2
        %around the y-axis
        rot_mat = [cos(theta)  0  -sin(theta)
                   0           1       0
                   sin(theta)  0  cos(theta)];
    case 3
        %around the z-axis
        rot_mat = [cos(theta)     sin(theta)  0
                   -sin(theta)    cos(theta)  0
                   0                  0       1];
 end
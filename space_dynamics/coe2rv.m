function [r,v] = coe2rv(coe, mu)
% Converts orbital elements to cartesian coordinates
% for position and velocity
% Accepts an array for orbital elements as seen below
%
% Inputs: coe = [a, e, i, omega, w, theta]'
%         mu  = gravitational parameter of the Earth
% Outputs: r - position vector
%          v - velocity vector

% separate elements into variables for ease of use and readibility
a = coe(1);
e = coe(2);
i = coe(3);
omega = coe(4);
w = coe(5);
theta = coe(6);

% r vector defininition and specific angular momentum
r_mag = (a*(1-e^2))/(1+e*cos(theta));
h = sqrt(mu*a*(1-e^2));
r = [r_mag*cos(theta),r_mag*sin(theta),0]';

% Rotation matrix creation
% Note: Line 22 uses a separate function where rotations around the X,Y,Z
% axes is done via a switch case
coe_rv_mat = (rot_mat(w,3)*rot_mat(i,1)*rot_mat(omega,3))';

r = coe_rv_mat*r; %New reference frame R after rotation

v = coe_rv_mat*[-mu/h*sin(theta), mu/h*(cos(theta)+e), 0]'; %New v vector
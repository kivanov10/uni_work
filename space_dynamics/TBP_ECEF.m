function [dXdt] = TBP_ECEF(t,X,mu)
% ode45 function that accepts ECI coordinates and velocity data and gives
% back the ECEF equivalent. Uses the transport theorem to go from one
% reference frame to the next.
% Inputs: t - time period of integration [s] Example: [0, 1000]
%         X - initial state vector guess (1x6 shape): Example: [r_i; v_i]
%         mu - gravitational parameter of the planet
% Output: dXdt - the derivates of the position and velocity in the ECEF
%                frame. The configuration should be the same as the intial
%                guess state vector [vel_j; accel_j]

%ECI frame data
r_i = X(1:3,:); 
v_i = X(4:6,:);

r_norm_i = norm(r_i);
a_i = -mu/r_norm_i^3 * r_i;

earth_rot = (2*pi)/(24*60*60); %rads/s
omega_j = [0;0;earth_rot];

%%%
% rotate the r vector to ECEF by doing a euler rotation around the 3/z axis
r_j = eye(3)*r_i;

%%% velocity
v_j  = v_i - cross(omega_j,r_j);
% a_i_j = a_i - cross(omega_j,v_i);


%%% acceleration
% euler = cross(a_i,r_i);
% coriolis = 2*cross(omega_j,v_i);
centrif = cross(omega_j,cross(omega_j,r_j));

coriolis = 2*cross(omega_j,v_j);

a_j = a_i - coriolis - centrif;

dXdt = [v_j;a_j];

%vdot_j = vdot_i + vdot_i x r_i + 2*omega_j x v_j - omega_j x omega_j x r_i
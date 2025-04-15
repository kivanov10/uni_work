function [dXdt] = TBP_ECI(t, X, mu)
% This function is made with the intent of being supplied to ode45. 
% The inputs and the calculations below represent the equations of motion
% of the satellite.
% Inputs: t - time period for the differential solver
%         X - state vector of the satellite. Rows 1:3 are position; Rows
%         4:6 are velocity
%         mu - gravitational parameter of the Earth
% Outputs: dXdt - this is the first order derivative of the state vector.
%                 When passed into the ode, the output of the
%                 differentiation will be position and velocity data


r = X(1:3,:);
v = X(4:6,:);
r_norm = norm(r);
% acceleration calculation
v_dot = -mu/r_norm^3 * r;

dXdt = [v;v_dot];

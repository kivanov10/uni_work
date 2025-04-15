function [dXdt] = AttittudeDynamics(t, X, I)
% ode 45 function. Calculates the torque-free motion of a rigid body. This
% is done specifically for a 3-1-3 rotation in Euler rotations.
% Inputs: t - time period of integration
%         X - state vector of rotational elements of the body. The
%         structure takes the form [theta_1; theta_2; theta_3; omega_x;
%         omega_y; omega_z], where the theta angles are the Euler rotations
%         and the omegas are the angular velocities in the three axes.
%         I - inertia matrix of the body. This assumes no products of
%         inertia, i.e principle axes.
% Outputs: dXdt - derivative values of the inputs. This is done using the
%                 Euler equations for change in theta angles and angular
%                 velocity for the 3-1-3 rotation.

% Extracting the diagonal inertia terms
    I_1 = I(1,1);
    I_2 = I(2,2);
    I_3 = I(3,3);
% Extracting the Euler angles
    theta_1 = X(1);
    theta_2 = X(2);
    theta_3 = X(3);
% Extracting the angular velocities
    omega_x = X(4);
    omega_y = X(5);
    omega_z = X(6);

% Euler equations for torque-free motion (L=0)
    w_x_dot = ((I_2-I_3)*omega_y*omega_z)/I_1;
    w_y_dot = ((I_3-I_1)*omega_z*omega_x)/I_2;
    w_z_dot = ((I_2-I_1)*omega_x*omega_y)/I_3;
% B matrix for change in angular vel
    B_theta_mat = (1/sin(theta_2)) * [sin(theta_3), cos(theta_3), 0;
                                     sin(theta_2)*cos(theta_3), -sin(theta_2)*sin(theta_3),     0;
                                     -cos(theta_2)*sin(theta_3), -cos(theta_2)*cos(theta_3),  sin(theta_2)]; 
% Change in angular vel
    theta_dot = B_theta_mat*[omega_x;omega_y;omega_z];

% Output vector
    dXdt = [theta_dot(1);theta_dot(2);theta_dot(3);w_x_dot;w_y_dot;w_z_dot];
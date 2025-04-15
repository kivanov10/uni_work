clc;
close all;
clearvars;
clear global;
set(0, 'DefaultTextInterpreter', 'latex');

% This script covers Part 2 of the Coursework (Week 6-8)
I_xx = 0.07583; % [kg m^2]
I_yy = 0.05833; % [kg m^2]
I_zz = 0.02916; % [kg m^2]

inert_mat = [I_xx  0     0;
              0   I_yy   0;
              0    0   I_zz];

r_0 = [5300.64;
       17575.73;
       -138.50]; % [km]

v_0 = [-4.2880;
       -1.9373;
       -0.6026]; % [km]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_0 = cross(r_0,v_0);

T = v_0/norm(v_0);
W = h_0/norm(h_0);
N = cross(T,W);
NTW_DCM = [N, T, W]; %[IT]
NTW_DCM_inv = inv(NTW_DCM);
NTW_DCM_trans = NTW_DCM';

disp('NTW_DCM');
disp(NTW_DCM);
disp('NTW_DCM invervse');
disp(NTW_DCM_inv);
disp('NTW_DCM transpose');
disp(NTW_DCM_trans);

BT = [0.7146, 0.6131, -0.3368;
      -0.6337, 0.7713, 0.0594;
      0.2962, 0.1710, 0.9397];

BI = BT*NTW_DCM'; %[BT]*[TI]=[BI]
disp('BI matrix:');
disp(BI);


%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Princple angle
phi = acos((trace(BI)-1)/2);
% Principle axes
e_1 = (BI(2,3)-BI(3,2))/(2*sin(phi));
e_2 = (BI(3,1)-BI(1,3))/(2*sin(phi));
e_3 = (BI(1,2)-BI(2,1))/(2*sin(phi));

e_test = sqrt(e_1^2+e_2^2+e_3^2); % = 1 (good result)

% Quaternion conversion

q_0 = 0.5*sqrt(trace(BI)+1);
q_1 = (BI(2,3)-BI(3,2))/(4*q_0);
q_2 = (BI(3,1)-BI(1,3))/(4*q_0);
q_3 = (BI(1,2)-BI(2,1))/(4*q_0);

quat_test = sqrt(q_0^2 + q_1^2 + q_2^2 + q_3^2); % = 1 (good result)


% Getting the Euler angles 

theta_1 = atan2(BI(3,1),BI(3,2));
theta_2 = acos(BI(3,3));
theta_3 = atan2(BI(1,3),BI(2,3));


%%%%%%%%%%%%%%%%%%%%%%%%%%% WEEK 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega_BI = [-8.0862e-5;1.4258e-5;2.2559e-4];

X_omega = [theta_1, theta_2, theta_3, omega_BI(1), omega_BI(2), omega_BI(3)];

t_period = linspace(0,3600,360);
options = odeset(RelTol = 1e-10);
[t_period,X_ode_omega] = ode45(@(t_period,X_ode_omega) ...
          AttittudeDynamics(t_period,X_ode_omega, inert_mat) , ...
          t_period,X_omega,options);

figure()
hold on;
grid on;
plot(t_period, X_ode_omega(:,1), 'b');
plot(t_period, X_ode_omega(:,2), 'r');
plot(t_period, X_ode_omega(:,3), 'g');
yline(0, '--');
xlabel('Time [s]');
ylabel('Euler Angles [rads]');
title('Euler angle evolution');
legend('\theta_x', '\theta_y', '\theta_z');

figure()
hold on;
grid on;
plot(t_period, X_ode_omega(:,4), 'b');
plot(t_period, X_ode_omega(:,5), 'r');
plot(t_period, X_ode_omega(:,6), 'g');
xlabel('Time [s]');
ylabel('Angular velocity [rads/s]');
title('Angular velocity evolution');
legend('\omega_x', '\omega_y', '\omega_z');

% Angular momentum and rotational kinetic energy
H_store = zeros(1,length(t_period));
T_store = zeros(1,length(t_period));

for jj=1:length(t_period)
    H_x = I_xx*X_ode_omega(jj,4);
    H_y = I_yy*X_ode_omega(jj,5);
    H_z = I_zz*X_ode_omega(jj,6);
    H_store(jj) = norm([H_x;H_y;H_z]);
    
    T_1 = H_x*X_ode_omega(jj,4);
    T_2 = H_y*X_ode_omega(jj,5);
    T_3 = H_z*X_ode_omega(jj,6);
    T_store(jj) = 0.5*(T_1+T_2+T_3);
end

figure();
hold on;
grid on;
yyaxis left;
plot(t_period, H_store);
xlabel('Time [s]');
ylabel('Angular momentum [$kg^2~m/s$]');
title('Conservation of angular momentum and rotational kinetic energy')
ylim([-1e-5,1.5e-5]);
yyaxis right;
plot(t_period, T_store);
ylabel('Rotational Kinetic Energy [J]');
ylim([-1e-8,1e-8]); 
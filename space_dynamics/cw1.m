clc;
close all;
clearvars;
clear global;
set(0, 'DefaultTextInterpreter', 'latex');

% This script covers Part 1 of the coursework (Week 2-5)
global kepler_iter

%%%%%%%%%%%%%%% CW data %%%%%%%%%%%%%%%%%
mu = 398600.4418; % [km^3/s^2]

R_e = 6378.137; % [km]  
w_e = 7.2921e-5; % [rad/s]


%%%%%%%%%%%%%%%%%%%%%% Part 1 data %%%%%%%%%%%%%%%%%%%%
a = 19052.49; % [km]
e = 0.6516; % eccentricity

i = 10.02; % deg
i = deg2rad(i); % rad

omega = 250.77; % deg
omega = deg2rad(omega); % rad

w = 310.67; % deg
w = deg2rad(w); % rad

M_0 = 8.77; % deg
M_0=deg2rad(M_0); % rad

tol_Kepler = 10e-10;

%%%%%%%%%%%%%%%%%% Week 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
[E_final, theta_final] = Kepler(e,M_0,tol_Kepler);

% Inbuilt MATLAB fzero function testing for Kepler equation
% [E_fzero, theta_fzero] = Kepler_fzero(e,M_0,tol_Kepler);

E_final = mod(E_final, 2*pi);
% print_kepler = ['Kepler iterations:', kepler_iter];


% Old test values
disp('Kepler E and theta values (in deg):');
disp(rad2deg(E_final));
disp(rad2deg(theta_final)); %49.5deg => periapsis
disp('Kepler iterations needed:');
disp(kepler_iter);
% 
% disp('Fzero Kepler');
% disp(rad2deg(E_fzero));
% disp(rad2deg(theta_fzero));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Week 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coe = [a, e, i, omega, w, theta_final]';

period = 2*(pi/sqrt(mu))*a^(3/2);

orbits = period*10; % can easily expand the number of rotations around Earth

t_increments = 1000;

% uses orbits to make it easier to increase time
t_period = linspace(0,orbits,t_increments); 

n = sqrt(mu/a^3);

M_evolve = M_0+n*t_period;

E_evolve = zeros(1,t_increments);
theta_evolve = zeros(1,t_increments);
kepler_iter_store = zeros(1,t_increments);
%using the M_evolve values in Kepler equations to get E and theta for all t
%points
for k=1:length(t_period)
    [E_evolve(k),theta_evolve(k)] = Kepler(e,M_evolve(k),tol_Kepler);
    kepler_iter_store(k) = kepler_iter;
end

disp("Average iterations need for Kepler fzero:");
disp(mean(kepler_iter_store));

%%%%%%%% Modulus %%%%%%%%
%Making the data easier to read in a plot
E_evolve = mod(E_evolve,2*pi);
theta_evolve = mod(theta_evolve,2*pi);
M_evolve = mod(M_evolve,2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Using tiledlayout to save on space and windows opening, easier to read on
%separate graphs
tiles = tiledlayout(1,3);
tiles.TileSpacing = 'compact';
nexttile;
plot(t_period, E_evolve, "b", "LineWidth",2);
xlabel(tiles,"Time (s)");
ylabel("Eccentric anomaly (rad)");
title("Eccentric anomaly through the orbit");
grid on;
hold off;
nexttile;


plot(t_period, theta_evolve, "b", "LineWidth",2);
ylabel("True anomaly (rad)");
title("True anomaly through the orbit");
grid on;
hold off;
nexttile;

plot(t_period, M_evolve, "b", "LineWidth",2);
ylabel("Mean anomaly (rad)");
title("Mean anomaly through the orbit");
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% COE to RV %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Storage for r and v

r_evolve = zeros(3,length(period));
v_evolve = zeros(3,length(period));

for j=1:length(t_period)

    coe = [a, e, i, omega, w, theta_evolve(j)]';
    [r_evolve(:,j),v_evolve(:,j)] = coe2rv(coe,mu);

end

X = [r_evolve; %state vector creation; 6x1000 size
    v_evolve];
%%%%%%%%%%%%%%%%%%%% Plotting ECI coord %%%%%%%%%%%%%%%%%%%%%%%
disp('ECI first point in orbit state vector:');
disp(X(:,1));
disp('ECI last point in orbit state vector:');
disp(X(:,end));

figure % new figure, stops tiledlayout

% r vector coords plotting
plot3(X(1,:),X(2,:),X(3,:), 'DisplayName', 'ECI Orbit');
hold on;
grid on;
% first and final r vector points
plot3(X(1,1),X(2,1),X(3,1),'ok','MarkerFaceColor','y', 'DisplayName', ...
    'Start Point');
plot3(X(1,end),X(2,end),X(3,end),'pentagram','MarkerFaceColor','r', ...
    'DisplayName', 'End Point');
title('ECI 3D orbit coordinates visualisation');
L1 = legend;
L1.AutoUpdate = 'off'; 
make_earth;

% Orbit unit vectors
r_unit = X(1:3, 1);
v_unit = X(4:6, 1);

h_unit = cross(r_unit,v_unit);
e_unit = cross(v_unit,h_unit)/mu - r_unit/norm(r_unit);

% e,h,p unit vector triad calc
ie = e_unit/norm(e_unit);
ih = h_unit/norm(h_unit);
ip = cross(ih, ie)/norm(cross(ih, ie));

%plotting unit vector arrows for orbital plane
quiver3(0,0,0,ie(1),ie(2),ie(3),1e4, 'k');
text(ie(1)*1e4,ie(2)*1e4,ie(3)*1e4,'$\hat{e}$');
quiver3(0,0,0,ih(1),ih(2),ih(3),1e4, 'k');
text(ih(1)*1e4,ih(2)*1e4,ih(3)*1e4,'$\hat{h}$');
quiver3(0,0,0,ip(1),ip(2),ip(3),1e4, 'k');
text(ip(1)*1e4,ip(2)*1e4,ip(3)*1e4,'$\hat{p}$');


%%%%%%%%%%%%%%%%%%%%% WEEK 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ode
options = odeset(RelTol = 1e-10);
[t,X_ode] = ode45(@(t,X_ode) TBP_ECI(t,X_ode,mu) ,t_period,[X(:,1)],options);
X_ode = X_ode'; %from 1000x6 to 6x1000 to match original X vector

figure % this figure is to plot the orbits and compare to the coe one
% r vector coords plotting
plot3(X_ode(1,:),X_ode(2,:),X_ode(3,:), 'DisplayName', ...
    'Integrated Orbit Coordinates', 'LineWidth', 2); %integrated r values
hold on;
grid on;
plot3(X(1,:),X(2,:),X(3,:), '--', 'DisplayName', 'Original OE orbit', ...
    'LineWidth', 2); %original coe values
title('Comparison between the ode45 and coe orbits');
%Code below generates a legend and stops it from adding the axis vectors
%and earth mesh to the legend
L1 = legend;
L1.AutoUpdate = 'off'; 

% first and final r vector points (uncomment to add)
% plot3(X(1,1),X(2,1),X(3,1),'ok','MarkerFaceColor','b');
% plot3(X(1,end),X(2,end),X(3,end),'ok','MarkerFaceColor','r');

%Generate an Earth with ECI axis unit vectors
make_earth;


%%%% Error propagation between integrated and coe r and v
r_error = X(1:3,:) - X_ode(1:3,:); 
r_error_store = zeros(1,t_increments);

v_error = X(4:6,:) - X_ode(4:6,:);
v_error_store = zeros(1,t_increments);

for k=1:length(t_period)
    r_error_store(k) = norm(r_error(1:3,k));
    v_error_store(k) = norm(v_error(1:3,k));
end

%%%%%%%%%%%%%%%%% Plot for error comparison between ode45 and coe %%%%%%%%
figure();
yyaxis left;
plot(t_period, r_error_store);
xlabel('Time [s]');
ylabel('Radius error [m]');
title(['Error evolution between integrated and orbital elements extracted' ...
       'values for radius and velocity'])

yyaxis right;
plot(t_period, v_error_store);
ylabel('Velocity error [m/s]');

% adding more orbits makes the orbits to diverge more and more (orbit
% determination problem, needs batch or sequential processors to fix)

h_ode_norm = zeros(1,t_increments);

for k=1:length(t_period)
    h_ode = cross(X_ode(1:3,k),X_ode(4:6,k));
    h_ode_norm(k) = norm(h_ode);

end

figure()
static_h = sqrt(mu*a*(1-e^2));
plot(t_period,h_ode_norm,'r','LineWidth', 2, 'DisplayName', 'Norm Value');
yline(static_h, 'b--','LineWidth', 2, 'DisplayName', 'Static Value');
xlabel('Time [s]');
ylabel('Angular momentum [$kg~m^2/s$]')
ylim([6.5e4,6.7e4]);
xlim([0,t(end)]);
legend();
title(['Angular momentum comparison between static theory value and ' ...
       'integrated one from orbital data'])

%%%%%%%%%%%%%%%%%%%%% WEEK 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset(RelTol = 1e-10);
[t,X_ode_j] = ode45(@(t,X_ode_j) TBP_ECEF(t,X_ode_j,mu), ...
              t_period,[X(:,1)],options);
X_ode_j = X_ode_j'; %from 1000x6 to 6x1000 to match original X vector

figure % this figure is to plot the orbits and compare to the coe one

%%%% r vector coords plotting %%%%%%%%%

%integrated ECI values
plot3(X_ode(1,:),X_ode(2,:),X_ode(3,:),'DisplayName', 'ECI','LineWidth', 2); 
hold on;
grid on;
%original coe values
plot3(X(1,:),X(2,:),X(3,:), 'r--','DisplayName', 'COE','LineWidth', 2); 
%integrated ECEF values
plot3(X_ode_j(1,:),X_ode_j(2,:),X_ode_j(3,:), 'g--','DisplayName', ...
      'ECEF','LineWidth', 2); 
title('Comparison between the three orbital representations')

% first and final r vector points
plot3(X_ode_j(1,1),X_ode_j(2,1),X_ode_j(3,1),'ok','MarkerFaceColor','y', ...
      'DisplayName','ECEF Start Point');
plot3(X_ode_j(1,end),X_ode_j(2,end),X_ode_j(3,end),'pentagram', ...
      'MarkerFaceColor','r','DisplayName','ECEF End Point ');

%skips make_earth elements from being added to the legend
L1 = legend;
L1.AutoUpdate = 'off'; 
%Generate an Earth with ECI axis unit vectors
make_earth;

%Makes a topdown (2D) view of the figure when generated, comment in to use
% view(2);
disp('ECEF first point in orbit state vector:');
disp(X_ode_j(:,1));
disp('ECEF last point in orbit state vector:');
disp(X_ode_j(:,end));




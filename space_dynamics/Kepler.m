function [E_final,theta_final] = Kepler(e, M, tol)
% Uses Newton-Ralphson method to calculate the eccentric and
% true anomaly based on supplied eccentricity and initial mean anomaly.
% Inputs: e - Eccentricity
%         M - Initial Mean Anomaly
%         tol - tolerance for the Newton-Ralphson method
% Outputs: E_final - Eccentric anomaly
%          theta_final - Mean anomaly


% Just outputting the iter value out the function is the better solution.
% The assignment however does not include kepler_iter as an output so
% obvious solution is to export it as a global
global kepler_iter


% first guess for the eccentric anomaly is the mean anomaly's value
E_guess = M;

% first iteration of the function with all variables on one side of the
% equals sign
func_loop = E_guess-e*sin(E_guess)-M;
kepler_iter = 0;
% using absolute value for the function as it can go into negatives
while abs(func_loop)>tol
    func_loop = E_guess-e*sin(E_guess)-M;
    func_p_loop = 1-e*cos(E_guess); %first derivative of func with respect to E
    delta = - (func_loop/func_p_loop);
    E_guess = E_guess + delta;
    kepler_iter = kepler_iter + 1;
end

% Outputs calculation
E_final = E_guess;
theta_final = 2*atan2(sqrt(1+e)*tan(E_final/2), sqrt(1-e));
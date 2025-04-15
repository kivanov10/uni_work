function [E_final,theta_final] = Kepler_fzero(e, M, tol)
% This is a validation function used to test the custom Kepler.m. It uses
% the inbuilt fzero function of Matlab to do so

E_guess = M;
func = @(E) E-e*sin(E)-M;
options = optimset('TolX', tol);
E_final = fzero(func, E_guess, options);
E_final_1 = mod(E_final,2*pi);

theta_final = 2*atan(sqrt((1+e)/(1-e))*tan(E_final/2));
theta_final_1 = 2*atan(sqrt((1+e)/(1-e))*tan(E_final_1/2));
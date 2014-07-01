function y = expc(x,t)
% y = expc(x,t)
%
% Bi Exponential decay + constant
%
% ARGS :
% x = argument vector [S0 T2 C]
% t = time vector in seconds
%
% RETURNS :
% y = biexponential + constant function of t


y = x(1) * exp(-t / x(2)) + x(3) * exp(-t / x(4)) + x(5);

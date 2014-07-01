function [S0,T2,a,c] = fit_nonexp(t,s,mono_est)

% t = time in msec
% s = signal

% echo spacing

tau = (t(2)-t(1))/2;

% Parameter constraints
lb = [0   0   0];
ub = [Inf Inf Inf];

% Setup optimization parameters
options = optimset('lsqcurvefit');
options.Display = 'off';
options.TolFun = 1e-3;
options.TolX = 1e-3;
options.MaxIter = 100;


% Initial parameter estimates
S0_est = mono_est(1);
T2_est = mono_est(2);
a_est  = 0.01;
c_est  = mono_est(3);

% Initial parameter estimates
% S0_est = s(1);
% T2_est = 20;
% a_est  = 0.01;
% c_est  = 0;

% Initial parameter vector
x0 = [S0_est T2_est a_est c_est];
%x0 = [a_est];



data.tau = tau;
%data.S0 = S0_est;
%data.T2 = T2_est;

% Start optimization
[x_fit,Res] = lsqcurvefit('nonexp',x0,t,s,lb,ub,options,data);

% Calculate fitted function values
s_fit = nonexp(x_fit,t,data);

% figure;plot(t,s);
% hold on;
% plot(t,s_fit,'r');
% hold off;

% Calculate return values
S0 = x_fit(1);
T2 = x_fit(2);
a  = x_fit(3);
c = x_fit(4);



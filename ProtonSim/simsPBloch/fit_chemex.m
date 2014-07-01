function [k,t_ex,s_fit,a,chi_square] = fit_chemex(Tau,T2)

%%% Tau and T2 are in msec

% Initial parameter estimates
k_est = 1;
tau_ex_est = 0;

% Setup optimization parameters
options = optimset('lsqcurvefit');
options.Display = 'off';
options.TolFun = 1e-3;
options.TolX = 1e-3;
options.MaxIter = 100;

% Initial parameter vector
x0 = [k_est tau_ex_est];

% Parameter constraints
lb = [0   0   0];
ub = [Inf Inf Inf];

% Choose tau's < 3 ms, since beyond that t2 values become wierd

clear a;
a = find(Tau<3);

% Start optimization

[x_fit,Residual] = lsqcurvefit('chemex',x0,Tau(a),1./T2(a),lb,ub,options);

k = x_fit(1)
t_ex = x_fit(2)

% Calculate fitted function values
s_fit = chemex(x_fit,Tau(a));

% figure;plot(Tau(a),1./T2(a),'b');
% hold on;
% plot(Tau(a),s_fit,'r');
% xlabel('Tau');
% ylabel('R2');

% chi square value
% chi^2 = sum ( (observed-expected)^2 / expected )

chi_square = sum( ((s_fit*1000-1000./T2(a)).^2) ./ (1000./T2(a)) )


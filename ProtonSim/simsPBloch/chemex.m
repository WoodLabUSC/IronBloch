function y = chemex(x,t)

k = x(1);
tau_ex = x(2);
tau_cp = t;

y = (k * tau_ex).*(1 - (tau_ex*(1./tau_cp)).*tanh(tau_cp./tau_ex));
function y = nonexp(x,t,data)


%x(1) = 68;

%y = data.S0.*exp(-t./data.T2).*exp( -(x(1)^(3/4)) * (data.tau^(3/4)) * (t.^(3/8)) );
y = x(1).*exp(-t./x(2)).*exp( -(x(3)^(3/4)) * (data.tau^(3/4)) * (t.^(3/8)) ) + x(4);

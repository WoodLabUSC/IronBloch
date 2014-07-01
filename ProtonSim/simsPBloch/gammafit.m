function GammaDist = gammafit(x,xdata)
% 7 Dec. 2005
% Nilesh Ghugre, CHLA
% Fit a gamma distribution function

mygamma = x(1);     % shape parameter
beta    = x(2);     % scaling parameter
mu      = 0.05; %x(3);    % location parameter
amp_scale = x(4);
R       = xdata;    
GammaDist = amp_scale * real( (((R-mu)./beta).^(mygamma-1)) .* (exp(- ((R-mu)./beta) )) ./ ( beta * gamma(mygamma) ));
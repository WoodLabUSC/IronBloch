function nextPos = randomNextStep(currPos,step,D)

%% random initialization of rand

sigma=sqrt(2*step*D);   % sigma-square = 2*step*D yields linear r2 rise over time with slope 2D.
%dx=randn(1)*sigma;  % random x, y, and z fluctuations
%dy=randn(1)*sigma;
%dz=randn(1)*sigma;

nextPos.x = currPos.x + randn(1)*sigma; % dx;
nextPos.y = currPos.y + randn(1)*sigma; % dy;
nextPos.z = currPos.z + randn(1)*sigma; % dz;

%%%%%%%%%%%%%%%%%%%%%%%%
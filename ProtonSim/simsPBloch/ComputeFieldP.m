function [delBzGridOut] = ComputeFieldP(B0,delX,FE,wetToDryWtRatio,sphereInfo,sim_box_side,sphereVolFrac,fieldGridStep, spill)
%% computes field grid based on given sphere positions
%% 17 Oct 2006
%% Nilesh Ghugre, CHLA/USC

%% Computing field due to each sphere and then calculating effective field
%% by syperposition
% B0 = 1.5;       % Hz

% Computing for an iron content of 1 mg/g, other concentrations can be
% calculated by simple scaling
% FEgrid = 1.0;  % mg iron /g tissue, dry wt.
Fe_conc = FE/wetToDryWtRatio;

%% NOTE: SQUID uses a value of 1600E-6 (SI). For volume magnetic
%% susceptibility, (4pi) * SI = (1) cgs. Hence,
%% 1600E-6 (SI) == 1600E-6/(4pi) = 1.2732E-4 (cgs) ie  (g/g FE) =
%% 1.2732E-7 (g/ mg FE) == 12.7E-8 (g/mg FE)

%% Second (more correct explanation) 
%% mass specific magnetic susceptibility 
%% ferritin = 1.6 E-6 m^3/kgFE (SI units) --> cgs units = 1.6E-6/4pi =
%% 1.2732e-007 = 12.732E-8
%% Michaelis data, ferrtin = 1.33E-6 m^3/kgFE (SI units) = 10.58E-8 (cgs)
%%                 hemosiderin = 1.1E-6 m^3/kgFE (SI units)= 8.7535E-8 (cgs)


%% delX = 4.3e-8 * Fe_conc;

% delX = 9.74e-8 * Fe_conc;       %% 4:1 hemosiderin/ferritin mix
% delX = delX * Fe_conc;

[X,Y,Z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);
%[X,Y,Z] = meshgrid(linspace(-sim_box_side/2-spill,sim_box_side/2+spill,100));

delBzGrid = zeros(size(X,1),size(X,2),size(X,3));

multFactor = (B0) * (4*pi/3) * (Fe_conc * delX / sphereVolFrac );

%tic
% loop through set of radii
parfor p = 1:length(sphereInfo.radius)
   
    radius = sphereInfo.radius(p);
    %clear X Y Z  %I think this line is causing a transparency violation, replacing - eamon 20111007
	X=[];
	Y=[];
	Z=[];
	%end update eamon 20111007
    
    [X,Y,Z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);
    
    X = X-sphereInfo.x(p);
    Y = Y-sphereInfo.y(p);
    Z = Z-sphereInfo.z(p);
    %% superposition of fields
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Using the polar equation
    [temp_TH,temp_PHI,temp_R]=cart2sph(X,Y,Z);
    %clear X Y Z  %second instance - I think this line is causing a transparency violation, replacing - eamon 20111007
	X=[];
	Y=[];
	Z=[];
	%end update eamon 20111007
    temp_theta = pi/2-temp_PHI;
    %clear temp_TH temp_PHI; %20111007 eamon fix transparency error again
	temp_TH=[];
	temp_PHI=[];
    % polar
    delBzGrid = delBzGrid + multFactor * (radius./temp_R).^3 .* (3*(cos(temp_theta).^2) - 1);
    
    %clear temp_TH temp_PHI; %20111007 eamon fix transparency error again
	temp_TH=[];
	temp_PHI=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

delBzGridOut1 = gplus(delBzGrid);

if labindex == 1
    delBzGridOut = delBzGridOut1;
else
    delBzGridOut = [];
end




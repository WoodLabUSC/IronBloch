function [signalOut]  = simulateBlochCPMG(sphereInfo,sim_box_side,patientInfo,delBzGrid,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D,hepatocyte_dim,cellBoundaryFlag,cellBiasFlag,sinusoidBoundaryFlag,B0,B0_multiplier)
%% Simulation on a cluster (parallel computation)
%% called by MriSimP
%% 26 Jun 2012
%% Eamon Doye, CHLA/USC
%% This is used after prepareField.m, this performs collision detection,
%% interpolation and signal generation
%% Proton path computation is made a little intelligent. While checking for
%% collisions, only spheres within a given radial distance are checked instead of
%% all the spheres in the volume.

%% random initialization of rand
rand('twister',sum(100*clock));

%% Select patient from patientInfo, order is increasing FE

% FE = patientInfo(patientIndx).FE;

% sphereInfo = sphereInfo;
% delBzGrid  = delBzGrid * FE;

% numSpheres = length(sphereInfo.radius);
% R_select = sphereInfo.radius;
% total_vol = sum((4*pi/3) * (R_select.^3));


%% assuming X% vol frac obtained from lysosome analysis,
%% simlulation volume in um^3 and box side in um
% sphereVolFrac = (patientInfo(patientIndx).volFrac)/100;
% sim_volume = mean(total_vol) / (sphereVolFrac);
% sim_box_side =  (sim_volume)^(1/3);

[X,Y,Z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating proton path and taking care of collisions with spheres.

%% time interval
t = 0:step:interval;

% random initial position of a proton
a = -sim_box_side/2;
b =  sim_box_side/2;

signal = zeros(3,length(t));
sigOut = zeros(3,length(t));
%% Accounting for tissue intrinsic T2 relaxation ie. under absense
%% of magnetic inhomogeneities, taken here as 50 ms (St. Pierre Blood paper). Under such
%% ideal conditions, T2 ~ T2*.
intrinsicPhase = -t./50;

gamma = 0.267e9;%42.6e6;

% gyromagnetic ratio for protons in water molecules: 42.6 MHz/T
% or (2*pi*42.6e6) rad/sec/T = 0.267 rad/sec/nT
% gamma = 2.675e8 (1/sec/T) from Jensen, Chandra, Yu,
% Quantitative model for interecho time dependence of the cpmg
% relaxation rate in iron rich gray matter, MRM 46:159-165
% (2001)

for num = (1:numProtons)
%for num = 1:numProtons    
    rand('twister',sum(100*clock));
	%for n=1:10
    protonPath.x = zeros(1,length(t));
    protonPath.y = zeros(1,length(t));
    protonPath.z = zeros(1,length(t));
    dt = step*1e-3;

    if cellBoundaryFlag == 1
        [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath2(protonPath, sphereInfo, sim_box_side, step, D, hepatocyte_dim, cellBiasFlag, sinusoidBoundaryFlag);
    else
        [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath(protonPath, sphereInfo, sim_box_side, step, D, cellBiasFlag, sinusoidBoundaryFlag, hepatocyte_dim);
    end

    delBzInt = zeros(1,length(t));
    delBzInt = interp3(X,Y,Z,delBzGrid,protonPath.x,protonPath.y,protonPath.z,'cubic');
    delBzInt(isnan(delBzInt)) = 0;

	%%Start the Bloch simulation for one proton - Eamon
	df = 0;
	gz = zeros(1,length(t));
	z = 0;
	magnetization = zeros(3,length(t));
	magnetization(:,1) = [1,0,0];
	bTot = zeros(3,length(t));
	bTot(3,:) = delBzInt + B0*B0_multiplier;  %%this line probably needs to be modded to accept other field strengths
	T1 = 1; %this value is wrong, substitute for something correct
	T2 = .05; %see St. Pierre Blood Paper
	E1 = exp(-dt/T1); %t1 decay
	E2 = exp(-dt/T2); %t2 decay
    
  	A = [E2 0 0; 0 E2 0; 0 0 E1]; %decay matrix 1
  	b = [0;0;(1-E1)]; %decay matrix 2
	
    TECounter = 1;
	for iterator = 2:length(t)
		df = delBzInt(iterator)*gamma;
		if (TECounter<=length(TE) && (t(iterator)==TE(TECounter) || (t(iterator)>TE(TECounter) && t(iterator-1)<TE(TECounter))))
		%%Add 180* rotation or whatever the specified rotation is
			%how do we specify axes
            %presently assumes perfect refocusing, so no phase cycling is
            %not necessary.  Needs to be added - Eamon 20130524
			R = yrot(pi);
			magnetization(:,iterator) = R*magnetization(:,iterator-1);
			TECounter = TECounter+1;
		else
		R = zrot(-angle(complex(bTot(1,iterator),bTot(2,iterator))))*xrot(gamma*abs(complex(bTot(1,iterator),bTot(2,iterator)))*dt)*zrot(angle(complex(bTot(1,iterator),bTot(2,iterator))));
		magnetization(:,iterator) = R * magnetization(:,iterator-1);
		magnetization(:,iterator) = zrot(df*dt) * magnetization(:,iterator);
  		magnetization(:,iterator) = zrot(gamma*bTot(iterator)*z*dt) * magnetization(:,iterator);
  		magnetization(:,iterator) = A*magnetization(:,iterator) + b.*magnetization(:,1);
 % 		if iterator == 2
 % 		 	magnetization(:,2) = [1,0,0];
 %		 	end
 % 		if iterator == 5
 % 		 	magnetization(:,1:5)
%		end
        end  %End selection between the doing 180s and regular dephasing
    end  %End Bloch For Loop

%%End the bloch portion of the simulation
    protonPath = [];
    
    %% FID
	signal = magnetization;
    sigOut = sigOut + magnetization;
    %signal = complex(magnetization(1,:),magnetization(2,:));
%end
%figure;title('simulatePnew');subplot(5,1,1);plot(t,magnetization(1,:));subplot(5,1,2);plot(t,magnetization(2,:));subplot(5,1,3);plot(t,abs(complex(magnetization(1,:),magnetization(2,:))));subplot(5,1,4);plot(t,magnetization(3,:));subplot(5,1,5);plot(t,signal);

%{
    %% T2 Single spin echo experiment
%{  Uncomment when you have confirmed that the single echo experiment works - Eamon
    for m = 1:length(TE)
  
        totPhaseSE = zeros(1,length(t));    %% total phase accrual for T2 SE

        TE_indx = round(TE(m)/step);

        totPhaseSE = totPhase;
        totPhaseSE(TE_indx:interval/step) = totPhase(TE_indx:interval/step) - 2*totPhase(TE_indx);

        signalSE(m,:) = signalSE(m,:) + exp(complex(intrinsicPhase,totPhaseSE));

    end % for m
%}    
   	delPhase = [];
    totPhase = [];

    %end
end

signalOut = sigOut;  %this is how it returns on the HPCC w/ MATLAB 2011b

if 0 %this is for r2007b only, maybe add back later
signalOut1 = gplus(signal);

if labindex == 1
    signalOut = signalOut1;
else
    signalOut = [];
end
end

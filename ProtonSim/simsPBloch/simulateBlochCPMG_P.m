function [signalOut,magnetizationOut]  = simulateBlochCPMG_P(sphereInfo, ...
	patientInfo,delBzGrid,patientIndx,simParams)
%% Simulation on a cluster (parallel computation)
%% called by MriSimP
%% 26 Jun 2012
%% Eamon Doye, CHLA/USC
%% This is used after prepareField.m, this performs collision detection,
%% interpolation and signal generation
%% Proton path computation is made a little intelligent. While checking for
%% collisions, only spheres within a given radial distance are checked instead of
%% all the spheres in the volume.

%prep parameters from data structure
if isfield(simParams,'phaseCycle')
    phaseCycle = simParams.phaseCycle;
else
    phaseCycle = 'no'
end
TE = simParams.TE ./ 2;  %In this function, TE actually means the time that the pulse is applied
fieldGridStep = simParams.fieldGridStep;
sim_box_side = simParams.simBoxSide;
hepatocyte_dim = simParams.hepatocyteDim;
numProtons = simParams.numProtons;
cellBiasFlag = simParams.cellBiasFlag;
sinusoidBoundaryFlag = simParams.sinusoidBoundaryFlag;
cellBoundaryFlag = simParams.cellBoundaryFlag;
B0 = simParams.B0;
useInstantExcitation = simParams.useInstantExcitation;
%if useInstantExcitation %this is currently unnecessary but may change if
%the engine to specify flips or RF pulses changes
    flipAngleExcite = simParams.FAExcite;
    if isstr(simParams.FAEcho)
        if simParams.FAEcho(1) == 'x'
            flipAngleEcho = str2double(simParams.FAEcho(2:end)) * flipAngleExcite;
        end
    else
        flipAngleEcho = simParams.FAEcho;
    end
%end
D = simParams.D;
spill = simParams.spill;
step = simParams.step;
interval = simParams.interval;

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
81
[X,Y,Z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating proton path and taking care of collisions with spheres.

%% time interval
t = 0:step:interval;

% random initial position of a proton
a = -sim_box_side/2;
b =  sim_box_side/2;

signal = zeros(1,length(t));
signalOut = signal;
magnetization = zeros(3,length(t));
magnetizationOut = magnetization;
%% Accounting for tissue intrinsic T2 relaxation ie. under absense
%% of magnetic inhomogeneities, taken here as 50 ms (St. Pierre Blood paper). Under such
%% ideal conditions, T2 ~ T2*.
intrinsicPhase = -t./50;

gamma = 0.267e9;%42.6e6;

% gyromagnetic ratio for protons in water molecules: 42.6 MHz/T
% or (2*pi*42.6e6) rad/sec/T = 0.267 rad/sec/nT
% gamma = 2.675e8 (1/sec/T) from Jensen, Chandra, Yu,
% Quantitative model for interecho time dependence of the cpmg
% relaxation rate in iron rich gray m81atter, MRM 46:159-165
% (2001)

parfor num = (1:numProtons)
%for num = 1:numProtons    
    rand('twister',sum(100*clock));

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
	magnetization(:,1) = [0,0,1];
	bTot = zeros(3,length(t));
	bTot(3,:) = delBzInt;  %%this line probably needs to be modded to accept other field strengths
  
	T1 = .576; %Stanisz paper, 1.5T
	if (B0==3)
		T1 = .812;
	end
	T2 = .05; %see St. Pierre Blood Paper
	if (B0==3)
		T2 = .042;
	end
	E1 = exp(-dt/T1); %t1 decay
%	E1 = 0;
	E2 = exp(-dt/T2); %t2 decay
  	A = [E2 0 0; 0 E2 0; 0 0 E1]; %decay matrix 1
  	b = [0;0;(1-E1)]; %decay matrix 2
	TECounter = 1;
	for iterator = 2:length(t)
		df = delBzInt(iterator)*gamma;
		if (TECounter<=length(TE) && (t(iterator)>=TE(TECounter) && t(iterator-1)<TE(TECounter)))
		%%Add 180* rotation or whatever the specified rotation is
			%how do we specify axes
            
            %Turn phase cycle every other echo time
            phaseSign = 1;
            if strcmp(phaseCycle,'yes') && mod((TECounter-1),2)
                phaseSign = -1;
            end
            
			R = yrot((phaseSign*flipAngleEcho)*pi/180);
			magnetization(:,iterator) = R*magnetization(:,iterator-1);
			TECounter = TECounter+1;
        else
            R = zrot(-angle(complex(bTot(1,iterator),bTot(2,iterator))))*xrot(gamma*abs(complex(bTot(1,iterator),bTot(2,iterator)))*dt)*zrot(angle(complex(bTot(1,iterator),bTot(2,iterator))));
            magnetization(:,iterator) = R * magnetization(:,iterator-1);
        	magnetization(:,iterator) = zrot(df*dt) * magnetization(:,iterator);
      		magnetization(:,iterator) = zrot(gamma*bTot(iterator)*z*dt) * magnetization(:,iterator);
            magnetization(:,iterator) = A*magnetization(:,iterator) + b.*magnetization(:,1);
            if iterator == 2
                R = yrot(flipAngleExcite*pi/180);
                magnetization(:,2) = R*magnetization(:,1);
  		 	end
  		end
	end

%%End the bloch portion of the simulation
    protonPath = [];
    
    %% Summing the signals
	magnetizationOut = magnetizationOut + magnetization;
	signalOut = signalOut + complex(magnetization(1,:),magnetization(2,:));
%end
%figure;title('simulatePnew');subplot(5,1,1);plot(t,magnetization(1,:));subplot(5,1,2);plot(t,magnetization(2,:));subplot(5,1,3);plot(t,abs(complex(magnetization(1,:),magnetization(2,:))));subplot(5,1,4);plot(t,magnetization(3,:));subplot(5,1,5);plot(t,signal);

end

if labindex ~= 1
   signalOut = [];
   magnetizationOut = [];
end

%method of output used in R2007b.  It's unclear if this will work in later versions, so temporarily removing - Eamon
%signalOut1 = gplus(signal);
%magnetizationOut1 = gplus(magnetization);
%if labindex == 1
%    signalOut = signalOut1;
%    magnetizationOut = magnetizationOut1;
%else
%    signalOut = [];
%    magnetizationOut = [];
%end
%end
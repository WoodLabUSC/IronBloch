function [simR2, simR2s]  = simulateP2(B0_multiplier_vec,sphereInfo,sim_box_side,patientInfo,delBzGrid,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D,hepatocyte_dim,cellBoundaryFlag,cellBiasFlag,sinusoidBoundaryFlag)
%% Simulation on a cluster (parallel computation)
%% called by MriSimP
%% 20 Oct 2006
%% Nilesh Ghugre, CHLA/USC
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

signal = zeros(length(B0_multiplier_vec),length(t));

for k=1:length(B0_multiplier_vec)
    signalSE(k).signal = complex(zeros(length(TE),length(t)),zeros(length(TE),length(t)));
end


%% CPMG
% echo_time = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 1 2 3 4 5 6 10 20];
% signalCPMG = complex(zeros(length(echo_time),length(t)),zeros(length(echo_time),length(t)));

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

parfor num = 1:numProtons
    
    rand('twister',sum(100*clock));

    protonPath.x = zeros(1,length(t));
    protonPath.y = zeros(1,length(t));
    protonPath.z = zeros(1,length(t));

    if cellBoundaryFlag == 1
        [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath2(protonPath, sphereInfo, sim_box_side, step, D, hepatocyte_dim, cellBiasFlag, sinusoidBoundaryFlag);
    else
        [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath(protonPath, sphereInfo, sim_box_side, step, D, cellBiasFlag, sinusoidBoundaryFlag, hepatocyte_dim);
    end

    delBzInt = zeros(1,length(t));
    delBzInt = interp3(X,Y,Z,delBzGrid,protonPath.x,protonPath.y,protonPath.z,'cubic');
    delBzInt(isnan(delBzInt)) = 0;

    clear protonPath;

    
    for multiplierIndx = 1:length(B0_multiplier_vec)

        delPhase = zeros(1,length(t));              %% phase change
        totPhase = zeros(1,length(t));      %% total phase accrual

        delPhase = gamma * step * 1e-3 * delBzInt * B0_multiplier_vec(multiplierIndx);
        totPhase = cumsum(delPhase);

        %% FID

        signal(multiplierIndx,:) = signal(multiplierIndx,:) + exp(complex(intrinsicPhase,totPhase));

        %% T2 Single spin echo experiment

        for m = 1:length(TE)

            totPhaseSE = zeros(1,length(t));    %% total phase accrual for T2 SE

            TE_indx = round(TE(m)/step);

            totPhaseSE = totPhase;
            totPhaseSE(TE_indx:interval/step) = totPhase(TE_indx:interval/step) - 2*totPhase(TE_indx);

            signalSE(multiplierIndx).signal(m,:) = signalSE(multiplierIndx).signal(m,:) + exp(complex(intrinsicPhase,totPhaseSE));

        end % for m

        clear delPhase totPhase totPhaseSE;


    end
    
    
end

signalOut1 = gplus(signal);

for multiplierIndx = 1:length(B0_multiplier_vec)
    signalSEOut1(multiplierIndx).signal = gplus(signalSE(multiplierIndx).signal);
end


if labindex == 1
    
    clear signal signalSE;
    
    signal = signalOut1;
    signalSE = signalSEOut1;

    clear signalOut1 signalSEOut1
    
    %%%%% computing R2 and R2*
    

    % step  = 0.0005; %msec        very important factor !!!
    % interval = 60; %msec
    % t = 0:step:interval;
    % TE = logspace(log(0.1)/log(10),log(30)/log(10),15);

    for multiplierIndx = 1:length(B0_multiplier_vec)

        %%%%%%%%%%%%%%%
        %%%% RESULTS
        %% T2star
        totSignal  = (1/numProtons) * sum(signal(multiplierIndx,:),1);           % complex fid
        % figure;plot(t,abs(totSignal));

        T2s_est = [1 10 20 50];

        for u=1:length(T2s_est)
            [S0(u),T2s(u),C(u),Res_1(u)] = fitexp_mc(t,abs(totSignal),T2s_est(u));
            % s_fit = expc([S0 T2s(u) C],t);
            % figure;plot(t,abs(totSignal));
            % hold on;
            % plot(t,s_fit,'r');
            % hold off;
        end

        [P,Q] = min(Res_1);     % choose fit with least residual

        simR2s(multiplierIndx).R2s = 1000/(T2s(Q));
        simR2s(multiplierIndx).S0 = S0(Q);
        simR2s(multiplierIndx).C = C(Q);
        simR2s(multiplierIndx).Res = Res_1(Q);
        simR2s(multiplierIndx).t = t;
        simR2s(multiplierIndx).absSignal = abs(totSignal);
        
        %% T2 SE

        %% TE defined before main for-loop
        % signalSE = zeros(numProtons,length(t));
        % totPhaseSE = zeros(numProtons,length(t));    %% total phase accrual for T2 SE

        %% extracting echoes

        for m = 1:length(TE)

            TE_indx = round(TE(m)/step);

            totSignalSE = (1/numProtons) * signalSE(multiplierIndx).signal(m,:);

            if(m==length(TE))       %% checking if end of observation time frame is reached
                [signalSEecho(m),signalSEindx(m)] = max(real(totSignalSE(2*TE_indx-100:2*TE_indx)));
            else
                [signalSEecho(m),signalSEindx(m)] = max(real(totSignalSE(2*TE_indx-100:2*TE_indx+100)));
            end

        end

        % figure;plot(t,abs(totSignalSE));

        tEcho = TE * 2;

        % figure;plot(tEcho,signalSEecho);

        T2_est = [1 10 20 50 80 120];


        for u = 1:length(T2_est)
            [S0(u),T2(u),C(u),Res_2(u)] = fitexp_mc(tEcho,signalSEecho,T2_est(u));
            %         s_fit = expc([S0 T2(u) C],tEcho);
            %         figure;plot(tEcho,signalSEecho,'o-');
            %         hold on;
            %         plot(tEcho,s_fit,'r');
            %         hold off;
        end

        [P,Q] = min(Res_2);

        simR2(multiplierIndx).R2 = 1000/T2(Q);
        simR2(multiplierIndx).S0 = S0(Q);
        simR2(multiplierIndx).C = C(Q);
        simR2(multiplierIndx).Res = Res_2(Q);
        simR2(multiplierIndx).tEcho = tEcho;
        simR2(multiplierIndx).signalSEecho = signalSEecho;

        clear totSignal totSignalSE signalSEecho signalSEindx
    end    

    
else
    simR2 = [];
    simR2s = [];
end



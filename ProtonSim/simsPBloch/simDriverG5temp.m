%% Simulation driver file
%% 3 November 2006
%% Nilesh Ghugre, CHLA/USC
%% Modified for computed gamma function

%% Patients created by PreparePatientInfo.m
% patientDir = 'G:\Nilesh\Monte Carlo\DCT-Oct2006\initial simsP 3 Results\FE 0-20';
% load patientInfo
%%%%%%%%%%%%%%%%%%

%%%% Dummy numbers created PreparePatientInfo.m
load dummyInfo

%% Choose cases.

%% For FE = 0.5, 1:20
% patientDir = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 initial simsP 3 Results/FE 0-20 Sigma 4 Hep 20um';
% k=1;
% patientInfo(k) = dummyInfo(1);
% k=k+1;
% for i = 2:2:40
%     patientInfo(k) = dummyInfo(i);
%     k=k+1;
% end

% %% For FE = 21:40
% patientDir = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 initial simsP 3 Results/FE 20-40 Sigma 4 Hep 20um';
% k=1;
% for i = 42:2:80
%     patientInfo(k) = dummyInfo(i);
%     k=k+1;
% end

%% For FE = 0.5, 1:40
patientDir = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 initial simsP 3 Results/FE 0-40 Sigma 4 Hep 20um';
k=1;
%patientInfo(k) = dummyInfo(1);
k=k+1;
%for i = 2:2:80
%    patientInfo(k) = dummyInfo(i);
%    k=k+1;
%end

clear dummyInfo i k
%%%%%%%%%%%%%%%%%%


sim_box_side = 40;     % um
hepatocyte_dim = 10;    % in um

fieldGridStep = 0.5;    % um for every grid point, 0.5 is nyquist for single sphere, however
spill = 2;      % um, to avoid NAN's during interpolation of boundary points

numWorkersField = 32;

cellBiasFlag = 1;   %% 0 for uniform, 1 for gaussian bias

%% Select patient from patientInfo, order is increasing FE
for patientIndx = 1:length(patientInfo)
    
    prepareFieldP
    % SphereVisualize(sphereInfo.radius,sphereInfo.x,sphereInfo.y,sphereInfo.z,sim_box_side);
    patientID = patientInfo(patientIndx).id;
    mkdir(sprintf('%s/BiasedDistribution/Restricted/%s',patientDir,num2str(patientID)));
    
    if cellBiasFlag == 1
     
        save(sprintf('%s/BiasedDistribution/Restricted/%s/sphereInfo.mat',patientDir,num2str(patientID)), 'sphereInfo');
        save(sprintf('%s/BiasedDistribution/Restricted/%s/delBzGridP.mat',patientDir,num2str(patientID)), 'delBzGridP');
        save(sprintf('%s/BiasedDistribution/Restricted/%s/params.mat',patientDir,num2str(patientID)), 'sim_box_side','hepatocyte_dim');
    else
        save(sprintf('%s/UniformDistribution/Restricted/%s/sphereInfo.mat',patientDir,num2str(patientID)), 'sphereInfo');
        save(sprintf('%s/UniformDistribution/Restricted/%s/delBzGridP.mat',patientDir,num2str(patientID)), 'delBzGridP');
        save(sprintf('%s/UniformDistribution/Restricted/%s/params.mat',patientDir,num2str(patientID)), 'sim_box_side','hepatocyte_dim');

    end

    clear sphereInfo delBzGridP;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NOTE the same sphereInfo and delBzGridP can be used for restricted an
%% unrestricted motion

step  = 0.0005; %msec        very important factor !!!
interval = 60; %msec
D = 0.76; % micron^2/msec, this is the value for human liver

numProtons = 1000;
%% flags to decide restricted (1) or unrestricted (0) proton motion
% cellBoundaryFlag = 1;

simResultsFlag = 1;
CPMGsimResultsFlag = 0;

numWorkersSim = 32;

%% Select patient from patientInfo
for patientIndx = 1:length(patientInfo)
    
    patientID = patientInfo(patientIndx).id;

    for cellBiasFlag = 1

        if cellBiasFlag == 1
            load(sprintf('%s/BiasedDistribution/Restricted/%s/sphereInfo.mat',patientDir,num2str(patientID)));
            load(sprintf('%s/BiasedDistribution/Restricted/%s/delBzGridP.mat',patientDir,num2str(patientID)));
        else
            load(sprintf('%s/UniformDistribution/Restricted/%s/sphereInfo.mat',patientDir,num2str(patientID)));
            load(sprintf('%s/UniformDistribution/Restricted/%s/delBzGridP.mat',patientDir,num2str(patientID)));

        end

        for cellBoundaryFlag = 1

            MriSimP
%             save signal signal
%             save signalCPMG signalCPMG

            if cellBiasFlag == 1

                if cellBoundaryFlag == 1
                    save(sprintf('%s/BiasedDistribution/Restricted/%s/signal.mat',patientDir,num2str(patientID)), 'signal');
                    save(sprintf('%s/BiasedDistribution/Restricted/%s/signalSE.mat',patientDir,num2str(patientID)), 'signalSE');
                    %save(sprintf('%s/BiasedDistribution/Restricted/%s/signalCPMG.mat',patientDir,num2str(patientID)), 'signalCPMG');
                    save(sprintf('%s/BiasedDistribution/Restricted/%s/numProtons.mat',patientDir,num2str(patientID)), 'numProtons');
                else
                    save(sprintf('%s/BiasedDistribution/Unrestricted/%s/signal.mat',patientDir,num2str(patientID)), 'signal');
                    save(sprintf('%s/BiasedDistribution/Unrestricted/%s/signalSE.mat',patientDir,num2str(patientID)), 'signalSE');
                    %save(sprintf('%s/BiasedDistribution/Unrestricted/%s/signalCPMG.mat',patientDir,num2str(patientID)), 'signalCPMG');
                    save(sprintf('%s/BiasedDistribution/Unrestricted/%s/numProtons.mat',patientDir,num2str(patientID)), 'numProtons');

                end
            else
                if cellBoundaryFlag == 1
                    save(sprintf('%s/UniformDistribution/Restricted/%s/signal.mat',patientDir,num2str(patientID)), 'signal');
                    save(sprintf('%s/UniformDistribution/Restricted/%s/signalSE.mat',patientDir,num2str(patientID)), 'signalSE');
                    %save(sprintf('%s/UniformDistribution/Restricted/%s/signalCPMG.mat',patientDir,num2str(patientID)), 'signalCPMG');
                    save(sprintf('%s/UniformDistribution/Restricted/%s/numProtons.mat',patientDir,num2str(patientID)), 'numProtons');
                else
                    save(sprintf('%s/UniformDistribution/Unrestricted/%s/signal.mat',patientDir,num2str(patientID)), 'signal');
                    save(sprintf('%s/UniformDistribution/Unrestricted/%s/signalSE.mat',patientDir,num2str(patientID)), 'signalSE');
                    %save(sprintf('%s/UniformDistribution/Unrestricted/%s/signalCPMG.mat',patientDir,num2str(patientID)), 'signalCPMG');
                    save(sprintf('%s/UniformDistribution/Unrestricted/%s/numProtons.mat',patientDir,num2str(patientID)), 'numProtons');

                end

            end

            % clear signal signalSE signalCPMG;
            clear signal signalSE;

        end     % cellBoundaryFlag

    end     % cellBiasFlag

end     % patientIndx

save(sprintf('%s/BiasedDistribution/Restricted/simR2s.mat',patientDir), 'simR2s');
save(sprintf('%s/BiasedDistribution/Restricted/simR2.mat',patientDir), 'simR2');


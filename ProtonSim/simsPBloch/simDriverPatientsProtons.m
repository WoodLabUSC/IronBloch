%% Simulation driver file
%% 1 October 2007
%% Nilesh Ghugre, CHLA/USC
%% Modified for computed gamma function
%% and inter-cell iron distribution
%% Called by BatchScriptSim.m


%% NOTE the same sphereInfo and delBzGridP can be used for restricted an
%% unrestricted motion

simResultsFlag = 1;
CPMGsimResultsFlag = 0;


%% Select patient from patientInfo
%for patientIndx = selectedPatientIndx %1:length(patientInfo)
    
    % patientID = patientInfo(patientIndx).id

    load(sprintf('%s/sphereInfo.mat',patientDirRuns), 'sphereInfo');
    load(sprintf('%s/delBzGridP.mat',patientDirRuns), 'delBzGridP');
    
    %% Using specified multipliers

    delBzGridP = delBzGridP*(B0_multiplier*delX_multiplier/wetToDryWtRatio_multiplier);
    
    MriSimP

    %% parameters for saving results
    B0_sim = num2str(B0*B0_multiplier);
    D_sim  = num2str(D);
    wetToDryWtRatio_sim = num2str(wetToDryWtRatio*wetToDryWtRatio_multiplier);
    delX_sim = num2str((delX*delX_multiplier)*1e8);
    numProtons_sim = num2str(numProtons);
    
    
    if cellBoundaryFlag == 1
        patientDirSim = [patientDirRuns '/B0-' B0_sim '/D-' D_sim '/WTD-' wetToDryWtRatio_sim '/delX-' delX_sim '/Protons-' numProtons_sim '/Restricted' '/Run' num2str(runsIndxSim)];
    else
        patientDirSim = [patientDirRuns '/B0-' B0_sim '/D-' D_sim '/WTD-' wetToDryWtRatio_sim '/delX-' delX_sim '/Protons-' numProtons_sim '/Unrestricted' '/Run' num2str(runsIndxSim)];
    end

    mkdir(patientDirSim);

    save(sprintf('%s/signal.mat',patientDirSim), 'signal');
    save(sprintf('%s/signalSE.mat',patientDirSim), 'signalSE');
    save(sprintf('%s/numProtons.mat',patientDirSim), 'numProtons');
    save(sprintf('%s/simR2.mat',patientDirSim), 'simR2');
    save(sprintf('%s/simR2s.mat',patientDirSim), 'simR2s');

    disp(sprintf('FE = %f',patientInfo.FE));
    disp(sprintf('simR2s = %f',simR2s));
    disp(sprintf('simR2 = %f',simR2));
    disp('--------------------------------');
   
    
    clear signal signalSE simR2 simR2s;


% end     % patientIndx



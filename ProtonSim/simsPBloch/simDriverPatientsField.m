%% Simulation driver file
%% 1 October 2007
%% Nilesh Ghugre, CHLA/USC
%% Modified for computed gamma function
%% and inter-cell iron distribution
%% Called by BatchScriptSim.m



% %% Select patient from patientInfo, order is increasing FE
% for patientIndx = 1:length(patientInfo)
    
    prepareFieldP
    % SphereVisualize(sphereInfo.radius,sphereInfo.x,sphereInfo.y,sphereInfo.z,sim_box_side);

    mkdir(patientDirRuns);
    save(sprintf('%s/sphereInfo.mat',patientDirRuns), 'sphereInfo');
    save(sprintf('%s/delBzGridP.mat',patientDirRuns), 'delBzGridP');
    save(sprintf('%s/params.mat',patientDirRuns), 'sim_box_side','hepatocyte_dim');
    save(sprintf('%s/patientInfo.mat',patientDirRuns), 'patientInfo');


    % clear sphereInfo delBzGridP;
    
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




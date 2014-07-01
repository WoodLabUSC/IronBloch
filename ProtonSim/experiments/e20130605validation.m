%% Eamon Doyle, CHLA/USC
%% Slim Driver Script Incorporating Better Job Setup Mechanism

[jobParams simParams] = simInitParams();

generateFieldFlag = 0;

FE = [2:6:14];
jobParams.numWorkers = 6;
jobParams.emailUpdatesTo = [];
simParams.simSelect = 'CPMG';
simParams.useInstantExcitation = 1;

%%%%% Display variable settings before execution %%%%%%

disp('-------------SETTINGS-------------------');
disp(sprintf('wetToDryWtRatio =             %s',num2str(simParams.wetToDryWtRatio)));
disp(sprintf('delX =                        %s',num2str(simParams.delX)));
disp(sprintf('B0 =                          %s',num2str(simParams.B0)));
disp(sprintf('cellBiasFlag =                %s',num2str(simParams.cellBiasFlag)));
disp(sprintf('distributionVariabilityType = %s',num2str(simParams.ironDist.distributionVariabilityType)));
disp(sprintf('numWorkersField =             %s',num2str(jobParams.numWorkers)));
%disp(sprintf('numWorkersSim =               %s',num2str(jobParams.numWorkersSim)));
disp(sprintf('numProtons =                  %s',num2str(simParams.numProtons)));
disp(sprintf('cellBoundaryFlag =            %s',num2str(simParams.ironDist.cellBoundaryFlag)));
disp(sprintf('D =                           %s',num2str(simParams.D)));
disp(sprintf('B0_multiplier =               %s',num2str(simParams.B0_multiplier)));
disp(sprintf('wetToDryWtRatio_multiplier =  %s',num2str(simParams.wetToDryWtRatio_multiplier)));
disp(sprintf('delX_multiplier =             %s',num2str(simParams.delX_multiplier)));
disp(sprintf('useInstantExcitation =        %s',num2str(simParams.useInstantExcitation)));


%% 28 September  2007
%% Nilesh Ghugre, CHLA/USC
%% Batch script for simulation

whichMachine = 1;    %%  1 : Zaphod
                     %%  2 : G5

%% output directory
% patientDir = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 initial simsP 3 Results/Patients3 FE Hep 20um NN GRP-10000 M9P2 WDR-4.1/ProtonDependence'


switch whichMachine

    case 1      %% for Zaphod
        resultsDir = 'G:\Nilesh-MC-Results\G5 simsP patients final Results';
        FileDependenciesString = 'G:\Nilesh\Monte Carlo\DCT-Oct2006\G5 sims\G5 simsP patients final';    

    case 2      %% for G5
        resultsDir = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 simsP patients final Results';
        FileDependenciesString = '/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 simsP patients final';    
end


%%%%%%%%%%%%%%%%%%%%% Nominal parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Field will be generated and saved with these parameters
%% Alterations can be done during proton simulation since they are only
%% multiplying factors.

wetToDryWtRatio = 4.1;  %% WDR, from Zuyderhoudt et al.
%% previous simulations were performed at 3.5, unless specified.

%% To get the iron concentration, keeping approximately the same multiplying 
%% factor as the patients for similar R2 range, (FE*delX/volFrac).
delX = 9.74e-8;       %% 4:1 hemosiderin/ferritin mix


%% All fields will be generated for 1T, appropriate field multipliers will
%% be used during MRI sim below.
B0 = 1.5;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% Parameters that can be altered for geometry generation

cellBiasFlag = 4;   %% 0 for uniform,
%% 1 for inter-cellular gaussian bias,
%% 2 for inter-cellular patient specific anisotropy,
%% 3 for patient specific inter-cellular anisotropy
%% along with inter-particle distance measures.
%% 4, same as 3 but with sinusoids.
%% 5, same as 4 but here only sinusoids are filled, when spheres cannot be
%% accommodated, they are spilled into the hepatocyes. (simulating Type 4
%% hemochromatosis).

if cellBiasFlag == 0
    cellSigma = [];     %% specified but not used !
elseif cellBiasFlag == 1
    cellSigma = 10;     %% std dev of cell-iron pdf
elseif cellBiasFlag == 2
    cellSigma = [];     %% specified but not used !
elseif cellBiasFlag == 3
    cellSigma = [];     %% specified but not used !
elseif cellBiasFlag == 4
    cellSigma = [];     %% specified but not used !    
elseif cellBiasFlag == 5
    cellSigma = [];     %% specified but not used !        
end


%% Dont change for now...
NNfactor = 1;   %% nearest neighbor conversion factor, distance = NNfactor* (from NN histogram)
NNfactorS = 1;    %% sinusoids


%% Distribution Variability Type

distributionVariabilityType = 1;
%% 3, 4, 5, 6 are available but are not going to be used.

%% Variability is std. dev or RMS error
switch distributionVariabilityType

    case 1
        %% 1 : Mean distribution parameters
        Size_Spread_FE_Variability = 0;
        Anisotropy_Spread_FE_Variability = 0;
        Anisotropy_Shape_FE_Variability = 0;
        NN_Shape_FE_Variability = 0;

        
    case 2
        %% 2 : Random
        Size_Spread_FE_Variability = 0.013;      %% stdev
        Anisotropy_Spread_FE_Variability = 0.064;
        Anisotropy_Shape_FE_Variability = 0.412;
        NN_Shape_FE_Variability = 1.221;
        

    case 3
        %% 3 : Object Size Variability, Spread
        Size_Spread_FE_Variability = 0.013;      %% multiples of stdev = 0.012896
        Anisotropy_Spread_FE_Variability = 0;
        Anisotropy_Shape_FE_Variability = 0;
        NN_Shape_FE_Variability = 0;


    case 4
        %% 4 : Nearest Neighbor (NN) Variability, Shape
        NN_Shape_FE_Variability = 1.221;
        Size_Spread_FE_Variability = 0;
        Anisotropy_Spread_FE_Variability = 0;
        Anisotropy_Shape_FE_Variability = 0;


    case 5
        %% 5 : Anisotropy Variability, Spread
        Anisotropy_Spread_FE_Variability = 0.064;
        Size_Spread_FE_Variability = 0;
        Anisotropy_Shape_FE_Variability = 0;
        NN_Shape_FE_Variability = 0;


    case 6
        %% 6 : Anisotropy Variability, Shape
        Anisotropy_Shape_FE_Variability = 0.412;
        Size_Spread_FE_Variability = 0;
        Anisotropy_Spread_FE_Variability = 0;
        NN_Shape_FE_Variability = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% patientID = {'rw';'gh';'he';'sc';'gz';'vm';...
%              'vl';'bb';'mf';'mb';'lq';'yp';...
%              'bs';'bj';'ts';'na';'bl';'cl';'sp';'es'};
% FE = [1.3 1.4 4.4 4.6 5.9 7.8 ...
%       12.7 13.4 14.8 16.2 16.6 19.2 ... 
%       23.1 25.5 29.0 29.6 30.0 32.9 35.4 57.8];  
  

FE_1 = [0.5 1 2.5:2.5:12.5 15:5:60];
FE_2 = [2:4 6:9 11:14 16:19 21:24 26:29 31:34 36:39 41:44 46:49 51:54 56:59];

FE = sort([FE_1 FE_2]);



%% new run...

FE = [61:80];



sim_box_side = 80;     % um
hepatocyte_dim = 20;    % in um

fieldGridStep = 0.5;    % um for every grid point, 0.5 is nyquist for single sphere, however
spill = 2;      % um, to avoid NAN's during interpolation of boundary points




%%%%%%%%%%  MRI SIMULATION  %%%%%%%%%%

%jobManager = 'MRISIM-JM';
numWorkersSim = 36;

numProtons = 5000;

step  = 0.0005; %msec        very important factor !!!
interval = 60; %msec, if interval needs to be changed, please have a look 
% at MriSimP and simulateP also since constants related to it may be hard 
% coded when echos are collected. This may be specially true for CPMG sims. 

patientIndx = 1;        %% dummy

%%%%%%%%%%%% Parameters that can be altered for proton motion simulation

%% Multipliers for field, wetToDryWtRatio and delX (from nominal)
%% Change the term in numerator to whatever value needs to be interrogated.
%% For example, for 3T, B0_multiplier = 3/1.5 = 2.
B0_multiplier= (B0)/B0;
wetToDryWtRatio_multiplier= (wetToDryWtRatio)/wetToDryWtRatio;
delX_multiplier  = (delX)/delX;

Dfactor = 1/2;
D = 0.76*Dfactor; % micron^2/msec, this is the value for human liver

cellBoundaryFlag = 1;   %% 0:off (unrestricted diffusion), 1:on (restricted diffusion)

%% 0:off (sinusoid unrestricted), proton can enter sinusoid
%% 1:on (sinusoid restricted), proton cannot enter sinusoid
sinusoidBoundaryFlag = 0;   %% there is no need to turn on boundaries.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%  GEOMETRY/FIELD GENERATION  & SIMULATION %%%%%%%%%%

%jobManager = 'MRISIM-JM';
numWorkersField = 40;

generateFieldFlag = 0;  %% 0: dont generate, 1: generate


%%%%% Display variable settings before execution %%%%%%

disp('-------------SETTINGS-------------------');
disp(sprintf('wetToDryWtRatio =             %s',num2str(wetToDryWtRatio)));
disp(sprintf('delX =                        %s',num2str(delX)));
disp(sprintf('B0 =                          %s',num2str(B0)));
disp(sprintf('cellBiasFlag =                %s',num2str(cellBiasFlag)));
disp(sprintf('distributionVariabilityType = %s',num2str(distributionVariabilityType)));
disp(sprintf('numWorkersField =             %s',num2str(numWorkersField)));
disp(sprintf('numWorkersSim =               %s',num2str(numWorkersSim)));
disp(sprintf('numProtons =                  %s',num2str(numProtons)));
disp(sprintf('cellBoundaryFlag =            %s',num2str(cellBoundaryFlag)));
disp(sprintf('D =                           %s',num2str(D)));
disp(sprintf('B0_multiplier =               %s',num2str(B0_multiplier)));
disp(sprintf('wetToDryWtRatio_multiplier =  %s',num2str(wetToDryWtRatio_multiplier)));
disp(sprintf('delX_multiplier =             %s',num2str(delX_multiplier)));

if generateFieldFlag == 0
    
    disp(sprintf('Generate field?               No'));
    
else
    
    disp(sprintf('Generate field?               Yes')); 
    
end

disp('----------------------------------------');

answer = input('Are settings OK? (y/n):','s')
disp('----------------------------------------');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(answer,'y')

    disp('------------- Beginning Execution ---------------------------');

    for FEselect = 1:length(FE)

        [patientInfo,patientDir] = PreparePatientInfo(resultsDir,cellBiasFlag,cellSigma,FE(FEselect),distributionVariabilityType,Size_Spread_FE_Variability,Anisotropy_Spread_FE_Variability,Anisotropy_Shape_FE_Variability,NN_Shape_FE_Variability);
        patientDir
        % patientInfo

        for runsIndx = 1

            patientDirRuns = [patientDir '/Run' num2str(runsIndx)];

            if generateFieldFlag == 1
                simDriverPatientsField
            end


            for runsIndxSim = 1     %% for sims

                %% individual multipliers for B0, wetToDryWtRatio and delX will
                %% be evaluated
                simDriverPatientsProtons


                %             %% combined multiplier will be evaluated, we will evaluate for
                %             %% multiple fields simultaneously. B0, wetToDryWtRatio and delX
                %             %% are only multipliers, so evaluating for any one will be the
                %             %% same as evaluating for the others. Multiple fields will be
                %             %% computed simultaneously for the same proton paths.
                %             %% Enter the fields that need to be examined in the numerator.
                %             %             B0_multiplier_vec = [0.1 0.25 0.4 0.5 0.6 0.75...
                %             %                                  0.9 1.0 1.1 1.25 1.5 1.75 2.0 ...
                %             %                                  2.25 2.5 2.75 3.0 5.0 7.0 9.4 ...
                %             %                                  ]./B0;
                %
                %             B0_multiplier_vec = [0.25 0.5 0.75...
                %                 1.0 1.25 1.5 1.75 2.0 ...
                %                 2.5 3.0 5.0 7.0 ...
                %                 ]./B0;
                %
                %             %             B0_multiplier_vec = [0.5 0.75...
                %             %                                  1.0 1.25 1.5 2.0 ...
                %             %                                  3.0 7.0 ...
                %             %                                  ]./B0;
                %
                %             %             B0_multiplier_vec = [0.5 ...
                %             %                                  1.0 1.5 ...
                %             %                                  3.0 7.0 ...
                %             %                                  ]./B0;
                %             simDriverPatientsProtons2




            end

        end

    end

else
    
    disp('Execution aborted.');
    
end







% %%%%%%%%%%  MRI SIMULATION for sepatate runs %%%%%%%%%%
% 
% %jobManager = 'MRISIM-JM';
% numWorkersSim = 30;
% 
% numProtons = 5000;
% 
% step  = 0.0005; %msec        very important factor !!!
% interval = 60; %msec, if interval needs to be changed, please have a look 
% % at MriSimP and simulateP also since constants related to it may be hard 
% % coded when echos are collected. This may be specially true for CPMG sims. 
% 
% patientIndx = 1;        %% dummy
% 
% %%%%%%%%%%%% Parameters that can be altered for proton motion simulation
% 
% %% Multipliers for field, wetToDryWtRatio and delX (from nominal)
% %% Change the term in numerator to whatever value needs to be interrogated.
% %% For example, for 3T, B0_multiplier = 3/1.5 = 2.
% B0_multiplier= (B0)/B0;
% wetToDryWtRatio_multiplier= (wetToDryWtRatio)/wetToDryWtRatio;
% delX_multiplier  = (delX)/delX;
% 
% D = 0.76; % micron^2/msec, this is the value for human liver
% 
% cellBoundaryFlag = 1;   %% 0:off (unrestricted diffusion), 1:on (restricted diffusion)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for FEselect = 1:length(FE)
% 
%     [patientInfo,patientDir] = PreparePatientInfo(resultsDir,cellBiasFlag,cellSigma,FE(FEselect),distributionVariabilityType,Size_Spread_FE_Variability,Anisotropy_Spread_FE_Variability,Anisotropy_Shape_FE_Variability,NN_Shape_FE_Variability);
%     patientDir
%   
%     for runsIndx = 1        %% for field
% 
%         patientDirRuns = [patientDir '/Run' num2str(runsIndx)];
% 
%         for runsIndxSim = 1     %% for sims
%         
%             simDriverPatientsProtons
%         
%         end
%     end
%     
%     
% end






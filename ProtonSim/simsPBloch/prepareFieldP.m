%%% Prepare magnetic field grid from a collection of spheres placed based
%%% on patient lysosome distribution. This is a modification of
%%% prepareField.m and uses DCT on the cluster.
%%% 17 Oct 2006
%%% Nilesh Ghugre, CHLA/USC

% sphereVolFrac = (patientInfo(patientIndx).volFrac)/100;

%% Specified number of spheres are generated, given the volume fraction.
%% Simulation box size depends on number of spheres.
% numSpheres = 500;
% [sphereInfo, sim_box_side, sphereVolFrac] = GetSpheres1(patientInfo, patientIndx, numSpheres)
% fieldGridStep = 0.2;    % um for every grid point, for single use 0.03
% spill = 2;      % um, to avoid NAN's during interpolation of boundary points, for single use 0.5


%% For a given simulation box and sphere volume fraction, generate sphere
%% positions. Number of spheres depends on volume fraction.
disp('----------getting spheres------------');
tic;[sphereInfo, numSpheres, sphereVolFrac] = GetSpheres2(NNfactor,NNfactorS,cellSigma,patientInfo, patientIndx, sim_box_side,hepatocyte_dim,cellBiasFlag);toc
numSpheres
disp('-------------------------------------');
% SphereVisualize(sphereInfo.radius,sphereInfo.x,sphereInfo.y,sphereInfo.z,sim_box_side);
%% numSpheres = length(sphereInfo.radius);

%% Regular matlab computation
% tic;[delBzGrid] = ComputeField(sphereInfo,sim_box_side,sphereVolFrac);toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison with previous sims in patients

% load('G:\Nilesh\Monte Carlo\MonteCarlo_Aug2006\sims\patient-5\sphereInfo');
% load('G:\Nilesh\Monte Carlo\MonteCarlo_Aug2006\sims\patient-5\delBzGrid');
% load('G:\Nilesh\Monte Carlo\MonteCarlo_Aug2006\sims\patient-5\patientInfo');
% 
% %% Select patient from patientInfo, order is increasing FE
% patientIndx = 5;
% 
% numSpheres = length(sphereInfo.radius);
% R_select = sphereInfo.radius;
% total_vol = sum((4*pi/3) * (R_select.^3));
% 
% %% assuming X% vol frac obtained from lysosome analysis,
% %% simlulation volume in um^3 and box side in um
% sphereVolFrac = (patientInfo(patientIndx).volFrac)/100;
% sim_volume = mean(total_vol) / (sphereVolFrac);
% sim_box_side =  (sim_volume)^(1/3)
% 
% fieldGridStep = 0.2;    % um for every grid point, for single use 0.03
% spill = 2;      % um, to avoid NAN's during interpolation of boundary points, for single use 0.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computation using DCT
% tic;[delBzGrid] = ComputeFieldP(sphereInfo,sim_box_side,sphereVolFrac);toc
%%%%%%% parallel computing
disp('-----------------Computing field-----------------');
tic
%sched = findResource('scheduler','type','jobmanager','name','MRISIM-JM','LookupURL','cluster1');
sched = findResource('scheduler','type','jobmanager','name',getJobmanagerInfo('jm')','LookupURL',getJobmanagerInfo('host'));
pjob = createParallelJob(sched);
set(pjob, 'FileDependencies', {'ComputeFieldP.m'});
% set(pjob,'FileDependencies',{'/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 initial simsP 3 FE 0-20'});

set(pjob, 'MaximumNumberOfWorkers', numWorkersField);

% set(pjob, 'Timeout', 60);   % sec

fieldTask = createTask(pjob, @ComputeFieldP, 1, {B0,delX,patientInfo.FE,wetToDryWtRatio,sphereInfo,sim_box_side,sphereVolFrac,fieldGridStep, spill});

submit(pjob);
waitForState(pjob,'finished');

outputP = getAllOutputArguments(pjob);
toc
%eamon added to debug%
errmsgs = get(pjob.Tasks, {'ErrorMessage'});
nonempty = ~cellfun(@isempty, errmsgs);
celldisp(errmsgs(nonempty));
%end eamon debug%
delBzGridP = outputP{1};

destroy(pjob);
clear outputP fieldTask pjob sched;

disp('--------------------------------------------------');


%%%%%%%%%%

%% optional use of a function
% delBzGrid = generateField(B0,sim_box_side,fieldGridStep,sphereInfo,totalSphereVol,sphereVolFrac,delX); 
%% later could try to write in C and convert to a mex function, not sure
%% whether it will save time...
                                                                                                                                                                                                                              
% % % % % %% Looking at 2D field profile
% fieldGridStep = 0.2;    % um for every grid point, for single use 0.03
% spill = 2;      % um, to avoid NAN's during interpolation of boundary points, for single use 0.5
% [X1,Z1] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);
% layerindx = 12;
% figure; mesh(X1,Z1,squeeze(delBzGrid(:,layerindx,:)));view(180,160); %axis off;
% figure; mesh(X1,Z1,squeeze(delBzGridP(:,layerindx,:)));view(180,160); %axis off;
% figure; mesh(X1,Z1,squeeze(delBzGrid(:,layerindx,:))-squeeze(delBzGridP(:,layerindx,:)));view(180,160); %axis off;
% figure; mesh(X1,Z1, (squeeze(delBzGrid(:,layerindx,:))-squeeze(delBzGridP(:,layerindx,:)))./ (squeeze(delBzGrid(:,layerindx,:))*100) );view(180,160); %axis off;
% 
% % max and min % difference
% max(max(max((delBzGrid-delBzGridP)./(delBzGrid*100))))
% min(min(min((delBzGrid-delBzGridP)./(delBzGrid*100))))


%% are of the order of E-012%, extremely negligible, probably caused by
%% floating point operation errors.

% %% Cartesian, there seems to be some conversion factor compared to polar
% tic
% delBzGrid = delBzGrid + (currSphereVol/totalSphereVol) * (1/3)*  B0 * (delX / sphereVolFrac ) * ((radius).^3) * ( (2*Z.^2 - X.^2 -Y.^2) ./ ((X.^2 + Y.^2 + Z.^2).^(5/2)) );
% toc











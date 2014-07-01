function [sphereInfo, numSpheres, sphereVolFrac] = GetSpheres2(NNfactor,NNfactorS,cellSigma,patientInfo, patientIndx, sim_box_side,hepatocyte_dim,cellBiasFlag)
%% 24 Oct 2006
%% Nilesh Ghugre, CHLA/USC
%% For a given simulation box and sphere volume fraction, generate sphere
%% positions. Number of spheres depends on volume fraction.
%% Sphere volume fraction will be approximately equal to that specified in
%% patientInfo depending on the order in which the simulation volume is
%% filled.
%% 16 April, 2007, modified for varying cellBiasFlag

%% random initialization of rand
rand('twister',sum(100*clock));

R_highres = linspace(min(patientInfo(patientIndx).r),max(patientInfo(patientIndx).r),100);

sphereVolFrac = (patientInfo(patientIndx).volFrac)/100;     %% required

totalSphereVol =0;
simVolume = sim_box_side^3;
currSphereVolFrac = totalSphereVol/simVolume;

R_select = [];
indx = 1;
while currSphereVolFrac < sphereVolFrac

    rand_num = rand(1,1);
    [Y,I] = min(abs(patientInfo(patientIndx).dist_est_cumsum-rand_num));
    R_select(indx) = R_highres(I);

    totalSphereVol = totalSphereVol + (4*pi/3) * (R_select(indx)^3);
    currSphereVolFrac = totalSphereVol/simVolume;
        
    indx = indx +1;
end

%% from generated spheres
sphereVolFrac = currSphereVolFrac;
numSpheres = length(R_select);

% figure; hist(R_select,100);

%% Sorting sphere radii so that largest is the first.
% R_select = sort(R_select,'descend');

%% memory allocations for sphere centers
xTemp = zeros(1,numSpheres);
yTemp = zeros(1,numSpheres);
zTemp = zeros(1,numSpheres);


if cellBiasFlag == 0
    %% M2 uniformly distributes spheres in the medium
    [xTemp,yTemp,zTemp] = placeSpheresM2(xTemp,yTemp,zTemp,sim_box_side,R_select);
    cellIndx = [];
    cellFillCounter = [];
    R_select_orig = [];
    cellSpheres = [];
    cellInfo = [];    
elseif cellBiasFlag == 1
    %% M3 is used for gaussian-based cell anisotropy, it is not based on
    %% patient Light microscopy study.  
    [xTemp,yTemp,zTemp,cellIndx,cellFillCounter] = placeSpheresM3(cellSigma,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim);
    R_select_orig = [];
    R_select = [];
    cellSpheres = [];
    cellInfo = [];    
elseif cellBiasFlag == 2    
    %% M5 considers inter-cellular iron anisotropy based on patient Light
    %% microscopy analysis. Within cells, iron is distributed in random
    %% fashion.
    [xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select] = placeSpheresM5(patientInfo,patientIndx,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim);
    cellSpheres = [];
    cellInfo = [];
elseif cellBiasFlag == 3
    clear xTemp yTemp zTemp;       % used in a different way
    
    % M6 accounts for inter-cellular iron anisotropy (from LM) as well as
    % inter-paticle distances (nearest neighbor from EM) to distribute
    % iron.
    %tic
    [cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM6(NNfactor,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
    %toc
    % cellSpheres
    
%     %% in M8, modification of M6, with a 10um diameter cylindrical sinusoidal region created in the center of
%     %% the simulation box. 1/3rd of the iron goes into the sinusoid.
%     %tic
%     [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM8(NNfactor,NNfactorS,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
%     %toc
%     sphereInfo.cellSpheresSinusoids = cellSpheresSinusoids;

    
%     %% second version of M6 which uses a function to distribute spheres in
%     %% each cell, incase it has to be employed over the cluster
%     tic
%     [cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select] = placeSpheresM6_2(patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
%     toc

    %% in M7, new sphere is placed with respect to random sphere already
    %% placed, hence not forming a long chain but a circular blob, extreme
    %% case, not used in simulation
    %tic
    %%[cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM7(NNfactor,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
    %toc
    
elseif cellBiasFlag == 4
    
%     %% in M9, modification of M8, simulates 18 small cylindrical sinusoidal
%     %% regions with a 10um diameter in the simulation box. the are arranged in a non overlapping manner 
%     %% 1/3rd of the iron goes into the sinusoids equally.
%     %tic
%     [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM9(NNfactor,NNfactorS,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
%     %toc
%     sphereInfo.cellSpheresSinusoids = cellSpheresSinusoids;

    
%     %% M9P is the cluster version of M9. Sinusoidal iron deposition is
%     %% performed on the cluster.
%     %tic
%     [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM9P(NNfactor,NNfactorS,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
%     %toc
%     sphereInfo.cellSpheresSinusoids = cellSpheresSinusoids;
    

    %% In M9P2, both sinusoidal and hepatocyte iron deposition are
    %% performed over the cluster
    %tic
    [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM9P2(NNfactor,NNfactorS,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
    %toc
    sphereInfo.cellSpheresSinusoids = cellSpheresSinusoids;        
 
elseif cellBiasFlag == 5
    
    %% In M10P2, both sinusoidal and hepatocyte iron deposition are
    %% performed over the cluster
    tic
    [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM10P2(NNfactor,NNfactorS,patientInfo,patientIndx,[],[],[],sim_box_side,R_select,hepatocyte_dim);
    toc
    sphereInfo.cellSpheresSinusoids = cellSpheresSinusoids;        
    
    
end


numSpheres = length(R_select);

%% 3D visualization of sphere positions
% SphereVisualize(R_select,xTemp,yTemp,zTemp,sim_box_side);

%% Store in a structure, radius and center coordinate

%% for all
sphereInfo.x = xTemp;
sphereInfo.y = yTemp;
sphereInfo.z = zTemp;

%% for M3, M5, M6
sphereInfo.cellIndx = cellIndx;
sphereInfo.cellFillCounter = cellFillCounter;
sphereInfo.cellSigma = cellSigma;

%% for M5, M6
sphereInfo.radiusOriginal = R_select_orig;
sphereInfo.radius = R_select;

%% for M6
sphereInfo.cellSpheres = cellSpheres;
sphereInfo.cellInfo = cellInfo;



% % Visualize sphere positions
% lims = [-sim_box_side/2 sim_box_side/2];
% figure; hold on;
% grid; 
% xlim(lims), ylim(lims), zlim(lims)
% plot3(sphereInfo.xTemp,yTemp,zTemp,'.','markersize',6);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% set(gca, 'ZTick', []);
% set(gca, 'XTickLabel', []);
% set(gca, 'YTickLabel', []);
% set(gca, 'ZTickLabel', []);
% box on;


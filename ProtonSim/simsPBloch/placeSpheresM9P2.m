function [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM9P2(NNfactor,NNfactorS,patientInfo,patientIndx,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim)
%% Spheres are placed in a box in a non-overlapping manner.
%% If overlap occurs, the sphere is placed at the surface of colliding
%% sphere and then checked again for overlap. The procedure is continued
%% for successive overlaps and if there is a deadlock, a new position is
%% generated and the procedure is repeated until no overlap.
%% Variant of placeSpheresM2.
%% Here a gaussian function is used to decide the probablity of a sphere
%% being placed in a hepatocyte (cell).
%% sim_box_side should be a multiple of (cubiod) cell dimensions.
%% 27 Oct. 2006
%% Nilesh Ghugre, CHLA/USC
%% 22 Dec. 2006, Optimized and modified with following changes
%% - sphere radii are sorted in descending order so that largest spheres
%% are placed first. As spheres get filled, it is easier to find a place
%% for a smaller sphere than for a large one.
%% - findCellSpaceCount > 50, try 50 times before a sphere is thrown into
%% an undesignated cell
%% - sphere overlap is checked only with spheres already present in the
%% same cell.
%% 2 Jan 2007, Corrected for spheres getting placed near the center of the
%% environment (0,0,0) since compensation for coordinate conversion was not
%% correctly done when sphere collisions occured.
%% 11 Jan 2007, Detected that when cell reshuffling is being done, cell
%% numbers were getting repeated. Corrected this to generate unique cell
%% id's in random order.
%% 12 Jan 2007, Used 'randperm' to generate randomly ordered cell numbers
%% 26 Feb. 2007, Changed cell iron distribution to be exponential similar
%% to that observed in patients using Light microscopy
%% 23 March, 2007, Cell Iron distribution is based on cellIronPDF from
%% individual patients. An iron-fraction is obtained for each cell which
%% when applied to total number of spheres (or total iron content) gives
%% the number of spheres that need to be placed in a given cell.
%% 11 April, 2007, Completely modified version of placeSpheresM5. Here
%% inter-cellular iron distribution and inter-paticle distances have been
%% accounted for.
%% 3 May, 2007, here new sphere is placed with respect to last sphere
%% placed.
%% 10 May, 2007, creating a separate cylindrical compartment for sinusoids
%% right in the center of the simulation volume. The diameter of the
%% cylinder is 10 um and it hold ~30% of the iron (spheres)
%% 18 July, 2007, modification of M8, 18 small cylindrical sinusoidal
%% regions are placed in a non-overlapping manner.
%% 18 July, 2007, Implementation of M9 over the cluster for sinusoidal iron
%% deposition
%% 19 July, 2007, Hepatocyte iron deposition is implemented over the
%% cluster here.
 
%% random initialization of rand
rand('twister',sum(100*clock));

cylinderDiameter = 10;          %% um, book reference from Ignacio
sinusoidIronFraction = patientInfo(patientIndx).sinusoidIronFraction; % 0.3;     %% Hermatz et. al. Blood, 1 July 2000, vol 96, no. 1, Figure 2


totalNumSpheres = length(R_select);

nnDist = patientInfo(patientIndx).nnDist;
nnCDF = patientInfo(patientIndx).nnCDF;
nnPDF = patientInfo(patientIndx).nnPDF;

%%%%%%% centers of cells define their boundary

numNucleiInX = sim_box_side/hepatocyte_dim; % in single dim

[x,y,z] = meshgrid(-sim_box_side/2:sim_box_side/2);

[centerX,centerY,centerZ] = meshgrid(-(sim_box_side/2)+(sim_box_side/numNucleiInX)/2:sim_box_side/numNucleiInX:sim_box_side/2);

centerX = centerX(:);
centerY = centerY(:);
centerZ = centerZ(:);


numSpheresSinusoid = floor(totalNumSpheres*sinusoidIronFraction);
numSpheresHepatocyte = totalNumSpheres - numSpheresSinusoid;

R_select_orig = R_select;
clear R_select;

R_select_S = R_select_orig(1:numSpheresSinusoid);


%%%%%% Organizing data for sinusoids %%%%%%%%%%%%%

%% sinusoids comprise of ~6% of the volume of the liver lobe.
% totalSinusoidalVolume = sim_box_side^3 * 0.06;

%%% total cylindrical length of sinusoids
% totalSinusoidalLength = totalSinusoidalVolume/(pi*(cylinderDiameter/2)^2);

%% single sinusoid from M8 with diameter 10um and length 80um has a volume
%% of pi*5^2*80 = 6.2832e+003, which is 6.2832e+003/80^3 = 1.23% of
%% simulation volume. Now, 1.23*5 = 6.15 which is approximately equal to
%% sinusoidal volume in liver lobe.

%% So, we can have 5 cylindrical sinusoids running parallel at the
%% intersections of the hepatocytes.
%% However, this will be a different geometry and there will be the
%% infinite cylinder problem. Hence we will have sinusoids with length =
%% hepatocyte dimensions and placed in a non-continuous manner.
%% 18 sinusoids can be placed in such a configuration. In this case,
%% sinusoid volume % = pi*25*(18*20)/80^3 = 5.52% (~6%)

numberOfSinusoids = 18;

numSpheresSinusoidEach(1:numberOfSinusoids) = floor(numSpheresSinusoid/numberOfSinusoids);

%% fill in excess spheres in the first sinusoid.
numSpheresSinusoidEach(1) = numSpheresSinusoidEach(1) + (numSpheresSinusoid - sum(numSpheresSinusoidEach));


R_select(1).r = R_select_S(1:numSpheresSinusoidEach(1));
xTemp(1).x = zeros(1,numSpheresSinusoidEach(1));
yTemp(1).y = zeros(1,numSpheresSinusoidEach(1));
zTemp(1).z = zeros(1,numSpheresSinusoidEach(1));
for k = 2:numberOfSinusoids

    R_select(k).r = R_select_S(1+sum(numSpheresSinusoidEach(1:k-1)):sum(numSpheresSinusoidEach(1:k)));
    %figure;plot(R_select(k).r-R_select_S(1+sum(numSpheresSinusoidEach(1:k-1)):sum(numSpheresSinusoidEach(1:k))));

    xTemp(k).x = zeros(1,numSpheresSinusoidEach(k));
    yTemp(k).y = zeros(1,numSpheresSinusoidEach(k));
    zTemp(k).z = zeros(1,numSpheresSinusoidEach(k));

end

%%% testing if correct
% R_select_temp = [];
% for k = 1:numberOfSinusoids
%     R_select_temp = [ R_select_temp R_select(k).r];
% end
% figure;plot(R_select_temp-R_select_S);

%% offset coordinates for the sinusoids from the center of the simulation
%% volume.

%% level 1

xTemp(1).Offset = 0;
yTemp(1).Offset = 0;

xTemp(2).Offset = -hepatocyte_dim;
yTemp(2).Offset = hepatocyte_dim;

xTemp(3).Offset = hepatocyte_dim;
yTemp(3).Offset = hepatocyte_dim;

xTemp(4).Offset = hepatocyte_dim;
yTemp(4).Offset = -hepatocyte_dim;

xTemp(5).Offset = -hepatocyte_dim;
yTemp(5).Offset = -hepatocyte_dim;

%% z limits are same for all in the same level
for k = 1:5
    zTemp(k).limits = [hepatocyte_dim 2*hepatocyte_dim];
end

%% level 3, same as level 1

xTemp(6).Offset = 0;
yTemp(6).Offset = 0;

xTemp(7).Offset = -hepatocyte_dim;
yTemp(7).Offset = hepatocyte_dim;

xTemp(8).Offset = hepatocyte_dim;
yTemp(8).Offset = hepatocyte_dim;

xTemp(9).Offset = hepatocyte_dim;
yTemp(9).Offset = -hepatocyte_dim;

xTemp(10).Offset = -hepatocyte_dim;
yTemp(10).Offset = -hepatocyte_dim;

%% z limits are same for all in the same level
for k = 6:10
    zTemp(k).limits = [0 -hepatocyte_dim];
end

%% level 2

xTemp(11).Offset = 0;
yTemp(11).Offset = -hepatocyte_dim;

xTemp(12).Offset = 0;
yTemp(12).Offset = hepatocyte_dim;

xTemp(13).Offset = hepatocyte_dim;
yTemp(13).Offset = 0;

xTemp(14).Offset = -hepatocyte_dim;
yTemp(14).Offset = 0;

%% z limits are same for all in the same level
for k = 11:14
    zTemp(k).limits = [0 hepatocyte_dim];
end

%% level 4, same as level 2

xTemp(15).Offset = 0;
yTemp(15).Offset = -hepatocyte_dim;

xTemp(16).Offset = 0;
yTemp(16).Offset = hepatocyte_dim;

xTemp(17).Offset = hepatocyte_dim;
yTemp(17).Offset = 0;

xTemp(18).Offset = -hepatocyte_dim;
yTemp(18).Offset = 0;

%% z limits are same for all in the same level
for k = 15:18
    zTemp(k).limits = [-hepatocyte_dim -hepatocyte_dim*2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% testing function
% placeSinusoidalSpheresM9P(1,numSpheresSinusoidEach,xTemp,yTemp,zTemp,R_select,cylinderDiameter,nnCDF,nnDist,NNfactorS)

%%%%% Send to cluster %%%%%%%%%%
jm = findResource('scheduler','type','jobmanager','name','default_jobmanager','LookupURL','localhost');
job = createJob(jm);
set(job, 'MaximumNumberOfWorkers', numberOfSinusoids);
set(job,'FileDependencies',{'placeSinusoidalSpheresM9P.m'});

for sinusoidCount = 1:numberOfSinusoids

    createTask(job,@placeSinusoidalSpheresM9P,3,{sinusoidCount,numSpheresSinusoidEach,xTemp,yTemp,zTemp,R_select,cylinderDiameter,nnCDF,nnDist,NNfactorS});

end


%tic
submit(job)
waitForState(job)

results = getAllOutputArguments(job);
%toc

% get(get(job,'Tasks'),'ErrorMessage')

destroy(job);



xTempS = [];
yTempS = [];
zTempS = [];
for k = 1:18
    
    myXTemp = results{k,1};
    myYTemp = results{k,2};
    myZTemp = results{k,3};
    
    xTempS = [xTempS myXTemp(k).x];
    yTempS = [yTempS myYTemp(k).y];
    zTempS = [zTempS myZTemp(k).z];
    
    clear myXTemp myYTemp myZTemp

end

% lims = [-sim_box_side/2 sim_box_side/2];
% figure; hold on;
% grid; 
% xlim(lims), ylim(lims), zlim(lims)
% plot3(xTempS,yTempS,zTempS,'.','markersize',4);
% set(gca, 'XTick', [-20 0 20]);
% set(gca, 'YTick', [-20 0 20]);
% set(gca, 'ZTick', [-20 0 20]);
% set(gca, 'XTickLabel', []);
% set(gca, 'YTickLabel', []);
% set(gca, 'ZTickLabel', []);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear results
clear xTemp yTemp zTemp;

%%%%%%% Finding the cells in which the sinusoidal spheres lie in.

sinusoidSphereCell = zeros(1,numSpheresSinusoid);

for sphereCount = 1:numSpheresSinusoid

    for p = 1:length(centerX)

        %% cell bounds

        aX = (centerX(p)-hepatocyte_dim/2);
        bX = (centerX(p)+hepatocyte_dim/2);

        aY = (centerY(p)-hepatocyte_dim/2);
        bY = (centerY(p)+hepatocyte_dim/2);

        aZ = (centerZ(p)-hepatocyte_dim/2);
        bZ = (centerZ(p)+hepatocyte_dim/2);

        %% position has to be checked if sphere is going outside cell
        %% boundary
        nextPos.x = xTempS(sphereCount);
        nextPos.y = yTempS(sphereCount);
        nextPos.z = zTempS(sphereCount);
        [nextPos,environCrossFlag] = environCross3(nextPos,aX,bX,aY,bY,aZ,bZ);
        
        if environCrossFlag == 0
            sinusoidSphereCell(sphereCount) = p;
            break;
        end

    end

end

for p = 1:length(centerX)

    cellSpheresSinusoids(p) = sum(sinusoidSphereCell == p);

end
             
disp('Finished placing sinusoidal spheres...');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hepatocyte deposits

%%%%%%%%%% cell spheres distribution

%% Get from PatientInfo
normFE = patientInfo(patientIndx).normFE;
cellIronCDF = patientInfo(patientIndx).cellIronCDF;
% cellIronPDF = patientInfo(patientIndx).cellIronPDF;

% figure;plot(normFE,cellIronCDF);

%% first randomize cell selection
cellNumRand = randperm(length(centerX));
% figure;plot(cellNumRand);

cellIronCoeff = zeros(1,length(centerX));
%% find iron content in each cell based on distribution
for p = 1:length(centerX)

    rand_num = rand(1,1);
    I = find(cellIronCDF>rand_num);
    cellIronCoeff(cellNumRand(p)) = normFE(I(1));

end

% figure;bar(totalNumSpheres*cellIronCoeff./sum(cellIronCoeff));
cellSpheres = floor(numSpheresHepatocyte*cellIronCoeff./sum(cellIronCoeff));

%% Number of spheres change slightly because of rounding...
numSpheresHepatocyte = sum(cellSpheres);

R_select_H = R_select_orig(numSpheresSinusoid+1:numSpheresSinusoid+numSpheresHepatocyte);

%% Rearranging sphere radius vector for each cell
for k = 1:length(cellSpheres)

    if k==1
        offset = 0;
    else
        offset = sum(cellSpheres(k-1));
    end

    cellSpheresRadius(k).radiusVec = R_select_H((offset+1):(offset+cellSpheres(k)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% testing function
% placeHepatocyteSpheresM9P(1,cellSpheres,cellSpheresRadius,centerX,centerY,centerZ,hepatocyte_dim,sinusoidSphereCell,xTempS,yTempS,zTempS,R_select_S,nnCDF,nnDist,NNfactor)

%%%%% Send to cluster %%%%%%%%%%
jm = findResource('scheduler','type','jobmanager','name',getJobmanagerInfo('jm'),'LookupURL',getJobmanagerInfo('host'));
job = createJob(jm);
set(job, 'MaximumNumberOfWorkers', length(centerX)/2);
set(job,'FileDependencies',{'placeHepatocyteSpheresM9P.m' 'environCross3.m'});

%% consider each of the 64 cells for sphere placement
for p = 1:length(centerX)

    createTask(job,@placeHepatocyteSpheresM9P,1,{p,cellSpheres,cellSpheresRadius,centerX,centerY,centerZ,hepatocyte_dim,sinusoidSphereCell,xTempS,yTempS,zTempS,R_select_S,nnCDF,nnDist,NNfactor});

end


%tic
submit(job)
waitForState(job)

results = getAllOutputArguments(job);
%toc

% get(get(job,'Tasks'),'ErrorMessage')

destroy(job);


for p = 1:length(centerX)
    
    tempVar = results{p};
    cellInfo(p).xTemp = tempVar(p).xTemp;
    cellInfo(p).yTemp = tempVar(p).yTemp;
    cellInfo(p).zTemp = tempVar(p).zTemp;


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear results
disp('Finished placing hepatocyte spheres...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Returning variables in the required format.

R_select = [R_select_S R_select_H];
cellFillCounter = 0;    %% not used  
cellIndx = [];          %% not used


xTempH = [];
yTempH = [];
zTempH = [];

for p = 1:length(cellInfo);
    
   xTempH = [xTempH cellInfo(p).xTemp];
   yTempH = [yTempH cellInfo(p).yTemp];
   zTempH = [zTempH cellInfo(p).zTemp]; 
    
end

clear xTemp yTemp zTemp

xTemp = [xTempS xTempH];
yTemp = [yTempS yTempH];
zTemp = [zTempS zTempH];












function [cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select] = placeSpheresM6P(patientInfo,patientIndx,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim);
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
%% 27 April, 2007, Modified to be evaluated over the cluster

%% random initialization of rand
rand('twister',sum(100*clock));

% %%% nearest neighbor information for patients
% load patientNN;
% patientSelect = [2:12 14];
% patientNN = patientNN(patientSelect);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalNumSpheres = length(R_select);


% nnDist = patientNN(patientIndx).nnDist;
% nnCDF = patientNN(patientIndx).nnCDF;
% nnPDF = patientNN(patientIndx).nnPDF;

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
cellSpheres = floor(totalNumSpheres*cellIronCoeff./sum(cellIronCoeff));

%% Number of spheres change slightly because of rounding...
totalNumSpheres = sum(cellSpheres);

R_select_orig = R_select;
R_select = R_select(1:totalNumSpheres);
R_select_new = R_select;    %% for the purpose of returning

%% Rearranging sphere radius vector for each cell
for k = 1:length(cellSpheres)

    if k==1
        offset = 0;
    else
        offset = sum(cellSpheres(k-1));
    end

    cellSpheresRadius(k).radiusVec = R_select((offset+1):(offset+cellSpheres(k)));

end


clear R_select;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for whichCell = 1:length(cellSpheres)
    whichCell;
    [cellInfo(whichCell)] = getCellInfo(whichCell,cellSpheres,cellSpheresRadius,centerX,centerY,centerZ,sim_box_side,hepatocyte_dim,nnCDF,nnDist);
end


%%%% Returning variables in the required format.

R_select = R_select_new;
cellFillCounter = 0;    %% not used  
cellIndx = [];          %% not used


xTemp = [];
yTemp = [];
zTemp = [];


for p = 1:length(cellInfo);
    
   xTemp = [xTemp cellInfo(p).xTemp];
   yTemp = [yTemp cellInfo(p).yTemp];
   zTemp = [zTemp cellInfo(p).zTemp]; 
    
end





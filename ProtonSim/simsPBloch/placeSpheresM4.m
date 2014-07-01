function [xTemp,yTemp,zTemp,cellIndx,cellFillCounter] = placeSpheresM4(patientInfo,patientIndx,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim);
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

%% random initialization of rand
rand('twister',sum(100*clock));

numSpheres = length(xTemp);

%%%%%%% centers of cells define their boundary

numNucleiInX = sim_box_side/hepatocyte_dim; % in single dim

[x,y,z] = meshgrid(-sim_box_side/2:sim_box_side/2);

[centerX,centerY,centerZ] = meshgrid(-(sim_box_side/2)+(sim_box_side/numNucleiInX)/2:sim_box_side/numNucleiInX:sim_box_side/2);

centerX = centerX(:);
centerY = centerY(:);
centerZ = centerZ(:);

%% Get inter-cell iron distribution from patientInfo.

% %% exponential for testing
% normFE = linspace(0,1,1000);
% cellIronPDF = exp(-normFE*5);
% % figure;plot(normFE,cellIronPDF);
% cellIronCDF = cumsum(cellIronPDF);

%% Get from PatientInfo
normFE = patientInfo(patientIndx).normFE;
cellIronCDF = patientInfo(patientIndx).cellIronCDF;

cellIronCDF = length(centerX)* cellIronCDF./max(cellIronCDF);    % normalize to number of cells
% figure;plot(normFE,cellIronCDF);

%% first randomize cell selection
cellNumRand = randperm(length(centerX));
% figure;plot(cellNumRand);

cellIndx = zeros(1,numSpheres);

%% place the first sphere randomly in a chosen cell
%% choose a cell based on distribution
p = 1;

rand_num = min(normFE) + (max(normFE)-min(normFE))*rand(1,1);
cellNum = cellIronCDF(length(find(normFE<rand_num)));
cellIndx(p) = cellNumRand(ceil(cellNum));

a = (centerX(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
b =  (centerX(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
xTemp(p) = a + (b-a) * rand(1);

a = (centerY(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
b =  (centerY(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
yTemp(p) = a + (b-a) * rand(1);

a = (centerZ(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
b =  (centerZ(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
zTemp(p) = a + (b-a) * rand(1);


%% Counter to note number of times a different cell is chosen due to space
%% filled.
cellFillCounter = 0;


for p=2:numSpheres


    overlapFlag = 1;

    %% choose a cell based on distribution
    rand_num = min(normFE) + (max(normFE)-min(normFE))*rand(1,1);
    cellNum = cellIronCDF(length(find(normFE<rand_num)));
    cellIndx(p) = cellNumRand(ceil(cellNum));

    findCellSpaceCount = 0;     %% count for tracking if cell has space for a new sphere

    while overlapFlag == 1

        %% if new sphere cannot find space within X trys, a new cell is
        %% selected at random
        if findCellSpaceCount > 50
            randCell = round(1+ (length(centerX)-1)*rand(1,1));
            cellIndx(p) = randCell;
            findCellSpaceCount = 0;
            cellFillCounter = cellFillCounter + 1;
        end

        %% new box coordinates such that sphere does not lie outside box dims
        a = (centerX(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
        b =  (centerX(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
        xTemp(p) = a + (b-a) * rand(1);

        a = (centerY(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
        b =  (centerY(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
        yTemp(p) = a + (b-a) * rand(1);

        a = (centerZ(cellIndx(p))-hepatocyte_dim/2)+R_select(p);
        b =  (centerZ(cellIndx(p))+hepatocyte_dim/2)-R_select(p);
        zTemp(p) = a + (b-a) * rand(1);

        numPlaceOnSurface = 0;

        %% Check for overlap with previously placed spheres in the
        %% currently chosen cell
        currCellSpheres = find(cellIndx == cellIndx(p));
        currCellSpheres = currCellSpheres(currCellSpheres < p); % so that sphere is not compared with itself


        while numPlaceOnSurface < 10


            if isempty(currCellSpheres)
                overlapFlag = 0;
            else

                %% Check for overlap with previously placed spheres
                for j = 1:length(currCellSpheres)

                    centerDist = sqrt( (xTemp(p)-xTemp(currCellSpheres(j))).^2 + (yTemp(p)-yTemp(currCellSpheres(j))).^2 + (zTemp(p)-zTemp(currCellSpheres(j))).^2 );

                    if centerDist > (R_select(p)+R_select(currCellSpheres(j)))
                        overlapFlag = 0;
                    else
                        overlapFlag = 1;
                        insideWhichSphere = currCellSpheres(j);
                        break;
                    end
                end     % j, overlap check count

            end


            if overlapFlag == 1
                %% convert position of current sphere to polar wrt to collided
                %% sphere, r_1 will be distance between their centers
                [theta_1,phi_1,r_1]=cart2sph(xTemp(p)-xTemp(insideWhichSphere),yTemp(p)-yTemp(insideWhichSphere),zTemp(p)-zTemp(insideWhichSphere));

                %% place current sphere on the surface of collided sphere
                r_1 = R_select(p)+R_select(insideWhichSphere);

                %% reconvert
                [xTemp(p),yTemp(p),zTemp(p)]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
                %% getting back into original coordinate system
                xTemp(p) = xTemp(p)+xTemp(insideWhichSphere);
                yTemp(p) = yTemp(p)+yTemp(insideWhichSphere);
                zTemp(p) = zTemp(p)+zTemp(insideWhichSphere);

                numPlaceOnSurface = numPlaceOnSurface + 1;

            else
                break;
            end

        end     % while numPlaceOnSurface

        findCellSpaceCount = findCellSpaceCount + 1;

    end     % while overlap
end     % for p

cellFillCounter

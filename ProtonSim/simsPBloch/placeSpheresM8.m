function [cellSpheresSinusoids,cellSpheres,cellInfo,xTemp,yTemp,zTemp,cellIndx,cellFillCounter,R_select_orig,R_select,nnDist,nnPDF,nnCDF] = placeSpheresM8(NNfactor,NNfactorS,patientInfo,patientIndx,xTemp,yTemp,zTemp,sim_box_side,R_select,hepatocyte_dim)
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

cylinderDiameter = 10;          %% um, book from Ignacio
sinusoidIronFraction = 0.3;     %% Hermatz et. al., Blood, Fig.2
%% sinusoidal volume is 6% of liver lobe, The Clinical Electron Microscopy
%% Society of Japan 2004

%% random initialization of rand
rand('twister',sum(100*clock));

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
R_select = R_select_S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Placing Sinusoidal deposits

for sphereCount = 1:numSpheresSinusoid

    %% placing the first sphere
    if sphereCount == 1
        %% place the sphere anywhere in the 10um diameter cylindrical
        %% region with long axis passing thru origin and parallel to z
        %% axis.
        a = -sim_box_side/2+R_select(sphereCount);
        b =  sim_box_side/2-R_select(sphereCount);
        zTemp(sphereCount) = a + (b-a) * rand(1);

        axialR = 0 + (cylinderDiameter/2-0)*rand(1);
        axialTheta = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians
        xTemp(sphereCount) = axialR*cos(axialTheta);
        yTemp(sphereCount) = axialR*sin(axialTheta);

    elseif ~mod(sphereCount,100000)         %% for forming groups of spheres, not used for now

        overlapFlag = 1;
        while overlapFlag == 1
            a = -sim_box_side/2+R_select(sphereCount);
            b =  sim_box_side/2-R_select(sphereCount);
            zTemp(sphereCount) = a + (b-a) * rand(1);

            axialR = 0 + (cylinderDiameter/2-0)*rand(1);
            axialTheta = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians
            xTemp(sphereCount) = axialR*cos(axialTheta);
            yTemp(sphereCount) = axialR*sin(axialTheta);


            %% Check for overlap with previously placed spheres
            for j = 1:sphereCount-1

                centerDist = sqrt( (xTemp(sphereCount)-xTemp(j)).^2 + (yTemp(sphereCount)-yTemp(j)).^2 + (zTemp(sphereCount)-zTemp(j)).^2 );

                if centerDist > (R_select(sphereCount)+R_select(j))
                    overlapFlag = 0;
                else
                    overlapFlag = 1;
                    insideWhichSphere = j;
                    break;
                end

            end     % j, overlap check count

        end     %% while overlap

    else

        %% choose any sphere already placed.
        % selectedSphIndx = round(1 + ((sphereCount - 1)-1) * rand(1));

        %% choose previous sphere placed
        selectedSphIndx = sphereCount - 1;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% place a sphere next to the selected sphere at a distance obtained
        %% from the patient distribution
        %% this distance has to be greater than the sphere
        %% radius,else collision has occured !
        %% hence the loop...
        currDist = 0;
        %while currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
        rand_num = rand(1,1);
        I = find(nnCDF>rand_num);
        currDist = nnDist(I(1)) * NNfactorS;    % distance

        if currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
            currDist = (R_select(selectedSphIndx)+R_select(sphereCount))*(1+0.05);   %% 5% more than sticking distance
        end
        %end

        overlapFlag = 1;
        environCrossFlag = 1;

        repeatDistanceCount = 1;
        %% counter for the use of currDist. Use it for a specified
        %% number of times then regenerate. Doing this will prevent
        %% a collision deadlock if the number of spheres is large.

        newSeedCount = 1;
        %% if new position is not found near the selected sphere,
        %% change the selected sphere to a random sphere already
        %% placed.

        while overlapFlag == 1 || environCrossFlag == 1

            %z = z+1;
            if repeatDistanceCount > 10
                %disp('in repeatDistanceCount')
                if newSeedCount > 4
                    %disp('in newSeedCount')
                    %% randomly select new starting sphere which
                    %% has been already placed.
                    selectedSphIndx = round(1 + ((sphereCount - 1)-1) * rand(1));

                    newSeedCount = 1;
                    %repeatDistanceCount = 1;
                end

                currDist = 0;
                % while currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
                rand_num = rand(1,1);
                I = find(nnCDF>rand_num);
                currDist = nnDist(I(1)) * NNfactorS;    % distance

                if currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
                    currDist = (R_select(selectedSphIndx)+R_select(sphereCount))*(1+0.05);   %% 5% more than sticking distance
                end

                % end

                overlapFlag = 1;
                environCrossFlag = 1;
                repeatDistanceCount = 1;
                newSeedCount = newSeedCount + 1;

            end

            %% for the next sphere to be placed, r_2 = currDist and choose random theta and phi.
            %% this position is with respect to origin(0,0,0).
            r_1 = currDist;
            theta_1 = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians
            phi_1 = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians

            %% convert to cartesian
            [displacementX,displacementY,displacementZ]=sph2cart(theta_1,phi_1,r_1);

            %% to get real position , use addition of vectors.
            xTemp(sphereCount) = xTemp(selectedSphIndx) + displacementX;
            yTemp(sphereCount) = yTemp(selectedSphIndx) + displacementY;
            zTemp(sphereCount) = zTemp(selectedSphIndx) + displacementZ;

            %% position has to be checked if sphere is going outside
            %% cylindrical sinusoidal boundary

            if ((sqrt(xTemp(sphereCount)^2+yTemp(sphereCount)^2)) > (cylinderDiameter/2)) || (zTemp(sphereCount)<(-sim_box_side/2)) || (zTemp(sphereCount)>(sim_box_side/2))
                environCrossFlag = 1;
            else
                environCrossFlag = 0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%

            if environCrossFlag == 0

                numPlaceOnSurface = 0;
                while numPlaceOnSurface < 10

                    %% Check for overlap with previously placed spheres
                    for j = 1:sphereCount-1

                        centerDist = sqrt( (xTemp(sphereCount)-xTemp(j)).^2 + (yTemp(sphereCount)-yTemp(j)).^2 + (zTemp(sphereCount)-zTemp(j)).^2 );

                        if centerDist > (R_select(sphereCount)+R_select(j))
                            overlapFlag = 0;
                        else
                            overlapFlag = 1;
                            insideWhichSphere = j;
                            break;
                        end

                    end     % j, overlap check count

                    if overlapFlag == 1
                        %% convert position of current sphere to polar wrt to collided
                        %% sphere, r_1 will be distance between their centers
                        [theta_1,phi_1,r_1]=cart2sph(xTemp(sphereCount)-xTemp(insideWhichSphere),yTemp(sphereCount)-yTemp(insideWhichSphere),zTemp(sphereCount)-zTemp(insideWhichSphere));

                        %% place current sphere on the surface of collided sphere
                        r_1 = R_select(sphereCount)+R_select(insideWhichSphere);

                        %% reconvert
                        [xTemp(sphereCount),yTemp(sphereCount),zTemp(sphereCount)]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
                        %% getting back into original coordinate system
                        xTemp(sphereCount) = xTemp(sphereCount)+xTemp(insideWhichSphere);
                        yTemp(sphereCount) = yTemp(sphereCount)+yTemp(insideWhichSphere);
                        zTemp(sphereCount) = zTemp(sphereCount)+zTemp(insideWhichSphere);

                        numPlaceOnSurface = numPlaceOnSurface + 1;

                    else
                        break;
                    end


                end     %% while numPlaceOnSurface


            end     %% if environCrossFlag

            repeatDistanceCount = repeatDistanceCount + 1;

        end     %% while overlapFlag

    end     % if sphereCount


%     if ~mod(sphereCount,500)
%         sphereCount
%     end


end     % for sphereCount


xTempS = xTemp;
yTempS = yTemp;
zTempS = zTemp;

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
%% consider each of the 64 cells for sphere placement
for p = 1:length(centerX)

    p;
    numSpheres = cellSpheres(p);

    R_select = cellSpheresRadius(p).radiusVec;

    if numSpheres ~= 0

        %% cell bounds

        aX = (centerX(p)-hepatocyte_dim/2);
        bX = (centerX(p)+hepatocyte_dim/2);

        aY = (centerY(p)-hepatocyte_dim/2);
        bY = (centerY(p)+hepatocyte_dim/2);

        aZ = (centerZ(p)-hepatocyte_dim/2);
        bZ = (centerZ(p)+hepatocyte_dim/2);


        %%%%%% Method 2

        clear xTemp yTemp zTemp;

        groupSize = 10000;     % number of spheres in a given group, started with 50
        %% for a very large number, no groups will be created, particles
        %% will be placed in a long chain. The loop for remaining spheres
        %% will only be executed and not the 'groupIndx' for loop.

        numGroups = floor(numSpheres/groupSize);
        remainingSpheres = rem(numSpheres,groupSize);

        sphereCount = 1;

        for groupIndx = 1:numGroups

            %% removed contents since for very large groupSize, this loop
            %% will not be executed
                
        end     %% for groupIndx

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% INSERT LAST GROUP FOR REMAINING SPHERES...

        if remainingSpheres ~= 0

            groupCount = 1;

            overlapFlag = 1;
            overlapFlagS = 1;

            while overlapFlag == 1 || overlapFlagS == 1

                %% place the first group sphere at random in the cell.
                %% new box coordinates such that sphere does not lie outside
                %% cell dims
                a = (centerX(p)-hepatocyte_dim/2)+R_select(sphereCount);
                b =  (centerX(p)+hepatocyte_dim/2)-R_select(sphereCount);
                xTemp(sphereCount) = a + (b-a) * rand(1);

                a = (centerY(p)-hepatocyte_dim/2)+R_select(sphereCount);
                b =  (centerY(p)+hepatocyte_dim/2)-R_select(sphereCount);
                yTemp(sphereCount) = a + (b-a) * rand(1);

                a = (centerZ(p)-hepatocyte_dim/2)+R_select(sphereCount);
                b =  (centerZ(p)+hepatocyte_dim/2)-R_select(sphereCount);
                zTemp(sphereCount) = a + (b-a) * rand(1);

                
                %%%%%%%% Overlap with hepatocyte spheres
                if sphereCount > 1
                    %% Check for overlap with previously placed spheres
                    for j = 1:sphereCount-1

                        centerDist = sqrt( (xTemp(sphereCount)-xTemp(j)).^2 + (yTemp(sphereCount)-yTemp(j)).^2 + (zTemp(sphereCount)-zTemp(j)).^2 );

                        if centerDist > (R_select(sphereCount)+R_select(j))
                            overlapFlag = 0;
                        else
                            overlapFlag = 1;
                            insideWhichSphere = j;
                            break;
                        end

                    end
                else        %% first sphere, no need for overlap check
                    overlapFlag = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%% Overlap with sinusoid spheres

                selectedSpheres = find(sinusoidSphereCell == p);

                if ~isempty(selectedSpheres)

                    for j = 1:length(selectedSpheres)

                        centerDist = sqrt( (xTemp(sphereCount)-xTempS(selectedSpheres(j))).^2 + (yTemp(sphereCount)-yTempS(selectedSpheres(j))).^2 + (zTemp(sphereCount)-zTempS(selectedSpheres(j))).^2 );

                        if centerDist > (R_select(sphereCount)+R_select_S(selectedSpheres(j)))
                            overlapFlagS = 0;
                        else
                            overlapFlagS = 1;
                            insideWhichSphereS = j;
                            break;
                        end

                    end
                else
                    overlapFlagS = 0;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end     % j, overlap check count

            sphereCount = sphereCount + 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            while groupCount < remainingSpheres

                selectedSphIndx = sphereCount - 1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %% place a sphere next to the selected sphere at a distance obtained
                %% from the patient distribution
                %% this distance has to be greater than the sphere
                %% radius,else collision has occured !
                %% hence the loop...
                currDist = 0;
                while currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
                    rand_num = rand(1,1);
                    I = find(nnCDF>rand_num);
                    currDist = nnDist(I(1)) * NNfactor;    % distance
                end

                overlapFlag = 1;
                overlapFlagS = 1;
                environCrossFlag = 1;

                repeatDistanceCount = 1;
                %% counter for the use of currDist. Use it for a specified
                %% number of times then regenerate. Doing this will prevent
                %% a collision deadlock if the number of spheres is large.

                newSeedCount = 1;
                %% if new position is not found near the selected sphere,
                %% change the selected sphere to a random sphere already
                %% placed.
                
                %z = 0;
                while overlapFlag == 1 || environCrossFlag == 1 || overlapFlagS == 1

                    %z = z+1;
                    if repeatDistanceCount > 50
                        %disp('in repeatDistanceCount')
                        if newSeedCount > 20
                            %disp('in newSeedCount')
                            %% randomly select new starting sphere which
                            %% has been already placed.
                            selectedSphIndx = round(1 + ((sphereCount - 1)-1) * rand(1));
                            
                            newSeedCount = 1;
                            %repeatDistanceCount = 1;
                        end                       
                        
                        currDist = 0;
                        while currDist < (R_select(selectedSphIndx)+R_select(sphereCount))
                            rand_num = rand(1,1);
                            I = find(nnCDF>rand_num);
                            currDist = nnDist(I(1)) * NNfactor;    % distance
                        end

                        overlapFlag = 1;
                        environCrossFlag = 1;
                        repeatDistanceCount = 1;
                        newSeedCount = newSeedCount + 1;

                    end


                    %% for the next sphere to be placed, r_2 = currDist and choose random theta and phi.
                    %% this position is with respect to origin(0,0,0).
                    r_1 = currDist;
                    theta_1 = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians
                    phi_1 = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians

                    %% convert to cartesian
                    [displacementX,displacementY,displacementZ]=sph2cart(theta_1,phi_1,r_1);

                    %% to get real position , use addition of vectors.
                    xTemp(sphereCount) = xTemp(selectedSphIndx) + displacementX;
                    yTemp(sphereCount) = yTemp(selectedSphIndx) + displacementY;
                    zTemp(sphereCount) = zTemp(selectedSphIndx) + displacementZ;

                    %% position has to be checked if sphere is going outside cell
                    %% boundary
                    nextPos.x = xTemp(sphereCount);
                    nextPos.y = yTemp(sphereCount);
                    nextPos.z = zTemp(sphereCount);
                    [nextPos,environCrossFlag] = environCross3(nextPos,aX,bX,aY,bY,aZ,bZ);
                    %%%%%%%%%%%%%%%%%%%%%%%

                    if environCrossFlag == 0

                        numPlaceOnSurface = 0;
                        while numPlaceOnSurface < 10


                            %%%%%%%% Overlap with sinusoid spheres
                            selectedSpheres = find(sinusoidSphereCell == p);

                            if ~isempty(selectedSpheres)

                                for j = 1:length(selectedSpheres)

                                    centerDist = sqrt( (xTemp(sphereCount)-xTempS(selectedSpheres(j))).^2 + (yTemp(sphereCount)-yTempS(selectedSpheres(j))).^2 + (zTemp(sphereCount)-zTempS(selectedSpheres(j))).^2 );

                                    if centerDist > (R_select(sphereCount)+R_select_S(selectedSpheres(j)))
                                        overlapFlagS = 0;
                                    else
                                        overlapFlagS = 1;
                                        insideWhichSphereS = j;
                                        break;
                                    end

                                end
                            else
                                overlapFlagS = 0;
                            end

                            if overlapFlagS == 1
                                %% convert position of current sphere to polar wrt to collided
                                %% sphere, r_1 will be distance between their centers
                                [theta_1,phi_1,r_1]=cart2sph(xTemp(sphereCount)-xTempS(selectedSpheres(insideWhichSphereS)),yTemp(sphereCount)-yTempS(selectedSpheres(insideWhichSphereS)),zTemp(sphereCount)-zTempS(selectedSpheres(insideWhichSphereS)));

                                %% place current sphere on the surface of collided sphere
                                r_1 = R_select(sphereCount)+R_select_S(selectedSpheres(insideWhichSphereS));

                                %% reconvert
                                [xTemp(sphereCount),yTemp(sphereCount),zTemp(sphereCount)]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
                                %% getting back into original coordinate system
                                xTemp(sphereCount) = xTemp(sphereCount)+xTempS(selectedSpheres(insideWhichSphereS));
                                yTemp(sphereCount) = yTemp(sphereCount)+yTempS(selectedSpheres(insideWhichSphereS));
                                zTemp(sphereCount) = zTemp(sphereCount)+zTempS(selectedSpheres(insideWhichSphereS));

                                numPlaceOnSurface = numPlaceOnSurface + 1;

                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            
                            %% Check for overlap with previously placed spheres
                            for j = 1:sphereCount-1

                                centerDist = sqrt( (xTemp(sphereCount)-xTemp(j)).^2 + (yTemp(sphereCount)-yTemp(j)).^2 + (zTemp(sphereCount)-zTemp(j)).^2 );

                                if centerDist > (R_select(sphereCount)+R_select(j))
                                    overlapFlag = 0;
                                else
                                    overlapFlag = 1;
                                    insideWhichSphere = j;
                                    break;
                                end

                            end     % j, overlap check count

                            if overlapFlag == 1
                                %% convert position of current sphere to polar wrt to collided
                                %% sphere, r_1 will be distance between their centers
                                [theta_1,phi_1,r_1]=cart2sph(xTemp(sphereCount)-xTemp(insideWhichSphere),yTemp(sphereCount)-yTemp(insideWhichSphere),zTemp(sphereCount)-zTemp(insideWhichSphere));

                                %% place current sphere on the surface of collided sphere
                                r_1 = R_select(sphereCount)+R_select(insideWhichSphere);

                                %% reconvert
                                [xTemp(sphereCount),yTemp(sphereCount),zTemp(sphereCount)]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
                                %% getting back into original coordinate system
                                xTemp(sphereCount) = xTemp(sphereCount)+xTemp(insideWhichSphere);
                                yTemp(sphereCount) = yTemp(sphereCount)+yTemp(insideWhichSphere);
                                zTemp(sphereCount) = zTemp(sphereCount)+zTemp(insideWhichSphere);

                                numPlaceOnSurface = numPlaceOnSurface + 1;

                            else
                                break;
                            end


                        end     %% while numPlaceOnSurface


                    end     %% if environCrossFlag

                    repeatDistanceCount = repeatDistanceCount + 1;
                    
                end     %% while overlapFlag

                %z
                sphereCount = sphereCount + 1;
                groupCount = groupCount + 1;

            end     %% while groupCount

        end     %% if remainingSpheres
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        cellInfo(p).xTemp = xTemp;
        cellInfo(p).yTemp = yTemp;
        cellInfo(p).zTemp = zTemp;

    else

        cellInfo(p).xTemp = [];
        cellInfo(p).yTemp = [];
        cellInfo(p).zTemp = [];
    end

end     % for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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













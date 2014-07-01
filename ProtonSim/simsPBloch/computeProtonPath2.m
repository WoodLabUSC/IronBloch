function  [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath2(protonPath, sphereInfo, sim_box_side, step, D, hepatocyte_dim, cellBiasFlag, sinusoidBoundaryFlag)
%% 30 October, 2006
%% Nilesh Ghugre, CHLA/USC
%% Restricted proton motion bounded by hepatocyte walls

%% random initialization of rand
rand('twister',sum(100*clock));

a = -sim_box_side/2;
b =  sim_box_side/2;

%tic
insideSphereFlag = 1;
while insideSphereFlag == 1


    % This point should not be inside a sphere.
    init.x = a + (b-a) * rand(1);
    init.y = a + (b-a) * rand(1);
    init.z = a + (b-a) * rand(1);

    [insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(init,sphereInfo);

end

%% if proton is inside sinusoid, place it on the surface
    if ( (cellBiasFlag == 4) & (sinusoidBoundaryFlag == 1) )       %% sinusoids present, prevent protons from entering them.
        [init] = environCrossSinusoid(init,hepatocyte_dim);
    end


%toc
%figure;plot(sphereInfo.radius,surfaceDist);

%% store initial position
protonPath.x(1) = init.x;
protonPath.y(1) = init.y;
protonPath.z(1) = init.z;

%% Find out which hepatocyte the proton is in. Proton motion will be
%% restricted by the walls of this hepatocyte.

%%%%%%% centers of cells define their boundary

numNucleiInX = sim_box_side/hepatocyte_dim; % in single dim

[x,y,z] = meshgrid(-sim_box_side/2:sim_box_side/2);

[centerX,centerY,centerZ] = meshgrid(-(sim_box_side/2)+(sim_box_side/numNucleiInX)/2:sim_box_side/numNucleiInX:sim_box_side/2);

centerX = centerX(:);
centerY = centerY(:);
centerZ = centerZ(:);

%%%%% Find proton distance from all the cell centers. Proton lies in cell
%%%%% whose center is nearest.

cellCenterDist = (sqrt( (protonPath.x(1)-centerX).^2 + (protonPath.y(1)-centerY).^2 + (protonPath.z(1)-centerZ).^2 ));

[minCellDist,minCellDistIndx] = min(cellCenterDist);
% proton lies in the minCellDistIndx'th cell with center
% centerX(minCellDistIndx),centerY(minCellDistIndx),centerZ(minCellDistIndx)

protonCellIndx = minCellDistIndx;   %% renaming

%% cell bounds

aX = (centerX(protonCellIndx)-hepatocyte_dim/2);
bX = (centerX(protonCellIndx)+hepatocyte_dim/2);

aY = (centerY(protonCellIndx)-hepatocyte_dim/2);
bY = (centerY(protonCellIndx)+hepatocyte_dim/2);

aZ = (centerZ(protonCellIndx)-hepatocyte_dim/2);
bZ = (centerZ(protonCellIndx)+hepatocyte_dim/2);

% plotProtonPath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sort spheres based on distance of proton to surface
[currSurfaceDist,currSurfaceDistIndx] = sort(surfaceDist);

%% find which spheres are within the bounded proton path, limited by 
%% hepatocyte dimensions
nearCellSpheresIndx = find(surfaceDist < (2*hepatocyte_dim));

sphereInfoNew.radius = sphereInfo.radius(nearCellSpheresIndx);
sphereInfoNew.x = sphereInfo.x(nearCellSpheresIndx);
sphereInfoNew.y = sphereInfo.y(nearCellSpheresIndx);
sphereInfoNew.z = sphereInfo.z(nearCellSpheresIndx);


%% consider sub-interval for expected proton motion
subTimeInterval = 0.1;     % ms
subTimeIntervalIndx = subTimeInterval/step;
exptProtonMotion = sqrt(6*D*subTimeInterval)*2;     % um
%% consider 2 times the motion for worst case
exptProtonMotion = exptProtonMotion * 2;

[insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(init,sphereInfoNew);
%% find which spheres are within the proton path during the sub-interval
nearSpheresIndx = find(surfaceDist < exptProtonMotion);

% initialize flags
insideSphereFlagTimes = 0;
environCrossFlagTimes = 0;

for tIndx = 2:length(protonPath.x)

    %% calculate next step
    currPos.x = protonPath.x(tIndx-1);
    currPos.y = protonPath.y(tIndx-1);
    currPos.z = protonPath.z(tIndx -1);
    nextPos = randomNextStep(currPos,step,D);

    %% Check if environment boundary has been crossed, if so, place the
    %% proton on the surface of the cell wall
    [nextPos] = environCross2(nextPos,aX,bX,aY,bY,aZ,bZ);
    
    
    if (~mod(tIndx,100))    %% only check for every 100 steps to reduce computation time.
        
        if ( (cellBiasFlag == 4) & (sinusoidBoundaryFlag == 1) )       %% sinusoids present, prevent protons from entering them.
            [nextPos] = environCrossSinusoid(nextPos,hepatocyte_dim);
        end

    end
    
    if (rem(tIndx,subTimeIntervalIndx) == 0)

        %% Calculate distances to all spheres and check if inside any sphere
        [insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(nextPos,sphereInfoNew);
        [currSurfaceDist,currSurfaceDistIndx] = sort(surfaceDist);
        %% find which spheres are within the proton path during the sub-interval
        nearSpheresIndx = find(surfaceDist < exptProtonMotion);
    else

        if (isempty(nearSpheresIndx))

            insideSphereFlag = 0;

        else

            [insideSphereFlag,insideWhichSphere] = calculateNearSphereDist(nextPos,sphereInfoNew,nearSpheresIndx);

        end

    end


    if insideSphereFlag == 1

        insideSphereFlagTimes = insideSphereFlagTimes+1;

        %% placing protons inside-sphere to on-sphere
        [theta_1,phi_1,r_1]=cart2sph(nextPos.x-sphereInfoNew.x(insideWhichSphere),nextPos.y-sphereInfoNew.y(insideWhichSphere),nextPos.z-sphereInfoNew.z(insideWhichSphere));
        r_1(r_1<sphereInfoNew.radius(insideWhichSphere))=sphereInfoNew.radius(insideWhichSphere);
        [nextPos.x,nextPos.y,nextPos.z]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
        nextPos.x = nextPos.x + sphereInfoNew.x(insideWhichSphere);
        nextPos.y = nextPos.y + sphereInfoNew.y(insideWhichSphere);
        nextPos.z = nextPos.z + sphereInfoNew.z(insideWhichSphere);
        clear theta_1 phi_1 r_1;

    end

    protonPath.x(tIndx) = nextPos.x;
    protonPath.y(tIndx) = nextPos.y;
    protonPath.z(tIndx) = nextPos.z;

end


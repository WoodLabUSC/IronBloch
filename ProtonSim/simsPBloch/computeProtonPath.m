function  [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath(protonPath, sphereInfo, sim_box_side, step, D, cellBiasFlag, sinusoidBoundaryFlag, hepatocyte_dim)

%% random initialization of rand
rand('twister',sum(100*clock));

a = -sim_box_side/2;
b =  sim_box_side/2;

%% consider sub-interval for expected proton motion
subTimeInterval = 0.1;     % ms
subTimeIntervalIndx = subTimeInterval/step;
exptProtonMotion = sqrt(6*D*subTimeInterval)*2;     % um
%% consider 2 times the motion for worst case
exptProtonMotion = exptProtonMotion * 2;

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

%% sort spheres based on distance of proton to surface
[currSurfaceDist,currSurfaceDistIndx] = sort(surfaceDist);
%% find which spheres are within the proton path during the sub-interval
nearSpheresIndx = find(surfaceDist < exptProtonMotion);
% figure;plot(currSurfaceDist);

% initialize flags
insideSphereFlagTimes = 0;
environCrossFlagTimes = 0;

for tIndx = 2:length(protonPath.x)

    %% calculate next step
    currPos.x = protonPath.x(tIndx-1);
    currPos.y = protonPath.y(tIndx-1);
    currPos.z = protonPath.z(tIndx -1);
    nextPos = randomNextStep(currPos,step,D);

    %% Check if environment boundary has been crossed, if so wrap, else
    %% leave unchanged
    [nextPos,environCrossFlag] = environCross(nextPos,sim_box_side);
    
    
    if (~mod(tIndx,100))        %% only check for every 100 steps to reduce computation time.
        
        if ( (cellBiasFlag == 4) & (sinusoidBoundaryFlag == 1) )       %% sinusoids present, prevent protons from entering them.
            [nextPos] = environCrossSinusoid(nextPos,hepatocyte_dim);
        end

    end

    if environCrossFlag == 1

        environCrossFlagTimes = environCrossFlagTimes+1;

        %% Calculate distances to all spheres and check if inside any sphere
        [insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(nextPos,sphereInfo);
        [currSurfaceDist,currSurfaceDistIndx] = sort(surfaceDist);
        %% find which spheres are within the proton path during the sub-interval
        nearSpheresIndx = find(surfaceDist < exptProtonMotion);
    end

    if (rem(tIndx,subTimeIntervalIndx) == 0)

        %% Calculate distances to all spheres and check if inside any sphere
        [insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(nextPos,sphereInfo);
        [currSurfaceDist,currSurfaceDistIndx] = sort(surfaceDist);
        %% find which spheres are within the proton path during the sub-interval
        nearSpheresIndx = find(surfaceDist < exptProtonMotion);
    else

        if (isempty(nearSpheresIndx))

            insideSphereFlag = 0;

        else

            [insideSphereFlag,insideWhichSphere] = calculateNearSphereDist(nextPos,sphereInfo,nearSpheresIndx);

        end

    end


    if insideSphereFlag == 1

        insideSphereFlagTimes = insideSphereFlagTimes+1;

        %% placing protons inside-sphere to on-sphere
        [theta_1,phi_1,r_1]=cart2sph(nextPos.x-sphereInfo.x(insideWhichSphere),nextPos.y-sphereInfo.y(insideWhichSphere),nextPos.z-sphereInfo.z(insideWhichSphere));
        r_1(r_1<sphereInfo.radius(insideWhichSphere))=sphereInfo.radius(insideWhichSphere);
        [nextPos.x,nextPos.y,nextPos.z]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
        nextPos.x = nextPos.x + sphereInfo.x(insideWhichSphere);
        nextPos.y = nextPos.y + sphereInfo.y(insideWhichSphere);
        nextPos.z = nextPos.z + sphereInfo.z(insideWhichSphere);
        clear theta_1 phi_1 r_1;

    end

    protonPath.x(tIndx) = nextPos.x;
    protonPath.y(tIndx) = nextPos.y;
    protonPath.z(tIndx) = nextPos.z;

end


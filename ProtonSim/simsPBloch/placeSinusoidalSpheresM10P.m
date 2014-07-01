function [xTemp, yTemp, zTemp] = placeSinusoidalSpheresM10P(sinusoidCount,numSpheresSinusoidEach,xTemp,yTemp,zTemp,R_select,cylinderDiameter,nnCDF,nnDist,NNfactorS)


for sphereCount = 1:numSpheresSinusoidEach(sinusoidCount)

    sphereCount;
    %% placing the first sphere
    if sphereCount == 1


        %% place the sphere anywhere in the 10um diameter cylindrical
        %% region with long axis passing thru origin and parallel to z
        %% axis.
        a = zTemp(sinusoidCount).limits(1) + R_select(sinusoidCount).r(sphereCount);
        b = zTemp(sinusoidCount).limits(2) - R_select(sinusoidCount).r(sphereCount);
        zTemp(sinusoidCount).z(sphereCount) = a + (b-a) * rand(1);

        axialR = 0 + (cylinderDiameter/2-0)*rand(1);
        axialTheta = (pi/180)*(0+(360-0)*rand(1)); % random angle in degrees, convert to radians
        xTemp(sinusoidCount).x(sphereCount) = axialR*cos(axialTheta) + xTemp(sinusoidCount).Offset;
        yTemp(sinusoidCount).y(sphereCount) = axialR*sin(axialTheta) + yTemp(sinusoidCount).Offset;


    else

        %% choose previous sphere placed
        selectedSphIndx = sphereCount - 1;

        %% if selected sphere is zero coordinates, it has been marked to be placed in the hepatocyte, hence choose
        %% another random one.
        % selectedSphIndx = ceil(1 + ((sphereCount - 1)-1) * rand(1));

        while (xTemp(sinusoidCount).x(selectedSphIndx) == 0)
            %% choose any sphere already placed.
            selectedSphIndx = round(1 + ((sphereCount - 1)-1) * rand(1));
        end

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

        if currDist < (R_select(sinusoidCount).r(selectedSphIndx)+R_select(sinusoidCount).r(sphereCount))
            currDist = (R_select(sinusoidCount).r(selectedSphIndx)+R_select(sinusoidCount).r(sphereCount))*(1+0.05);   %% 5% more than sticking distance
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

        cannotAccommodateCount = 1;

        while (overlapFlag == 1 || environCrossFlag == 1) && (cannotAccommodateCount <= 20)

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

                if currDist < (R_select(sinusoidCount).r(selectedSphIndx)+R_select(sinusoidCount).r(sphereCount))
                    currDist = (R_select(sinusoidCount).r(selectedSphIndx)+R_select(sinusoidCount).r(sphereCount))*(1+0.05);   %% 5% more than sticking distance
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
            xTemp(sinusoidCount).x(sphereCount) = xTemp(sinusoidCount).x(selectedSphIndx) + displacementX;
            yTemp(sinusoidCount).y(sphereCount) = yTemp(sinusoidCount).y(selectedSphIndx) + displacementY;
            zTemp(sinusoidCount).z(sphereCount) = zTemp(sinusoidCount).z(selectedSphIndx) + displacementZ;

            %% position has to be checked if sphere is going outside
            %% cylindrical sinusoidal boundary




            %% offset needs to be subtracted to get back to center of
            %% volume.

            if sum(sinusoidCount == [1:5 11:14])      %% positive z

                if ((sqrt( (xTemp(sinusoidCount).x(sphereCount) - xTemp(sinusoidCount).Offset)^2 + (yTemp(sinusoidCount).y(sphereCount) - yTemp(sinusoidCount).Offset)^2 )) > (cylinderDiameter/2)) || (zTemp(sinusoidCount).z(sphereCount)<(zTemp(sinusoidCount).limits(1))) || (zTemp(sinusoidCount).z(sphereCount)>(zTemp(sinusoidCount).limits(2)))
                    environCrossFlag = 1;
                else
                    environCrossFlag = 0;
                end
            else                                    %% nogative z
                if ((sqrt( (xTemp(sinusoidCount).x(sphereCount) - xTemp(sinusoidCount).Offset)^2 + (yTemp(sinusoidCount).y(sphereCount) - yTemp(sinusoidCount).Offset)^2 )) > (cylinderDiameter/2)) || (zTemp(sinusoidCount).z(sphereCount)>(zTemp(sinusoidCount).limits(1))) || (zTemp(sinusoidCount).z(sphereCount)<(zTemp(sinusoidCount).limits(2)))
                    environCrossFlag = 1;
                else
                    environCrossFlag = 0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%

            if environCrossFlag == 0

                numPlaceOnSurface = 0;
                while numPlaceOnSurface < 10

                    %% Check for overlap with previously placed spheres
                    for j = 1:sphereCount-1

                        centerDist = sqrt( (xTemp(sinusoidCount).x(sphereCount)-xTemp(sinusoidCount).x(j)).^2 + (yTemp(sinusoidCount).y(sphereCount)-yTemp(sinusoidCount).y(j)).^2 + (zTemp(sinusoidCount).z(sphereCount)-zTemp(sinusoidCount).z(j)).^2 );

                        if centerDist > (R_select(sinusoidCount).r(sphereCount)+R_select(sinusoidCount).r(j))
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
                        [theta_1,phi_1,r_1]=cart2sph(xTemp(sinusoidCount).x(sphereCount)-xTemp(sinusoidCount).x(insideWhichSphere),yTemp(sinusoidCount).y(sphereCount)-yTemp(sinusoidCount).y(insideWhichSphere),zTemp(sinusoidCount).z(sphereCount)-zTemp(sinusoidCount).z(insideWhichSphere));

                        %% place current sphere on the surface of collided sphere
                        r_1 = R_select(sinusoidCount).r(sphereCount)+R_select(sinusoidCount).r(insideWhichSphere);

                        %% reconvert
                        [xTemp(sinusoidCount).x(sphereCount),yTemp(sinusoidCount).y(sphereCount),zTemp(sinusoidCount).z(sphereCount)]=sph2cart(theta_1,phi_1,r_1);  %% reconverting
                        %% getting back into original coordinate system
                        xTemp(sinusoidCount).x(sphereCount) = xTemp(sinusoidCount).x(sphereCount)+xTemp(sinusoidCount).x(insideWhichSphere);
                        yTemp(sinusoidCount).y(sphereCount) = yTemp(sinusoidCount).y(sphereCount)+yTemp(sinusoidCount).y(insideWhichSphere);
                        zTemp(sinusoidCount).z(sphereCount) = zTemp(sinusoidCount).z(sphereCount)+zTemp(sinusoidCount).z(insideWhichSphere);

                        numPlaceOnSurface = numPlaceOnSurface + 1;

                    else
                        break;
                    end


                end     %% while numPlaceOnSurface


            end     %% if environCrossFlag

            repeatDistanceCount = repeatDistanceCount + 1;
            cannotAccommodateCount = cannotAccommodateCount + 1;

        end     %% while overlapFlag

        if cannotAccommodateCount > 20

            xTemp(sinusoidCount).x(sphereCount) = 0;
            yTemp(sinusoidCount).y(sphereCount) = 0;
            zTemp(sinusoidCount).z(sphereCount) = 0;
        end

        % sphereCount = sphereCount + 1;


    end     % if sphereCount


    %             if ~mod(sphereCount,100)
    %                 sphereCount
    %             end


end     % for sphereCount

function [xTemp,yTemp,zTemp] = placeSpheresM2(xTemp,yTemp,zTemp,sim_box_side,R_select);
%% Spheres are placed in a box in a non-overlapping manner.
%% If overlap occurs, the sphere is placed at the surface of colliding
%% sphere and then checked again for overlap. The procedure is continued
%% for successive overlaps and if there is a deadlock, a new position is
%% generated and the procedure is repeated until no overlap.
%% 17 Oct. 2006
%% Nilesh Ghugre, CHLA/USC

%% random initialization of rand
rand('twister',sum(100*clock));

numSpheres = length(xTemp);

%% place the first sphere randomly
p = 1;
a = -sim_box_side/2+R_select(p);
b =  sim_box_side/2-R_select(p);

xTemp(p) = a + (b-a) * rand(1);
yTemp(p) = a + (b-a) * rand(1);
zTemp(p) = a + (b-a) * rand(1);

for p=2:numSpheres

    overlapFlag = 1;

    while overlapFlag == 1


        %% new box coordinates such that sphere does not lie outside box dims
        a = -sim_box_side/2+R_select(p);
        b =  sim_box_side/2-R_select(p);

        xTemp(p) = a + (b-a) * rand(1);
        yTemp(p) = a + (b-a) * rand(1);
        zTemp(p) = a + (b-a) * rand(1);


        numPlaceOnSurface = 0;

        while numPlaceOnSurface < 10

            %% Check for overlap with previously placed spheres
            for j = (p-1):-1:1

                centerDist = sqrt( (xTemp(p)-xTemp(j)).^2 + (yTemp(p)-yTemp(j)).^2 + (zTemp(p)-zTemp(j)).^2 );

                if centerDist > (R_select(p)+R_select(j))
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

    end     % while overlap
end     % for p



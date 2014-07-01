function [sphereInfo, sim_box_side, sphereVolFrac] = GetSpheres1(patientInfo, patientIndx, numSpheres)
%% 24 Oct 2006
%% Nilesh Ghugre, CHLA/USC
%% For given volume fraction and number of spheres, random positions for
%% spheres are generated within the simulation volume.
%% Specified number of spheres are generated, given the volume fraction.
%% Simulation box size depends on number of spheres.
%% Sphere volume fraction is the same as specified in patientInfo

R_highres = linspace(min(patientInfo(patientIndx).r),max(patientInfo(patientIndx).r),100);

% tic
for k = 1      %% how many times to generate

    rand_num = rand(1,numSpheres);

    % R_select contains the radii of chosen spheres
    for indx = 1:length(rand_num)
        [Y,I] = min(abs(patientInfo(patientIndx).dist_est_cumsum-rand_num(indx)));
        R_select(indx) = R_highres(I);
    end

    % figure; hist(R_select,100);
    min_R_select(k) = min(R_select);
    max_R_select(k) = max(R_select);
    total_vol_es(k) = sum((4*pi/3) * (R_select.^3));
    totalSphereVol = total_vol_es;
end
% toc

% figure; hist(R_select,100);

%% assuming X% vol frac obtained from lysosome analysis,
%% simlulation volume in um^3 and box side in um
sphereVolFrac = (patientInfo(patientIndx).volFrac)/100;
sim_volume = mean(total_vol_es) / (sphereVolFrac);
sim_box_side =  (sim_volume)^(1/3);

%% Sorting sphere radii so that largest is the first.
R_select = sort(R_select,'descend');

% choose a random poisition for each sphere inside sim volume 
% taking care that they do not protrude outside the sim volume.

%% memory allocations for sphere centers
xTemp = zeros(1,numSpheres);
yTemp = zeros(1,numSpheres);
zTemp = zeros(1,numSpheres);

%tic
[xTemp,yTemp,zTemp] = placeSpheresM2(xTemp,yTemp,zTemp,sim_box_side,R_select);
%toc

%% 3D visualization of sphere positions
% SphereVisualize(R_select,xTemp,yTemp,zTemp,sim_box_side);

%% Store in a structure, radius and center coordinate

sphereInfo.radius = R_select;
sphereInfo.x = xTemp;
sphereInfo.y = yTemp;
sphereInfo.z = zTemp;



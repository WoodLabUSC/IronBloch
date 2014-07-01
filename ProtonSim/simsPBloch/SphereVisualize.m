function SphereVisualize(R_select,xTemp,yTemp,zTemp,sim_box_side)
%% Visualize spheres in a box
%% based as in
%% Daneel D:\Nilesh-Backup\Daneel-Backup\MonteCarlo-2\MonteCarlo_April2006\tissue
%% model\test3.m
%% 17 Oct. 2006
%% Nilesh Ghugre


fieldGridStep = 0.02*sim_box_side;    % um for every grid point, for single use 0.03
spill = 0;      % um, to avoid NAN's during interpolation of boundary points, for single use 0.5
[X,Y,Z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);


bw = zeros(size(X));

for i=1:length(R_select)
    bw1 = zeros(size(X));
    bw1 = sqrt((X-xTemp(i)).^2 + (Y-yTemp(i)).^2 + (Z-zTemp(i)).^2) <= R_select(i);
    bw = bw | bw1;
    clear bw1;
end


lims = [-sim_box_side/2 sim_box_side/2];

figure, isosurface(X,Y,Z,bw,0), axis equal
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud
grid on;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_box_side = 100;
% radius_nuclei = 2;    % in um
% hepatocyte_dim = 20;    % in um
% numNucleiInX = sim_box_side/hepatocyte_dim;
% 
% [x,y,z] = meshgrid(-sim_box_side/2:sim_box_side/2);
% 
% [centerX,centerY,centerZ] = meshgrid(-(sim_box_side/2)+(sim_box_side/numNucleiInX)/2:sim_box_side/numNucleiInX:sim_box_side/2);
% 
% centerX = centerX(:);
% centerY = centerY(:);
% centerZ = centerZ(:);
% 
% fieldGridStep = 0.01*sim_box_side;    % um for every grid point, for single use 0.03
% spill = 0;      % um, to avoid NAN's during interpolation of boundary points, for single use 0.5
% [x,y,z] = meshgrid(-sim_box_side/2-spill:fieldGridStep:sim_box_side/2+spill);
% 
% bw = zeros(size(x));
% 
% 
% for i=1:length(centerX)
%     bw1 = zeros(sim_box_side,sim_box_side,sim_box_side);
%     bw1 = sqrt((x-centerX(i)).^2 + (y-centerY(i)).^2 + (z-centerZ(i)).^2) <= radius_nuclei;
%     bw = bw | bw1;
%     clear bw1;
% end
% 
% 
% lims = [-sim_box_side/2 sim_box_side/2];
% 
% figure, isosurface(x,y,z,bw,0.5), axis equal
% xlabel x, ylabel y, zlabel z
% xlim(lims), ylim(lims), zlim(lims)
% view(3), camlight, lighting gouraud

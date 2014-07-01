
sim_box_side = 80;
% Visualize sphere positions
lims = [-sim_box_side/2 sim_box_side/2];

figure; hold on;
grid; 
xlim(lims), ylim(lims), zlim(lims)
plot3(sphereInfo.x,sphereInfo.y,sphereInfo.z,'.','markersize',6);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', []);
box on;


% for num = 1:20
% 
%     rand('twister',sum(100*clock));
% 
%     protonPath.x = zeros(1,length(t));
%     protonPath.y = zeros(1,length(t));
%     protonPath.z = zeros(1,length(t));
% 
%     [protonPath, environCrossFlagTimes, insideSphereFlagTimes] = computeProtonPath2(protonPath, sphereInfo, sim_box_side, step, D, hepatocyte_dim, cellBiasFlag, sinusoidBoundaryFlag);
% 
%     plot3(protonPath.x,protonPath.y,protonPath.z,'.g','markersize',2);
% 
% end
% 
% 

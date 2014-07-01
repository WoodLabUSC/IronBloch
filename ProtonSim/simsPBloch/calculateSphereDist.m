function [insideSphereFlag,insideWhichSphere,surfaceDist] = calculateSphereDist(nextPos,sphereInfo)

xInit = nextPos.x;
yInit = nextPos.y;
zInit = nextPos.z;

insideWhichSphere = 0;
insideSphereFlag = 0;

% Compute distance of proton from surface of each sphere

surfaceDist = zeros(1,length(sphereInfo.radius));
surfaceDist = (sqrt( (xInit-sphereInfo.x).^2 + (yInit-sphereInfo.y).^2 + (zInit-sphereInfo.z).^2 )) - sphereInfo.radius;
surfaceDistInd = (surfaceDist<0);
surfaceDistInd = find(surfaceDistInd);
if min(size(surfaceDistInd)) > 0
    insideWhichSphere = surfaceDistInd(1);
    insideSphereFlag = 1;
else
    insideSphereFlag = 0;
end
%for i = 1:length(sphereInfo.radius)
    %surfaceDist(i) = (sqrt( (xInit-sphereInfo.x(i)).^2 + (yInit-sphereInfo.y(i)).^2 + (zInit-sphereInfo.z(i)).^2 )) - sphereInfo.radius(i);
%    if surfaceDist(i) < 0
%        % point is inside a sphere, break for loop
%        insideWhichSphere = i;
%        insideSphereFlag = 1;
%        break;
%    else
%        % point is not inside a sphere, continue comparison
%        insideSphereFlag = 0;
%    end
%end



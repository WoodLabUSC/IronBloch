function [insideSphereFlag,insideWhichSphere] = calculateNearSphereDist(nextPos,sphereInfo,nearSpheresIndx)

xInit = nextPos.x;
yInit = nextPos.y;
zInit = nextPos.z;

insideWhichSphere = 0;
insideSphereFlag = 0;

% Compute distance of proton from surface of nearest spheres

surfaceDist = zeros(1,length(nearSpheresIndx));
surfaceDist = (sqrt( (xInit-sphereInfo.x(nearSpheresIndx)).^2 + (yInit-sphereInfo.y(nearSpheresIndx)).^2 + (zInit-sphereInfo.z(nearSpheresIndx)).^2 )) - sphereInfo.radius(nearSpheresIndx);
surfaceDistInd = (surfaceDist<0);
nearSpheresIndx = nearSpheresIndx(surfaceDistInd);
if min(size(nearSpheresIndx)) > 0
    insideWhichSphere = nearSpheresIndx(1);
    insideSphereFlag = 1;
end
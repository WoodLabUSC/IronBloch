function [nextPos,environCrossFlag] = environCross3(nextPos,aX,bX,aY,bY,aZ,bZ)
%% Checking if a sphere placed is within a given cell.

environCrossFlag = 0;       %% no crossing

if (nextPos.x < aX)
    nextPos.x = aX;
    environCrossFlag = 1;
elseif (nextPos.x > bX)
    nextPos.x = bX;
    environCrossFlag = 1;
end

if (nextPos.y < aY)
    nextPos.y = aY;
    environCrossFlag = 1;
elseif (nextPos.y > bY)
    nextPos.y = bY;
    environCrossFlag = 1;
end

if (nextPos.z < aZ)
    nextPos.z = aZ;
    environCrossFlag = 1;
elseif (nextPos.z > bZ)
    nextPos.z = bZ;
    environCrossFlag = 1;
end



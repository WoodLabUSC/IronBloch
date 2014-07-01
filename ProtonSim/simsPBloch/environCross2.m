function [nextPos] = environCross2(nextPos,aX,bX,aY,bY,aZ,bZ)

environCrossFlag = 0;       %% no crossing

if (nextPos.x < aX)
    nextPos.x = aX;
elseif (nextPos.x > bX)
    nextPos.x = bX;
end

if (nextPos.y < aY)
    nextPos.y = aY;
elseif (nextPos.y > bY)
    nextPos.y = bY;
end

if (nextPos.z < aZ)
    nextPos.z = aZ;
elseif (nextPos.z > bZ)
    nextPos.z = bZ;
end



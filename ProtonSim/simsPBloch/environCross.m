function [nextPos,environCrossFlag] = environCross(nextPos,sim_box_side)

environCrossFlag = 0;       %% no crossing

if (nextPos.x < 0) & (nextPos.x < (-sim_box_side/2))
    % wrap
    nextPos.x = nextPos.x + sim_box_side;
    environCrossFlag = 1;
end
if (nextPos.x > 0) & (nextPos.x > (sim_box_side/2))
    % wrap
    nextPos.x = nextPos.x - sim_box_side;
    environCrossFlag = 1;
end

if (nextPos.y < 0) & (nextPos.y < (-sim_box_side/2))
    % wrap
    nextPos.y = nextPos.y + sim_box_side;
    environCrossFlag = 1;
end
if (nextPos.y > 0) & (nextPos.y > (sim_box_side/2))
    % wrap
    nextPos.y = nextPos.y - sim_box_side;
    environCrossFlag = 1;
end

if (nextPos.z < 0) & (nextPos.z < (-sim_box_side/2))
    % wrap
    nextPos.z = nextPos.z + sim_box_side;
    environCrossFlag = 1;
end
if (nextPos.z > 0) & (nextPos.z > (sim_box_side/2))
    % wrap
    nextPos.z = nextPos.z - sim_box_side;
    environCrossFlag = 1;
end


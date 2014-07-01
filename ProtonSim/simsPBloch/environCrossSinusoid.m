function [nextPos] = environCrossSinusoid(nextPos,hepatocyte_dim)
%% check if proton is inside a sinusoidal region, if yes, then place the
%% proton on it. the sinusoidal area bounds are approximated by rectangular
%% regions around the cylindrical sinusoids.

%% offset coordinates for the sinusoids from the center of the simulation
%% volume.
%% These are defined in placeSpheresM9P2.m

%% level 1

xTemp(1).Offset = 0;
yTemp(1).Offset = 0;

xTemp(2).Offset = -hepatocyte_dim;
yTemp(2).Offset = hepatocyte_dim;

xTemp(3).Offset = hepatocyte_dim;
yTemp(3).Offset = hepatocyte_dim;

xTemp(4).Offset = hepatocyte_dim;
yTemp(4).Offset = -hepatocyte_dim;

xTemp(5).Offset = -hepatocyte_dim;
yTemp(5).Offset = -hepatocyte_dim;

%% z limits are same for all in the same level
for k = 1:5
    zTemp(k).limits = [hepatocyte_dim 2*hepatocyte_dim];
end

%% level 3, same as level 1

xTemp(6).Offset = 0;
yTemp(6).Offset = 0;

xTemp(7).Offset = -hepatocyte_dim;
yTemp(7).Offset = hepatocyte_dim;

xTemp(8).Offset = hepatocyte_dim;
yTemp(8).Offset = hepatocyte_dim;

xTemp(9).Offset = hepatocyte_dim;
yTemp(9).Offset = -hepatocyte_dim;

xTemp(10).Offset = -hepatocyte_dim;
yTemp(10).Offset = -hepatocyte_dim;

%% z limits are same for all in the same level
for k = 6:10
    zTemp(k).limits = [0 -hepatocyte_dim];
end

%% level 2

xTemp(11).Offset = 0;
yTemp(11).Offset = -hepatocyte_dim;

xTemp(12).Offset = 0;
yTemp(12).Offset = hepatocyte_dim;

xTemp(13).Offset = hepatocyte_dim;
yTemp(13).Offset = 0;

xTemp(14).Offset = -hepatocyte_dim;
yTemp(14).Offset = 0;

%% z limits are same for all in the same level
for k = 11:14
    zTemp(k).limits = [0 hepatocyte_dim];
end

%% level 4, same as level 2

xTemp(15).Offset = 0;
yTemp(15).Offset = -hepatocyte_dim;

xTemp(16).Offset = 0;
yTemp(16).Offset = hepatocyte_dim;

xTemp(17).Offset = hepatocyte_dim;
yTemp(17).Offset = 0;

xTemp(18).Offset = -hepatocyte_dim;
yTemp(18).Offset = 0;

%% z limits are same for all in the same level
for k = 15:18
    zTemp(k).limits = [-hepatocyte_dim -hepatocyte_dim*2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numberOfSinusoids = 18;
cylinderDiameter = 10;          %% um, book reference from Ignacio

%% find out which sinusoid proton is inside.
for sinusoidCount = 1:numberOfSinusoids
   
    x_limit_low = xTemp(sinusoidCount).Offset - cylinderDiameter/2;
    x_limit_high = xTemp(sinusoidCount).Offset + cylinderDiameter/2;
    
    y_limit_low = yTemp(sinusoidCount).Offset - cylinderDiameter/2;
    y_limit_high = yTemp(sinusoidCount).Offset + cylinderDiameter/2;    
    
    z_limit_low = zTemp(sinusoidCount).limits(1);
    z_limit_high = zTemp(sinusoidCount).limits(2);
    
    if ( (nextPos.x > x_limit_low) & (nextPos.x < x_limit_high) & (nextPos.y > y_limit_low) & (nextPos.y < y_limit_high) & (nextPos.z > z_limit_low) & (nextPos.z < z_limit_high))
       
        %% Place proton on the surface of sinusoid at the nearest entry
        %% point.
        
        %% x
        if ( ( abs(nextPos.x - x_limit_low) ) < ( abs(nextPos.x - x_limit_high) ) )
            nextPos.x = x_limit_low;
        else
            nextPos.x = x_limit_high;
        end
        
        %% y
        if ( ( abs(nextPos.y - y_limit_low) ) < ( abs(nextPos.y - y_limit_high) ) )
            nextPos.y = y_limit_low;
        else
            nextPos.y = y_limit_high;
        end
        
        %% z
        if ( ( abs(nextPos.z - z_limit_low) ) < ( abs(nextPos.z - z_limit_high) ) )
            nextPos.z = z_limit_low;
        else
            nextPos.z = z_limit_high;
        end        
        
        break;      %% no need to iterate through other sinusoids since proton can be present in only one of them.
    end
        
end





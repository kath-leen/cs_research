% This function gets the coordinates (xa, ya) of the point "A" based on the
% coordinates (xb, yb) and (xc, yc) of the points "B" anc "C" and the "AB" 
% and "AC" distances.
% "ABC" is the triangle where the distances "AB and "AC" and the 
% coordinates of "B" anc "C" are known. All points are located on the 
% two-dimentional plane and therefore have two coordinates.
% Argument "point_a_desired" is optional and refers to the supposed point
% "A" coordinates.

function [ x, y, ok ] = get_coordinates_from_distances( dist_ab, dist_ac, point_b, point_c, point_a_desired )

[x, y, ok] = get_coordinates_from_distances_3d([point_b 0], [point_c 0], dist_ab, dist_ac, 0);

if length(x) == 1
    return;
elseif length(x) == 2
    
    if nargin < 5
        point_a_desired = [0, 0];
    end
    
    err1 = get_distance(point_a_desired, [x(1) y(1)]);
    err2 = get_distance(point_a_desired, [x(2) y(2)]);
    
    if (err1 < err2)
        x = x(1);
        y = y(1);
    else
        x = x(2);
        y = y(2);
    end
else
    return;
end

end


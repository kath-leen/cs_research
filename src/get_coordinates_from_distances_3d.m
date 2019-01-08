% This function gets the coordinates (xa, ya) of the point "A" based on the
% coordinates (xb, yb, zb) and (xc, yc, zc) of the points "B" anc "C" and 
% the "AB" and "AC" distances.
% "ABC" is the triangle where the distances "AB and "AC" and the 
% coordinates of "B" anc "C" are known. The third coordinate (za) of the 
% point "A" is also known. All points are located on the three-dimentional 
% plane and therefore have three coordinates.

function [xa, ya, ok] = get_coordinates_from_distances_3d(point_b, point_c, dist_ab, dist_ac, za)

if point_b == point_c
    xa = [];
    ya = [];
    ok = 0;
    return;
end

xb = point_b(1);
yb = point_b(2);
zb = point_b(3);
xc = point_c(1);
yc = point_c(2);
zc = point_c(3);

% This system of equations will be solved:
% (xb - xa)^2 + (yb - ya)^2 + (zb - za)^2 = dist_ab^2
% (xc - xa)^2 + (yc - 2)^2 + (zb - za)^2 = dist_ac^2

% If yb ~= yc, the system will be solved for x
if (yb ~= yc)
    A = (dist_ab^2 - dist_ac^2 - xb^2 + xc^2 - yb^2 + yc^2 - zb^2 + zc^2 - 2*za*(zb + zc)) / (-2*(yb - yc));
    B = (xb - xc) / (yb - yc);

    a = 1 + B^2;
    b = -2*xc + 2*B*yc - 2*A*B;
    c = xc^2 + yc - 2*yc*A + A^2 + za^2 - 2*za*zc + zc^2 - dist_ac^2;

    D = b^2 - 4*a*c;
    if (D < 0)
        xa = [];
        ya = [];
        ok = 0;
        return;
    end
    xa = [(-b + sqrt(D))/(2*a) (-b - sqrt(D))/(2*a)];
    ya = A - B*xa;
    if ~D
        xa = xa(1);
        ya = ya(1);
        ok = 1;
    else
        ok = 2;
    end
else
    % otherwise, it will be solved for y
    A = (dist_ab^2 - dist_ac^2 - xb^2 + xc^2 - yb^2 + yc^2 - zb^2 + zc^2 + 2*za*(zb + zc)) / (-2*(xb - xc));
    B = (yb - yc) / (xb - xc);

    a = 1 + B^2;
    b = -2*yb + 2*B*xb - 2*A*B;
    c = xb^2 + yb - 2*xb*A + A^2 + za^2 - 2*za*zc - dist_ab^2;

    D = b^2 - 4*a*c;
    if (D < 0) || imag(D)
        xa = [];
        ya = [];
        ok = 0;
        return;
    end
    ya = [(-b + sqrt(D))/(2*a) (-b - sqrt(D))/(2*a)];
    xa = A - B*ya;
    if ~D
        xa = xa(1);
        ya = ya(1);
        ok = 1;
    else
        ok = 2;
    end
end
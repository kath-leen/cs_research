function [x, y, ok] = get_coord_from_distances(coord_point1, coord_point2, R1, R2, h)

if coord_point1 == coord_point2
    x = [];
    y = [];
    ok = 0;
    return;
end

x1 = coord_point1(1);
y1 = coord_point1(2);
z1 = coord_point1(3);
x2 = coord_point2(1);
y2 = coord_point2(2);
z2 = coord_point2(3);
z = h;

% This system of equations will be solved:
% (x1 - x)^2 + (y1 - y)^2 + (z1 - z)^2 = R1^2
% (x2 - x)^2 + (y2 - 2)^2 + (z1 - z)^2 = R2^2

% If y1 ~= y2, the system will be solved for x
if (y1 ~= y2)
    A = (R1^2 - R2^2 - x1^2 + x2^2 - y1^2 + y2^2 - z1^2 + z2^2 - 2*z*(z1 + z2)) / (-2*(y1 - y2));
    B = (x1 - x2) / (y1 - y2);

    a = 1 + B^2;
    b = -2*x2 + 2*B*y2 - 2*A*B;
    c = x2^2 + y2 - 2*y2*A + A^2 + z^2 - 2*z*z2 + z2^2 - R2^2;

    D = b^2 - 4*a*c;
    if (D < 0)
        x = [];
        y = [];
        ok = 0;
        return;
    end
    x = [(-b + sqrt(D))/(2*a) (-b - sqrt(D))/(2*a)];
    y = A - B*x;
    if ~D
        x = x(1);
        y = y(1);
        ok = 1;
    else
        ok = 2;
    end
else
    % otherwise, it will be solved for y
    A = (R1^2 - R2^2 - x1^2 + x2^2 - y1^2 + y2^2 - z1^2 + z2^2 + 2*z*(z1 + z2)) / (-2*(x1 - x2));
    B = (y1 - y2) / (x1 - x2);

    a = 1 + B^2;
    b = -2*y1 + 2*B*x1 - 2*A*B;
    c = x1^2 + y1 - 2*x1*A + A^2 + z^2 - 2*z*z2 - R1^2;

    D = b^2 - 4*a*c;
    if (D < 0) || imag(D)
        x = [];
        y = [];
        ok = 0;
        return;
    end
    y = [(-b + sqrt(D))/(2*a) (-b - sqrt(D))/(2*a)];
    x = A - B*y;
    if ~D
        x = x(1);
        y = y(1);
        ok = 1;
    else
        ok = 2;
    end
end
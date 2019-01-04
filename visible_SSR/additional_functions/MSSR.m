function [ok, x, y] = MSSR(L, phi, h, rlz)
L = L/2;
x = 0;
y = 0;
if (L^2 < h^2)
    ok = 0;
    return;
end
ok = 1;
l = sqrt(L^2 - h^2);
if (phi < pi/2)
    x = -l*cos(phi);
    y = l*sin(phi);
elseif (phi < pi)
    x = l*cos(pi-phi);
    y = l*sin(pi-phi);
elseif (phi < 3*pi/2)
    x = l*cos(phi-pi);
    y = -l*sin(phi-pi);
else
    x = -l*cos(-phi);
    y = -l*sin(-phi);
end
x = x + rlz(1);
y = y + rlz(2);
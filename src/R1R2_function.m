% This function calculates the distance between the aircraft and the SSR
% (R1) and the distance between the aircraft and the receiver (R2) based on
% the path difference (L), the receiver-SSR-arcraft angle (phi) and the
% aircraft altitude (h)

function [ R1, R2 ] = R1R2_function( L, b, phi, h )

R1 = [];
R2 = [];
h = floor(h/375) * 375; % according to the ICAO standard
%h = h_in + rand()*375 - 375/2;
%h = h_in;

aa = 4*b^2*(cos(phi))^2 - 4*(L+b)^2;
bb = 4*(L + b)*(L^2 + 2*L*b);
cc = -4*b^2*(cos(phi))^2*h^2 - (L^2+2*b*L)^2;

D = bb^2 - 4*aa*cc;
if (D < 0)
    return;
end
R1_candidates = [(-bb + sqrt(D))/(2*aa) (-bb - sqrt(D))/(2*aa)];

R1_estimation = (L^2 + 2*L*b)/(2*L + 2*b - 2*b*cos(phi));
error = abs(R1_candidates - R1_estimation);
[~, ind] = min(error);
R1 = R1_candidates(ind);
R2 = L - R1 + b;

end
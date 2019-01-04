function L = get_distance( c_1, c_2)
%Calculate distance between (x1,y1,z1) and (x2,y2,z2)
L = sqrt((c_1(1) - c_2(1)).^2 + (c_1(2) - c_2(2)).^2 + (c_1(3) - c_2(3)).^2);
end


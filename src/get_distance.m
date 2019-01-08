% This function calculates distance between two points. Each point must
% have 3 coordinates (x, y, z)

function L = get_distance(point_1, point_2)

L = sqrt((point_1(1) - point_2(1)).^2 + (point_1(2) - point_2(2)).^2 + (point_1(3) - point_2(3)).^2);

end


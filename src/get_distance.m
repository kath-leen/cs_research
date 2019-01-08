% This function calculates distance between two points.

function L = get_distance(point_1, point_2)

if length(point_1) > length(point_2)
    point_2 = [point_2 zeros(1, length(point_1) - length(point_2))];
elseif length(point_1) < length(point_2)
    point_1 = [point_1 zeros(1, length(point_2) - length(point_1))];
end

L = sqrt(sum((point_1 - point_2).^2));

end

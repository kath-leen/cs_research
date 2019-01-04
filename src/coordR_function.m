function [ x, y, ok ] = coordR_function( R1_pr, R2_pr, pp, rlz, p_ideal )

[x, y, ok] = get_coord_from_distances([pp(1:2) 0], [rlz(1:2) 0], R1_pr, R2_pr, 0);

p_ideal_pr = [p_ideal(1:2) 0];
if length(x) == 1
    return;
elseif length(x) == 2
    r_p_p_1 = get_distance(p_ideal_pr, [x(1) y(1) 0]);
    r_p_p_2 = get_distance(p_ideal_pr, [x(2) y(2) 0]);
    if (r_p_p_1 < r_p_p_2)
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


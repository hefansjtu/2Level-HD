function [point, step] = bisection(x_init, dir, k ,w)
ubound = 1;
lbound = 0;
step = ubound;
point = x_init + ubound*dir;
while abs(fh_value(x_init, k, w) - fh_value(x_init + ubound*dir, k, w)) < 1e-5
    lbound = ubound;
    ubound = ubound*2;
    if ubound > 9999
        return;
    end
end

while abs(ubound - lbound) > 1e-5
    mid = (lbound + ubound)/2;
    if abs(fh_value(x_init, k, w) - fh_value(x_init + mid*dir, k, w)) < 1e-5
        lbound = mid;
    else
        ubound = mid;
    end
end
step = ubound;
point = x_init + ubound*dir;

end

function [val] = fh_value(x, k, w)
val = norm(x,1) - get_value(x, k, w);
end
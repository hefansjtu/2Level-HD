function [newpoints, steps] = gamma_extension(x_init, fh, k, w, directions)
N = length(x_init);
newpoints = x_init*ones(1,size(directions,2));
for iter = 1: size(directions,2)
    [point, step] = find_newpoint(x_init, fh, k, w, directions(:,iter));
%     [point, step] = find_newpoint_g(x_init, directions(:,iter));
    if norm(step*directions(:,iter),inf) > 1e-3*norm(x_init,inf)
        newpoints(:,iter) = point;
        steps(iter) = step;
    else
        steps(iter) = Inf;
    end
end
end

function [point, step] = find_newpoint(x_init, fh, k, w, dir)
ubound = 100;
lbound = 0;

while cal_fh(x_init + ubound*dir, k, w) < fh
    lbound = ubound;
    ubound = ubound*2;
end

if norm(ubound*dir,Inf) > 1000*norm(x_init, Inf)
    point = [];
    step = 0;
    return;
end

while abs(ubound - lbound) > 1e-10
    mid = (ubound + lbound)/2;
    point = x_init + mid*dir;
    if cal_fh(point, k, w) < fh
        lbound = mid;
    else
        ubound = mid;
    end
end
step = mid;
end

function [point, step] = find_newpoint_g(x_init,  dir)
ubound = 100;
lbound = 0;
fg = norm(x_init,1);
while norm(x_init + ubound*dir, 1) < fg
    lbound = ubound;
    ubound = ubound*2;
end

while abs(ubound - lbound) > 1e-10
    mid = (ubound + lbound)/2;
    point = x_init + mid*dir;
    if norm(point, 1) < fg
        lbound = mid;
    else
        ubound = mid;
    end
end
step = mid;
end

function [val] = cal_fh(x, k, w)
val = norm(x,1) - get_value(x, k, w);
end
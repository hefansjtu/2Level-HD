function [val] = norm0(x)
val = length(x(abs(x)>1e-3*norm(x,inf)));
end
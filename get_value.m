function [val] = get_value(x, k, w, weight)
if nargin == 3
[weight] = get_weight(x,k,w);
end
val = weight'*abs(x);
end
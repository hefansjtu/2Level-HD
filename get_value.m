function [val] = get_value(x, k, w, weight)
if nargin == 3
[weight] = get_weight(x,k,w);
end
val = weight'*abs(x);
end

function [weight] = get_weight(x,k,w)
[~, index] = sort(abs(x),'descend');
weight = ones(length(x),1);
weight(index(1:k)) = w*ones(k,1);
end
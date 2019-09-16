function [val] = objF(x)
global K
global w
global Phi 
global Y
global N
assert(length(x) == N, 'error! length = %d', length(x))

[~, index] = sort(abs(x),'descend');
weight = ones(size(x));
weight(index(1:K)) = w*ones(K,1);
if size(x,1) == size(Phi,2)
    temp = x;
else
    temp = x';
end
val = sum(weight.*abs(x)) + 1e4 * norm(Phi*temp - Y);
end
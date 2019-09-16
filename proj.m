function [y] = proj(x, B)
y = zeros(size(x));
for iter = 1:size(x,2)
    y(:,iter) = Proj(x(:,iter),B);
end

end

function [val] = Proj(x, B)

val= B*(inv(B'*B)*B'*x);
end
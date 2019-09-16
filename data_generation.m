function [x,phi,y] = data_generation(N,M,K)

index_zeros = randperm(N,N-K);

x =2*rand(N,1)-1;% randn(N,1);%sign(binornd(1, 0.5, N, 1)-0.5);%(0.5+rand(N,1)).*sign(2*rand(N,1)-1); 
x(index_zeros) = zeros(N-K,1);

phi = randn(M,N);


y = phi*x;

% rn = var(y)*0.01;
% y = y + rn*randn(M,1);

end
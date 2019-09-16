function [x] = solve_L1_noisefree(Phi,Y)
N = size(Phi, 2);
x =Phi'*Y;

getX = [eye(N) zeros(N,N)];
getT = [zeros(N,N) eye(N)];

A = [getX-getT;...
    -getX-getT];
b = zeros(2*N,1);

Aeq = Phi*getX;
beq = Y;

f = ones(1,N)*getT;

options = optimoptions('linprog','Algorithm','interior-point','Display','off');
[v,~,exitflag,~] = linprog(f,A,b,Aeq,beq,[],[],options);

if ~isempty(v) && exitflag == 1
    x = getX*v;
end

end
%  solving the following optimization problem
%                 find  x     
%         subject to Phi*x=Y
%                        \sum_{i \in I} |x_i| >= fh
%                        L_{two level}(x) < fg-fh-eps
%

function [x, flag, fval, x1, x2] = find_feasible(x_new, k, w, Phi, Y, fh, fg, eps, fix_idx)
N = length(x_new);
fix_len = length(fix_idx);
nfix_len = N-fix_len;
nfix_idx = setdiff(1:N,fix_idx);
x_fix = x_new;
x_fix(nfix_idx) = zeros(nfix_len,1);

[~, index] = sort(abs(x_new), 'descend');
Wh = eye(N);
Wh(index(k+1:end),:) = [];
Wg = eye(N);
Wg(index(1:k),:) = [];

Wnf = eye(N);
Wnf(:,fix_idx) = [];

fix_val = norm(x_fix,1) - norm(Wh*x_fix,1);

x1=x_new;
x2=x_new;
fval = Inf;
sign_mat = diag(sign(x_new));

getA = [eye(nfix_len) zeros(nfix_len,nfix_len)  zeros(nfix_len,nfix_len)  zeros(nfix_len,nfix_len)  zeros(nfix_len,1)  zeros(nfix_len,1)];
getB = [zeros(nfix_len,nfix_len) eye(nfix_len)  zeros(nfix_len,nfix_len)  zeros(nfix_len,nfix_len)  zeros(nfix_len,1)  zeros(nfix_len,1)];
getD = [zeros(nfix_len,nfix_len) zeros(nfix_len,nfix_len)  zeros(nfix_len,nfix_len)  eye(nfix_len)  zeros(nfix_len,1)  zeros(nfix_len,1)];
getT1 = [zeros(1, nfix_len) zeros(1, nfix_len) zeros(1, nfix_len) zeros(1, nfix_len) 1 0];
getT2 = [zeros(1, nfix_len) zeros(1, nfix_len) zeros(1, nfix_len) zeros(1, nfix_len) 0 1];

A=[getB-getD;...
    -getB-getD;...
    -Wnf'*sign_mat*Wnf*getA;...
    getA-getB-ones(nfix_len,1)*getT1;...
    getB-getA-ones(nfix_len,1)*getT1;...
    -Wh*sign_mat*Wnf*getA+ones(k,1)*getT2;...
    Wg*sign_mat*Wnf*getA-ones(N-k,1)*getT2;...
    -(1-w)*ones(1,k)*Wh*sign_mat*Wnf*getA;...
    ones(1,N)*Wnf*getD-(1-w)*ones(1,k)*Wh*sign_mat*Wnf*getA];

b=[zeros(nfix_len,1);...
    zeros(nfix_len,1);...
    zeros(nfix_len,1);...
    zeros(nfix_len,1);...
    zeros(nfix_len,1);...
    zeros(k,1);...
    zeros(N-k,1);...
    -fh;...
    fg-fh-eps-fix_val;    ];

Aeq = [Phi*Wnf*getB];

beq = [Y-Phi*x_fix];

f=getT1;

options = optimoptions('linprog','Algorithm','interior-point','Display','off');
[v,fval,exitflag,~] = linprog(f,A,b,Aeq,beq,[],[],options);

if ~isempty(v) && exitflag == 1
    x1 = Wnf*getA*v + x_fix;
    x2 = Wnf*getB*v + x_fix;
    if fval < 1e-4
        x_arr = [x1 x2 (x1+x2)/2];
        val_arr = [ get_value(x1,k,w)  get_value(x2,k,w)  get_value((x1+x2)/2,k,w)];
        [minV, idx] = min(val_arr);
        if minV < fg-fh
            x = x_arr(:,idx);
            flag = true;
            return;
        end
    end
end
x= x_new;
flag = false;
end
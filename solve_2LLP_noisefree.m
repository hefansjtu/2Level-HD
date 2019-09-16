function [x, k] = solve_2LLP_noisefree(k, Phi, Y, w, x_init, detect)
[M,N] = size(Phi);
weight_old = ones(N,1);
if nargin < 5     detect =1;  end
if nargin == 3      x_init = initial_point(Phi, Y);  end

% k = norm0(x_init);
f_init = get_value(x_init,k,w);


for iter = 1: 10
    if w~= 0
        w = 10^(-1*iter/2);
    end
    %     get_value(x_init,k,w)
    if detect
        k = support_detection(x_init, iter, M);
    end
    [~,index] = sort(abs(x_init),'descend');
    weight = ones(N,1);
    weight(index(1:k)) = w;
    if weight == weight_old
        break;
    else
        weight_old = weight;
    end
    
    %     end
    
    getX = [eye(N) zeros(N,N)];
    getT = [zeros(N,N) eye(N)];
    
    A = [getX-getT;...
        -getX-getT];
    b = zeros(2*N,1);
    
    Aeq = Phi*getX;
    beq = Y;
    
    f = weight'*getT;
    options = optimoptions('linprog','Algorithm','interior-point','Display','off');
    [v,~,exitflag,~] = linprog(f,A,b,Aeq,beq,[],[],options);
    
    if ~isempty(v) && exitflag == 1
        if norm(x_init-getX*v, inf) < 1e-4 %get_value(x_init,k,w)-get_value(getX*v,k,w)<1e-6
            break;
        else
            x_init = getX*v;
        end
    end
    
end
x = x_init;
% fprintf('Here 2L GD. INPUT: %f, OUTPUT: %f\n', f_init, get_value(x_init,k,w));
end


function [kk] = support_detection(x, iter, M )

[sortX,~] = sort(abs(x));

% tau = 8*norm(x,Inf)/M/iter;
% tau = 1*norm(x,Inf)/M/iter;
tau = 7.0*norm(x,Inf)/M/min(iter,6);
%
dd = diff(sortX);
large = find(dd > tau);
pos = large(find(large,1,'first'));
if isempty(pos)
    kk = length(x)-1;
else
    kk = length(x) - pos;
end

% kk = length(sortX(sortX>1e-3*norm(x,Inf))) - 1 ;

end

function [x] = initial_point(Phi, Y)
opts=[];
opts.maxit=1500;
opts.tol = 1e-1;
opts.weight = ones(size(Phi,2),1);

[x] = yall1_ext(Phi, Y, opts);

end
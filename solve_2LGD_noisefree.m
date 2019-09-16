%  The loss function is
%         min_{x}     L_{two-level}(x)
%         subject to Phi*x=Y
%     based on descent method and YALL1
%
%     Parameters
%     ----------
%     x_init : vector-like, length N, the initial solution given by
%     L1 algorithms
%
%     k : scalar, hyper-parameter in L_{two level}, the estimator of the true K.
%
%     Phi : array-like, shape (M, N), the measure matrix.
%
%     Y: vector-like, length M, the measurements.
%
%     w: scalar, hyper-parameter in L_{two level}, determines the weight.
%
%     detect: bool, decide whether to call iterative support detect.
% 
%     Returns
%     ----------
%     x : vector-like, the sparse signal.
%
%     k : scalar, the estimator of the true K.
%
%     References
%     ----------
%
%     Rule of thumb for tuning paramters:
%     w = 10^{-1*iteration/2};
%     k = norm0(x_init);
%


function [x, k] = solve_2LGD_noisefree(Phi, Y, w, x_init, detect, k)
maxIter = 15;
[M,N] = size(Phi);

if nargin < 5     detect =1;  end
if nargin <= 4  || isempty(x_init)    x_init = initial_point(Phi, Y);  end
if nargin < 6 && ~detect     k =norm0(x_init);  end

opts.tol = 1e-2;
opts.maxit=1500;

for iter = 1: maxIter
    if w~= 0
        w = 10^(-1*iter/2);
    end
    if iter == maxIter
        w = 0;
        opts.tol = 100*eps;
    end
    %     get_value(x_init,k,w)
    if detect
        k = support_detection(x_init, iter, M);
    end
    [~,index] = sort(abs(x_init),'descend');
    weight = ones(N,1);
    weight(index(1:k)) = w;

    opts.weights = weight;
    opts.x0=x_init;

    [x]=yall1(Phi, Y, opts);
    
    if norm(x - x_init, inf) < 1e-5
        if w == 0 && opts.tol == 100*eps
            break;
        else
            w = 0;
            opts.tol = 100*eps;
        end
    else
        x_init = x;
    end

end
x = x_init;
end


function [kk] = support_detection(x, iter, M )

[sortX,~] = sort(abs(x));

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

end

function [x] = initial_point(Phi, Y)
opts=[];
opts.maxit=1500;
opts.tol = 1e-1;
opts.weight = ones(size(Phi,2),1);

[x] = yall1(Phi, Y, opts);

end
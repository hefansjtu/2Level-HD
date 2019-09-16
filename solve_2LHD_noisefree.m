%  The loss function is
%         min_{x}     L_{two-level}(x)
%         subject to Phi*x=Y
%     based on hill detouring method
%
%     Parameters
%     ----------
%     x_init : vector-like, length N, the initial solution given by
%     local algorithms
%
%     k : scalar, hyper-parameter in L_{two level}, the estimator of the true K.
%
%     Phi : array-like, shape (M, N), the measure matrix.
%
%     Y: vector-like, length M, the measurements.
%
%     w: scalar, hyper-parameter in L_{two level}, determines the weight.
%
%     isCS: bool, true for CS problem else for optimization.
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
%     w = 0;
%     k = norm0(x_init);
%
function [x, k] = solve_2LHD_noisefree(x_init, Phi, Y, k, w, isCS, fixK)
if nargin == 6   fixK = ~isCS; end
if nargin <= 5   isCS = false; fixK = true; end
if nargin <= 4   w = 0; end
if nargin <= 3   k = norm0(x_init); end
if nargin <= 2   fprintf('error! require more inputs\n'); return; end

eps = 1e-6;
maxIter = 10;
maxII = 20;
N = length(x_init);
x=x_init;
k_init = k;

for iter = 1: maxIter
    inaccra_flag = 0;
    fg = norm(x_init, 1);
    f0 = get_value(x_init, k, w);
    %     if f0 <= eps*10
    %         return;
    %     end
    fh = fg-f0;
    directions = randn(N,2*N);
    [newpoints, steps] = gamma_extension(x_init, fh, k, w, directions);
    numofpoints = length(steps(steps<Inf));
    if numofpoints == 0
        continue;
    end
    
    for ii = 1: maxII
        f0 = get_value(x_init, k, w);
        fh = fg-f0;
        [minS, idx] = min(steps);
        if minS == Inf
            break;
        end
        x_new = newpoints(:,idx);
        
        [x_opt, flag, fval, x1, x2] = find_feasible(x_new, k, w, Phi, Y, fh, fg, eps, []);
        if flag
            if ~isCS
                if get_value(x_opt, k, w) < get_value(x_init,k,w)
                    x_init = x_opt;
                end
                [x_lp, k_lp] = solve_2LLP_noisefree(k, Phi, Y, w, x_opt, fixK);
                if norm(x_lp - x_init, inf ) > 1e-4 && get_value(x_lp, k, w) < get_value(x_init,k,w)
                    x_init = x_lp;
                    k = k_lp;
                    break;
                else
                    inaccra_flag = inaccra_flag + 1;
                    if inaccra_flag >= 4
                        ii = maxII;
                        break;
                    end
                end
            else
                [x_gd, k_gd] = solve_2LGD_noisefree(Phi, Y, w, x_opt, fixK, k); %K is given: 0 ; else: 1
                if norm(x_gd - x_init, inf ) > 1e-4 &&  norm0(x_gd) <= norm0(x_init)
                    x_init = x_gd;
                    k = k_gd;
                    break;
                else
                    [x_lp, k_lp] = solve_2LLP_noisefree(k, Phi, Y, w, x_opt, fixK);
                    if norm(x_lp - x_init, inf ) > 1e-4 &&  norm0(x_lp) <= norm0(x_init)
                        x_init = x_lp;
                        k = k_lp;
                        break;
                    else
                        inaccra_flag = inaccra_flag + 1;
                        if inaccra_flag >= 4
                            ii = maxII;
                            break;
                        end
                    end
                end
            end
        end
        Wh = (1-w)*ones(N,1);
        [~,index] = sort(abs(x_new), 'descend');
        Wh(index(k+1:end)) = zeros(N-k,1);
        dd = randn(N,1);
        
        ddd = dd - Wh*(Wh'*dd)/(norm(Wh))^2;
        [np, s] = bisection(x_new, ddd, k, w);
        newpoints(:,idx) = np;
        if norm(s*ddd,inf)<1e-2*norm(x_init,inf)
            steps(idx)=Inf;
        else
            steps(idx) = s;
        end
    end
    if ii == maxII
        if ~isCS
            break;
        else
            if k < floor(size(Y,1)/2)
                k = floor(size(Y,1)/2);
            elseif k == floor(size(Y,1)/2)
                break;
            else
                k = k - 1;
            end
            if  k/k_init <= 0.5
                break;
            end
        end
    end
    
    x=x_init;
end
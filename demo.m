clc; clear; close all;

iterNum = 50;
N = 128;
K_range = 25;
M_range = 50;%55-100
[KK, MM] = meshgrid(K_range,M_range);

isPhase = false;

result_x1 = zeros(size(KK));
result_x_GD = zeros(size(KK));
result_x_HD = zeros(size(KK));
result_x_ISD = zeros(size(KK));

for iter_m = 1: length(M_range)
    fprintf('M = %d\n', M_range(iter_m));
    for iter_k = 1: length(K_range)
        K = K_range(iter_k);    
        M = M_range(iter_m);
        %         fprintf('M = %d, K = %d...\n', M,K );
        if isPhase
            if M > K + 50 || K == 0  ||  iter_m>2 && result_x1(iter_m-2,iter_k) == iterNum && result_x1(iter_m-1,iter_k) == iterNum
                result_x1(iter_m, iter_k) = iterNum;
                result_x_GD(iter_m, iter_k) = iterNum;
                result_x_HD(iter_m, iter_k) = iterNum;
                result_x_ISD(iter_m, iter_k) = iterNum;30
                continue;
            elseif K > M
                continue;
            elseif iter_k > 2 && result_x_HD(iter_m,iter_k-2) == 0 && result_x_HD(iter_m,iter_k-1) == 0
                continue;
            end
        end
        for iter = 1: iterNum
            [X, Phi, Y] = data_generation(N,M,K);
            %% L1
            opts=[];
            opts.maxit=1500;
            opts.tol = 100*eps;
            opts.weight = ones(N,1);
            tic;
            [x1]=yall1_ext(Phi, Y, opts);
%             [x1]=solve_L1_noisefree(Phi,Y);
            if SNR(x1, X)>50
                result_x1(iter_m, iter_k) = result_x1(iter_m, iter_k) + 1;
            end
            
            %% 2L-GD
            w = 0.1;
            k = floor(norm0(x1));
            [x_GD] = solve_2LGD_noisefree(Phi, Y, w);
            if SNR(x_GD, X)>50
                result_x_GD(iter_m, iter_k) = result_x_GD(iter_m, iter_k) + 1;
            end
            
            %% ISD
            opts=[];
            opts.tol=100*eps;
            opts.maxit=10;
            opts.sigma = 0;
            [x_ISD,~] = Threshold_ISD_1D(Phi,Y,opts);
            if SNR(x_ISD, X)>50
                result_x_ISD(iter_m, iter_k) = result_x_ISD(iter_m, iter_k) + 1;
            end
            
            %% 2L-HD
            w = 0.00;
            k = floor((0.9*norm0(x_GD)+k)/2);
            [x_HD] = solve_2LHD_noisefree(x_GD, Phi, Y, k, w, 1);
            if SNR(x_HD, X)>50
                result_x_HD(iter_m, iter_k) = result_x_HD(iter_m, iter_k) + 1;
            end
        end
        
    end
end
result_x1 = result_x1./iterNum;
result_x_GD = result_x_GD./iterNum;
result_x_HD = result_x_HD./iterNum;
result_x_ISD = result_x_ISD./iterNum;

% filename = '0915N512M120-uniform.mat';
% save(filename);

function [val] = SNR(xx, X)
val = 10*log10(norm(X)^2/norm(X-xx)^2);
end

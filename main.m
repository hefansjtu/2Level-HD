clc; clear; close all;

iterNum = 50;
N = 128;
K = 50;
M = 80;
lalala = 0;
hahaha = 0;
wuwuwu = 0;
SNR_X1=zeros(1,iterNum);
SNR_X_GD=zeros(1,iterNum);
SNR_X_HD=zeros(1,iterNum);
SNR_X_ISD=zeros(1,iterNum);
counter_X1 = 0;
counter_X_GD = 0;
counter_X_HD = 0;
counter_X_ISD = 0;
norm0_X1=zeros(1,iterNum);
norm0_X_GD=zeros(1,iterNum);
norm0_X_HD=zeros(1,iterNum);
norm0_X_ISD=zeros(1,iterNum);
time_X1=zeros(1,iterNum);
time_X_GD=zeros(1,iterNum);
time_X_HD=zeros(1,iterNum);
time_X_ISD=zeros(1,iterNum);
value_X1=zeros(1,iterNum);
value_X_GD=zeros(1,iterNum);
value_X_HD=zeros(1,iterNum);
value_X_ISD=zeros(1,iterNum);
for iter = 1: iterNum
    
    fprintf('iter = %d...\n', iter);
    [X, Phi, Y] = data_generation(N,M,K);
    
    %% L1
    %     tic
    %     [x1] = solve_L1_noisefree(N,Phi,Y);
%     fprintf('begin L1...\n');
    opts=[];
    opts.maxit=1500;
    opts.tol = 100*eps;
    opts.weight = ones(N,1);
    tic;
    [x1] = yall1(Phi, Y, opts);
    time_X1(iter)  = toc;
    SNR_X1(iter) = SNR(x1, X);
    if SNR_X1(iter)  == 0
    end
    if SNR_X1(iter)>50
        counter_X1 = counter_X1 + 1;
    end
    norm0_X1(iter) = norm0(x1);
    %% 2L-GD
    w = 0.1;
    k = floor(norm0_X1(iter));
    tic
    [x_GD, k] = solve_2LGD_noisefree(k, Phi, Y, w);
    time_X_GD(iter)  = toc;
    SNR_X_GD(iter) = SNR(x_GD, X);
    if SNR_X_GD(iter)>50
        counter_X_GD = counter_X_GD + 1;
    else
        
    end
    norm0_X_GD(iter) = norm0(x_GD);
    
    %% ISD
    opts=[];
    %     opts.alpha=20;
    opts.tol=100*eps;
    opts.maxit=10;
    opts.sigma = 0;
    tic
    [x_ISD,~] = Threshold_ISD_1D(Phi,Y,opts);
    time_X_ISD(iter)  = toc;
    SNR_X_ISD(iter) = SNR(x_ISD, X);
    if SNR_X_ISD(iter)>50
        counter_X_ISD = counter_X_ISD + 1;
        if SNR_X_GD(iter) < 50
            wuwuwu = wuwuwu + 1;
        end
    elseif SNR_X_GD(iter) > 50
        lalala=lalala + 1;
    end
    norm0_X_ISD(iter) = norm0(x_ISD);
%     continue;
    
    %% 2L-HD
    w = 0.00;
    k = floor((0.9*norm0_X_GD(iter)+k)/2);
    tic
    [x_HD, ~] = solve_2LHD_noisefree(x_GD, k, Phi, Y, w);
    time_X_HD(iter)  = toc;
    SNR_X_HD(iter) = SNR(x_HD, X);
    if SNR_X_HD(iter)>50
        counter_X_HD = counter_X_HD + 1;
        if SNR_X_GD(iter) < 50 && SNR_X_ISD(iter) < 50
            hahaha = hahaha+1;
        end
    elseif SNR_X_GD(iter) > 50
        break;
    end
    norm0_X_HD(iter) = norm0(x_HD);
    
    
    value_X1(iter) = get_value(x1, k, w);    
    value_X_ISD(iter) = get_value(x_ISD, k, w);    
    value_X_GD(iter) = get_value(x_GD, k, w);    
    value_X_HD(iter) = get_value(x_HD, k, w);    
end

function [val] = norm0(x)
val = length(x(abs(x)>1e-3*norm(x,inf)));
end

function [val] = SNR(xx, X)
val = 10*log10(norm(X)^2/norm(X-xx)^2);
end


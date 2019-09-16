clc; clear; close all;
global K
global w
global Phi
global Y
global optimum
global N
iterNum = 10;
N = 128;
K_range = [8; 32; 64;];
M_range = [25; 60; 95;];
w = 0;
optimum = 0;

result_x_ReL1 = zeros(iterNum, length(K_range));
result_x_GA = zeros(iterNum, length(K_range));
result_x_ABCD = zeros(iterNum, length(K_range));
result_x_HD = zeros(iterNum, length(K_range));
time_x_ReL1 = zeros(iterNum, length(K_range));
time_x_GA = zeros(iterNum, length(K_range));
time_x_ABCD = zeros(iterNum, length(K_range));
time_x_HD = zeros(iterNum, length(K_range));

low = -1*ones(N,1);
up = 1*ones(N,1);

for iter_k = 1: length(K_range)
    K = K_range(iter_k);
    M = M_range(iter_k);
    
    for iter_m = 1: iterNum
        [X, Phi, Y] = data_generation(N,M,K);
        %% GA
        tic;
        [x_GA] = ga(@objF, N, [], [], [], [], low, up);
        result_x_GA(iter_m, iter_k) = objF(x_GA);
        time_x_GA(iter_m, iter_k) = toc;
        
        %% ABCD
        settings.NumFixed = 1;
        settings.SwitchWay  = 1;
        settings.SelectWay = 1;
        settings.SwitchWay2 = 0;
        settings.SwitchEps = 1e-2;
        settings.SwitchLength = 3;
        settings.SqpCheck = 0;
        setings.TMax = 20;
        %settings.MaxFeval = EvalFES;
        settings.ConvergeCheck = 1;
        settings.TimeCheck = 1;
        settings.NumFixed2 = 2;
        settings.SelectWay2 = 0;
        %---------------opts for DIRECT----------------
        opts.stop = 4;
        opts.epsilon = 1e-4;
        opts.converge = 5;
        opts.detail= 0;
        opts.accuracy = 1e-4;
        [ f_best, x_best, num_feval, RecordData,  running_t, exit_flag ] = ABCD_test( @objF,  low, up, settings, opts );
        result_x_ABCD(iter_m, iter_k) = f_best;
        time_x_ABCD(iter_m, iter_k) = running_t;
% %         
        %% 2L-GD
        fixK = true;
        tic
        [x_GD, ~] = solve_2LGD_noisefree(Phi, Y, w, [], ~fixK,  K);
        result_x_ReL1(iter_m, iter_k) = objF(x_GD);
        time_x_ReL1(iter_m, iter_k) = toc;
        %% 2L-HD
        isCS = false;
        tic
        [x_HD] = solve_2LHD_noisefree(x_GD, Phi, Y, K, w, isCS, fixK);
        result_x_HD(iter_m, iter_k) = objF(x_HD);
        time_x_HD(iter_m, iter_k) = toc;
        
    end
    
end

% filename = '0729N128-optimization.mat';
% save(filename);

clc;close all;clear;
state = 0;
rng('default');
addpath(genpath('src_HtGaP_DNM/'));
if 1
%% synthetic data
N = 65;
T = 6;
A = zeros(N, N, T);
B = cell(1, T);

bd = [0 10 30 65];
C = numel(bd) -1;
t = 1;
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
    A(dex, dex, t) = 1;
    A(:, :, t) = A(:, :, t) - diag(diag(A(:, :, t)));
    B{t} = sparse(A(:, :, t));
end
subplot(T, 4, (t-1)* 4 + 1);imagesc(A(:,:,t));colormap('jet');

for t = 2:T
    if t == 2
        bd = [0 15 25 35 65];
    elseif t == 3
        bd = [0 18 30 44 65];
    else
        bd = bd + [0 5 5 0];
    end
        
    C = numel(bd) -1;
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        A(dex, dex, t) = 1;
        A(:, :, t) = A(:, :, t) - diag(diag(A(:, :, t)));
        B{t} = sparse(A(:, :, t));
    end
    subplot(T, 4, (t-1)* 4 + 1);imagesc(A(:,:,t));
    
    if t == 2
        bd = [0 15 35 65];
    elseif t == 3
        bd = [0 20 40 65];
    end
end

options.burnin = 1000;
options.mcmcsamps = 100;
options.display = 'TRUE';

%% create training test data
Acore = B;
TrainRatio = 1; 
numslices = T;
idx_train = cell(numslices, 1); idx_test = cell(numslices, 1); BTrain_Mask = cell(numslices, 1);
figure(1);
for t = 1:numslices
[idx_train{t}, idx_test{t}] = Create_Mask_big_network(Acore{t}, TrainRatio);
% subplot(numslices, 4, (t-1)* 4 + 2);
% imagesc(~BTrain_Mask{t});
end
options.TrainRatio = TrainRatio;
options.idx_train = idx_train;
options.idx_test = idx_test;
options.Mex = 1;
options.eval = 1;
options.thinningfunc = 2;
options.display=1;
%%
% [result] = call_DRIFT(B, options);
numrep = 1; result = cell(numrep, 1);
options.samp_c_n = 0;
 
for rep = 1:numrep
%result{rep} = DGPPF_batch_Gibbs(B, options);
% result{rep} = HGPDR_batch_Gibbs(B, options);
result{rep} = HtGaP_DNM_batchGibbs(B, options);
end
end
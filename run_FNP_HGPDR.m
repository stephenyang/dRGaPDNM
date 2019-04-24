
clear;clc;addpath(genpath('src_HGPDR/'));
prompt = 'Type the seed number:\n';
seed = input(prompt);
filename = 'dblp_324.mat';
load(['dataset/',filename]);

options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 'False';

%% create training test data
% T = size(B, 3);
% options2.K = 15;
% Acore = cell(1, T);
% for t = 1:T
%     Acore{t} = B(:, :, t);
% end
numslices = numel(Acore);
TrainRatio = 0.8; 
options.TrainRatio = TrainRatio;
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);
for t = 1:numslices
[idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
end
options.idx_train = idx_train;
options.idx_test = idx_test;

TrainRatio = 1; options2.TrainRatio = TrainRatio;
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);
for t = 1:numslices
[idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
end
options2.idx_train = idx_train;
options2.idx_test = idx_test;
options2.burnin = 2000;
options2.mcmcsamps = 1000;
options2.display = 0;
%%
numrep = 5;
result = cell(numrep,1);
prompt = 'Type the rep number:\n';
rep = input(prompt);
    if 0
    %% missing link prediction
    options.eval = 1;
    result = HGPDR_batch_Gibbs(Acore(1:numslices), options);
    save(['result2/result_MLP_nips_HGPDR_rng',num2str(rep),'.mat'], 'options','result');
    end
    %% future network prediction
    if 1
    options2.eval = 0;
    for numslices = 3:numel(Acore)
        N = size(Acore{1}, 1);
        fprintf(['HGPDR_T:',num2str(numslices),', N:',num2str(N),'.\n']);
        samp = HGPDR_batch_Gibbs(Acore(1:(numslices-1)), options2);
        [pred] = increament_phi_HGPDR(Acore(numslices), samp);
        save(['result/HGPDR/result_FNP_pred',num2str(numslices),'_HGPDR_rng',num2str(rep),'.mat'], 'options','pred');
    end
    end

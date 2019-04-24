clear;clc;addpath(genpath('src_HtGaP_DNM/'));
prompt = 'Type the seed number:\n';
seed = input(prompt);
% prompt = 'Type the rep number:\n';
% rep = input(prompt);
s = RandStream('mt19937ar','Seed',seed); 
RandStream.setGlobalStream(s);
gs = RandStream.getGlobalStream();
fprintf('nips274_Random seed: %d\n', gs.Seed);
filename = 'dblp_324.mat';
load(['dataset/',filename]);
% T = size(B, 3);
% options2.K = 15;
% Acore = cell(1, T);
% for t = 1:T
%     Acore{t} = B(:, :, t);
% end
% Acore = yearly_collabs;
T = numel(Acore);
numslices = numel(Acore);
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
options2.display = 'False';
options2.Mex = 1;
%%
numrep = 10;
result = cell(numrep,1);
for rep = 1:numrep
    if 0
    %% missing link prediction
    options.eval = 1;
    result = DGPPF_batch_Gibbs(Acore(1:numslices), options);
    save(['result2/result_nips5K_DGPPF_rng',num2str(rep),'.mat'], 'options','result');
    end
    %% future network prediction
    options2.eval = 0;
    for numslices = 3:numel(Acore)
        N = size(Acore{1}, 1);
        fprintf(['DGPPF_T:',num2str(numslices),', N:',num2str(N),'.\n']);
        samp = DGPPF_batch_Gibbs(Acore(1:(numslices-1)), options2);
        [pred] = increament_DGPPF(Acore(numslices), samp);
        save(['result/result_FNP_pred',num2str(numslices),'_DGPPF_rng',num2str(rep),'.mat'],'pred','-v7.3');
    end
end
%%
clear;clc;addpath(genpath('src_HtGaP_DNM/'));
filename = 'dblp_324.mat';
load(['dataset/',filename]);
 
% T = numel(yearly_collabs);
% Acore = yearly_collabs;
%% create training test data
TrainRatio = 0.8;   
numslices = numel(Acore);
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);
options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 'False';
options.Mex = 1;
options.eval = 1;
run = 10;
seed = [1234];
fprintf('T: %d, N: %d.\n', [numel(Acore), size(Acore{1},1)]);
for rep = 1:run

%%
s = RandStream('mt19937ar','Seed', seed(rep)); % Good demo seed
RandStream.setGlobalStream(s);
gs = RandStream.getGlobalStream();
fprintf('Random seed: %d\n', gs.Seed);
    
for t = 1:numslices
[idx_train{t}, idx_test{t},~] = Create_Mask_network(full(Acore{t}), TrainRatio);
end

options.TrainRatio = TrainRatio;
options.idx_train = idx_train;  
options.idx_test = idx_test;
%%
% fprintf('NIPS_274_DGPPF %02d.\n', [rep]);
tic;
result = DGPPF_batch_Gibbs(Acore(1:numslices), options);
result.timecost = toc

save(['result/result_dblp324_DGPPF_rng',num2str(gs.Seed),'.mat'], 'options','result', '-v7.3');
end

clear;clc
prompt = 'Type the seed number:\n';
seed = input(prompt);
s = RandStream('mt19937ar','Seed',seed); 
RandStream.setGlobalStream(s);
gs = RandStream.getGlobalStream();
filename = 'dblp_324.mat';
load(['dataset/',filename]);
addpath(genpath('src_HGPDR/'));
options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 0;

%% create training test data
% T = size(B, 3);
% Acore = cell(1,T);
% for t = 1:T
%     Acore{t} = B(:, :, t);
% end
numslices = numel(Acore);
TrainRatio = 0.8; 
options.TrainRatio = TrainRatio;
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);
% for t = 1:numslices
% [idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
% end
% options.idx_train = idx_train;
% options.idx_test = idx_test;

% TrainRatio = 1; options2.TrainRatio = TrainRatio;
% idx_train = cell(numslices, 1); 
% idx_test = cell(numslices, 1); 
% BTrain_Mask = cell(numslices, 1);
% for t = 1:numslices
% [idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
% end
options2.idx_train = idx_train;
options2.idx_test = idx_test;
options2.burnin = 2000;
options2.mcmcsamps = 1000;
options2.display =0;
%%
numrep = 5;
result = cell(numrep,1);
seed = [4604 8365 6969 1234 5432 7575 8848 9999 1010];
for rep = 1:numrep
    if 1
        s = RandStream('mt19937ar','Seed', seed(rep)); % Good demo seed
        RandStream.setGlobalStream(s);
        gs = RandStream.getGlobalStream();
        fprintf('Random seed: %d\n', gs.Seed);
        
        for t = 1:numslices
            [idx_train{t}, idx_test{t},~] = Create_Mask_network(full(Acore{t}), TrainRatio);
        end
    %% missing link prediction
    options.eval = 1;
    options.idx_train = idx_train;
    options.idx_test = idx_test;
    options.TrainRatio = TrainRatio;
    result = HGPDR_batch_Gibbs(Acore(1:numslices), options);
    save(['result/result_HGPDR_rng',num2str(rep),'.mat'], 'options','result');
    end
    %% future network prediction
    if 0
    options2.eval = 0;
    for numslices = 3:numel(Acore)
        N = size(Acore{1}, 1);
        fprintf(['HGPDR_T:',num2str(numslices),', N:',num2str(N),'.\n']);
        samp = HGPDR_batch_Gibbs(Acore(1:(numslices-1)), options2);
        [pred] = increament_phi(Acore(numslices), samp);
        save(['result/result_pred',num2str(numslices),'_HGPDR_rng',num2str(rep),'.mat'], 'options','pred');
    end
    end
end
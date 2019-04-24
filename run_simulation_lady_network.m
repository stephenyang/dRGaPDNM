

clear;clc;close all;
addpath(genpath('src_HGPDR/'));
addpath(genpath('src_HtGaP_DNM/'));

for seed = [1234 5685 9999 4321 8848];
s = RandStream('mt19937ar','Seed',seed); 
state = seed;
RandStream.setGlobalStream(s);

% LADY_network;
% block_structured_data2;
%%
% Acore = cell(1,T);
% for t = 1:T
%     Acore{t} = Adj(:, :, t);
% end
%%
LADY_network2;
options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 1;
options.Mex = 1;
options.eval = 0;
options.thinningfunc = 1;
%% create training test data
Acore = cell(1,T);
for t = 1:T
    Acore{t} = EPr_t(:, :, t);
end
TrainRatio = 1; 
numslices = T;
idx_train = cell(numslices, 1); idx_test = cell(numslices, 1); BTrain_Mask = cell(numslices, 1);
figure(1);
for t = 1:numslices
[idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
% subplot(numslices, 4, (t-1)* 4 + 2);
% imagesc(~BTrain_Mask{t});
end
options.TrainRatio = TrainRatio;
options.idx_train = idx_train;
options.idx_test = idx_test;

% numrep = 3; 

% for rep = 1:numrep
result = DGPPF_batch_Gibbs(Acore, options);
% save(['result/result_DGPPF_', date, '_track',num2str(state),'.mat'], 'options','result');
% end
% 
% for rep = 1:numrep
result = HGPDR_batch_Gibbs(Acore, options);
% save(['result/result_HGPDR_', date, '_track',num2str(state),'.mat'], 'options','result');
% end
options.display = 1;
options.samp_c_n = 0;
% result = cell(1,3);
% for rep = 1:numrep
% fprintf('rand seed:%d',[state]);
result = HtGaP_DNM_batchGibbs(Acore, options);
% save(['result/result_newblock_HtGaP_', date, '_track',num2str(state),'.mat'], 'options','result');

% end
end
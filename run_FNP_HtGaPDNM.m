clear;clc;addpath(genpath('src_HtGaP_DNM/'));
prompt = 'Type the seed number:\n';
seed = input(prompt);
prompt = 'Type the rep number:\n';
rep = input(prompt);
s = RandStream('mt19937ar','Seed',seed); 
RandStream.setGlobalStream(s);
gs = RandStream.getGlobalStream();
fprintf('dblp324_Random seed: %d\n', gs.Seed);
filename = 'dblp_324.mat';
load(['dataset/',filename]);

% Acore = yearly_collabs;
T = numel(Acore);
numslices = numel(Acore);
TrainRatio = 1; options2.TrainRatio = TrainRatio;
idx_train = cell(T, 1); 
idx_test = cell(T, 1); 
BTrain_Mask = cell(T, 1);
for t = 1:T
[idx_train{t}, idx_test{t}, BTrain_Mask{t}] = Create_Mask_network(Acore{t}, TrainRatio);
end
options2.idx_train = idx_train;
options2.idx_test = idx_test;
options2.burnin = 2000;
options2.mcmcsamps = 1000;
options2.display =0;
options2.Mex = 1;
options2.thinningfunc =1;
options2.samp_c_n = 0;
%%
numrep = 5;
result = cell(numrep,1);
for rep = 1:numrep
    if 0
    %% missing link prediction
    options.eval = 1;
    result = HGPDR_batch_Gibbs(Acore(1:numslices), options);
    save(['result/result_MLP_enron_HGPDR_rng',num2str(rep),'.mat'], 'options','result');
    end
    %% future network prediction
    if 1
        options2.eval = 0;
        for numslices = 3:T;
            N = size(Acore{1}, 1);
            fprintf(['HtGaPDR_T:',num2str(numslices),', N:',num2str(N),'.\n']);
            % samp = HtGaP_DNM_large_network(Acore(1:(numslices-1)), options2);
            samp = HtGaP_DNM_batchGibbs(Acore(1:(numslices-1)), options2);
            [pred] = increament_phi(Acore(numslices), samp);
            fprintf('pred_slice:%02d, aucroc:%f, aucpr:%f.\n',[numslices, pred.aucroc, pred.aucprec]);
            save(['result/result_pred',num2str(numslices),'_HtGaPDNM_rng',num2str(rep),'.mat'],'pred','-v7.3');
        end
    end
end
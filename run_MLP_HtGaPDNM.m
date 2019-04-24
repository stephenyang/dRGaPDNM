
%%
clear;clc;addpath(genpath('src_HtGaP_DNM/'));
filename = 'dblp_324.mat';
load(['dataset/',filename]);
%% create training test data
TrainRatio = 0.8;   
numslices = numel(Acore);
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);
options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 0;
run = 10;
% seed = [4604 8365 6969 1234 5432  3691 7575 8848 9999 1010];
fprintf('T: %d, N: %d.\n', [numel(Acore), size(Acore{1},1)]);
for rep = 1:run

%%
% s = RandStream('mt19937ar','Seed', seed(rep)); % Good demo seed
% RandStream.setGlobalStream(s);
% gs = RandStream.getGlobalStream();
% fprintf('Random seed: %d\n', gs.Seed);
    
for t = 1:numslices
[idx_train{t}, idx_test{t},~] = Create_Mask_network(full(Acore{t}), TrainRatio);
end

options.TrainRatio = TrainRatio;
options.idx_train = idx_train;  
options.idx_test = idx_test;
%%
fprintf('dblp324_HtGaPDNM %02d.\n', [rep]);
tic;
result = HtGaP_DNM_large_network(Acore, options);
%% quantify results
if isfield(options,'TrainRatio')
    if options.TrainRatio == 1
        idx_test = idx_train;        
    end
        T = numslices;
        rate_test = cell(T,1);
        links = cell(T, 1);
        coll_rate = [];
        coll_links = [];
        for t = 1:T
            Acore{t} = full(Acore{t});
            rate_test{t} = result.ProbAve{t}(idx_test{t});
            links{t} = Acore{t}(idx_test{t});
            coll_rate = [coll_rate;rate_test{t}];
            coll_links = [coll_links;links{t}];
        end
        % coll_rate = double(coll_rate > 1e-3);
%         figure(234);
%         subplot(1,2,1);
        [XX, YY, TT, AUCroc] = perfcurve(coll_links, coll_rate, 1);
%         plot(XX,YY);
%         axis([0 1 0 1]), grid on, xlabel('FPR'), ylabel('TPR'), hold on;
%         x = [0:0.1:1];plot(x,x,'b--'), hold off; title(['AUCroc = ', num2str(AUCroc)])
        
%         subplot(1,2,2);
        [prec, tpr, fpr, thresh] = prec_rec(coll_rate, coll_links,  'numThresh',3000);
%         plot([0; tpr], [1 ; prec]); % add pseudo point to complete curve
%         xlabel('recall');
%         ylabel('precision');
%         title('precision-recall graph');
        AUCpr = trapz([0;tpr],[1;prec]);
        F1= max(2*tpr.*prec./(tpr+prec));
%         title(['AUCpr = ', num2str(AUCpr), '; F1 = ', num2str(F1)])
        fprintf('aucroc: %f, aucprec: %f.\n', [AUCroc, AUCpr]);
        result.aucroc = AUCroc;
        result.aucprec = AUCpr;
end
%%
result.timecost = toc

save(['result/result_HtGaPDNM_rng',num2str(gs.Seed),'.mat'], 'options','result', '-v7.3');
end
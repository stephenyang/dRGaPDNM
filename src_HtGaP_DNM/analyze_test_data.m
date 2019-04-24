
%%
clear;clc
filename = 'dataset_highschool180.mat';

load(['dataset/',filename]);

options.burnin = 2000;
options.mcmcsamps = 1000;
options.display = 'False';
T = size(B, 3);
Acore = cell(1,T);
for t = 1:T
    Acore{t} = B(:, :, t);
end
%% create training test data
TrainRatio = 0.8;   
numslices = numel(Acore);
idx_train = cell(numslices, 1); 
idx_test = cell(numslices, 1); 
BTrain_Mask = cell(numslices, 1);


seed = [8365 6969 12345 54321 4604 3691 7575 8848 9999 1010 12345 54321 1212];
run = numel(seed);
for rep = 1:run

%%
s = RandStream('mt19937ar','Seed', seed(rep)); % Good demo seed
RandStream.setGlobalStream(s);
gs = RandStream.getGlobalStream();
fprintf('Random seed: %d\n', gs.Seed);
    
for t = 1:numslices
[idx_train{t}, idx_test{t}] = Create_Mask_big_network(full(Acore{t}), TrainRatio);
end

coll_links = [];
for t = 1:T
    links{t} = Acore{t}(idx_test{t});
    coll_links = [coll_links;links{t}];
end
nnz(coll_links)

end
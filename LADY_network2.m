% clear;
close all;
N = 50;
R = 5;
EPr = zeros(N, N, R);
diagdex = sparse(1:N, 1:N, true);
%%
Rind = 1;
EPr_t = zeros(N) + 0.05 * rand(N);
a = 0.7; b = 1;
C = 5;bd = [0 12 20 32 40 50];
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc == 2 || cc==4
    EPr_t(dex, dex) = rand(numel(dex)) * 0.05;
else
    EPr_t(dex, dex) = a + (b-a).* rand(numel(dex)) ;
end
end
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
if 0
figure;imagesc(EPr_t);colorbar;% colormap('gray');
end
EPr(:,:, Rind) = EPr_t;
%%
Rind = 2;
EPr_t = zeros(N) + 0.05 * rand(N);
C = 6;bd = [0 15 24 34 41 48 50];
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc == 2 || cc==4
    EPr_t(dex, dex) = rand(numel(dex)) * 0.02;
else
    EPr_t(dex, dex) = a + (b-a).* rand(numel(dex)) ;
end
end
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
if 0
figure;imagesc(EPr_t);%colormap('gray');
end
EPr(:,:, Rind) = EPr_t;
%%
Rind = 3;
EPr_t = zeros(N) + 0.05 * rand(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);
C = 5;bd = [0 10 20 28 36 50];
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
% if cc == 2 || cc==4
%     EPr_t(dex, dex) = rand(numel(dex)) * 0.02;
% else
    EPr_t(dex, dex) = a + (b-a).* rand(numel(dex)) ;
% end
end
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
if 0
figure;imagesc(EPr_t);%colormap('gray');
end
EPr(:,:, Rind) = EPr_t;

%%
Rind =4;
EPr_t = zeros(N) + 0.05 * rand(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);

C = 5;bd = [0 10 20 28 36 50];
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
 if cc == 2%  || cc==4
%     EPr_t(dex, dex) = rand(numel(dex)) * 0.02;
 else
    EPr_t(dex, dex) = a + (b-a).* rand(numel(dex)) ;
 end
end
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
if 0
figure;imagesc(EPr_t);%colormap('gray');
end
EPr(:,:, Rind) = EPr_t;

%%
Rind = 5;
EPr_t = zeros(N) + 0.05 * rand(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);
C = 4;bd = [0 20 28 36 50];
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
% if cc == 2 || cc==4
%     EPr_t(dex, dex) = rand(numel(dex)) * 0.02;
% else
    EPr_t(dex, dex) = a + (b-a).* rand(numel(dex)) ;
% end
end
remove_id = [3 5 7 8 10 11 13 15 17 18 20 23 24 28];
for id = 1:numel(remove_id)
    EPr_t(remove_id(id), :) = 0;
    EPr_t(:, remove_id(id)) = 0;
end
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
if 0
figure;imagesc(EPr_t);%colormap('gray');
end
EPr(:,:, Rind) = EPr_t;

%%
T = 70;
% regime_seq_indices = zeros(1,T);
% regime_seq_indices(1:4)     = 1;
% regime_seq_indices(5:8)   = 6;
% regime_seq_indices(9:11)   = 3;
% regime_seq_indices(12:14) = 1;
% regime_seq_indices(15:19) = 2;
% regime_seq_indices(20:25) = 4;
% regime_seq_indices(26:30) = 2;
% regime_seq_indices(31:33) = 1;
% regime_seq_indices(34:38) = 3;
% regime_seq_indices(39:42) = 1;
% regime_seq_indices(47:50) = 4;
% regime_seq_indices(43:46) = 5;

regime_seq_indices(1:10) = 1;
regime_seq_indices(11:20) = 5;
regime_seq_indices(21:30) = 2;
regime_seq_indices(31:40) = 3;
regime_seq_indices(41:55) = 2;
regime_seq_indices(56:65) = 4;
regime_seq_indices(65:70) = 1;
% figure(2);plot(1:T, regime_seq_indices, 'bo');
% ylim([0 R+1])
% set(gca, 'YGrid', 'on', 'XGrid', 'off');

%%
EPr_t = zeros(N, N, T);
for t = 1:T
    Adj_t = rand(N) < EPr(:,:, regime_seq_indices(t));
    Adj_t = triu(Adj_t, 1) + triu(Adj_t, 1)';
    EPr_t(:,:,t) = Adj_t;
end

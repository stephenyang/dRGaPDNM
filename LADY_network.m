% clear;
close all;
N = 30;
R = 6;
EPr = zeros(N, N, R);
diagdex = sparse(1:N, 1:N, true);
%%
Rind = 1;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);
a = 0.7; b = 1;
EPr_t(1:10, 1:10) = a + (b-a).*rand(10, 10);
EPr_t(11:20, 11:20) = a + (b-a).*rand(10, 10);
EPr_t(21:30, 21:30) = a + (b-a).*rand(10, 10);
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;
%%
Rind = 2;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);
a = 0.9; b = 1;
subind = 1:7;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 8:15;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 16:23;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 24:30;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;
%%
Rind = 3;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);

a = 0.2; b = 0.7;
subind = 1:20;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

a = 0.7; b = 1;
subind = 1:10;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 11:20;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 21:30;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;
%%
Rind = 3;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);

a = 0.2; b = 0.7;
subind = 1:20;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

a = 0.7; b = 1;
subind = 1:10;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 11:20;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 21:30;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;
%%
Rind = 4;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);

a = 0.7; b = 1;
subind = [2:3,5:9,11,13:15];
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = [17:19,21:25,27,29:30];
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;
%%
Rind = 5;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);

a = 0.7; b = 1;
subind = [2:3,5:9,11,13:15];
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = [17:30];
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;

%%
Rind = 6;
EPr_t = zeros(N);
EPr_t(:, :) = EPr_t(:, :) + 0.02 * rand(N);
a = 0.7; b = 1;
subind = 1:13;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 14:22;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));
subind = 23:30;
EPr_t(subind, subind) = a + (b-a).*rand(numel(subind));

EPr_t =  triu(EPr_t, 1) + triu(EPr_t, 1)';
EPr_t(diagdex) = 0;
% figure;imagesc(1-EPr_t);colormap('gray');
EPr(:,:, Rind) = EPr_t;

figure(1);
for t = 1:R
    subplot(1,R,t);
    imagesc(1 - EPr(:,:,t));colormap('gray')
end

%%
T = 50;
regime_seq_indices = zeros(1,T);
regime_seq_indices(1:4)     = 1;
regime_seq_indices(5:8)   = 6;
regime_seq_indices(9:11)   = 3;
regime_seq_indices(12:14) = 1;
regime_seq_indices(15:19) = 2;
regime_seq_indices(20:25) = 4;
regime_seq_indices(26:30) = 2;
regime_seq_indices(31:33) = 1;
regime_seq_indices(34:38) = 3;
regime_seq_indices(39:42) = 1;
regime_seq_indices(47:50) = 4;
regime_seq_indices(43:46) = 5;
% figure(2);plot(1:T, regime_seq_indices, 'bo');
% ylim([0 R+1])
% set(gca, 'YGrid', 'on', 'XGrid', 'off');

%%
Adj = zeros(N, N, T);
for t = 1:T
    Adj_t = rand(N) < EPr(:,:, regime_seq_indices(t));
    Adj_t = triu(Adj_t, 1) + triu(Adj_t, 1)';
    Adj(:,:,t) = Adj_t;
end

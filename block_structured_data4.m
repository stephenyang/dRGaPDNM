
N = 50;
T =  30;
Adj = zeros(N, N, T);

%%
C = 5;
bd = [0 12 15 21 30 50];
for t = 1:10
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc == 2 || cc==4
    Adj(dex, dex, t) = rand(numel(dex)) < 0.02;
else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%%
C = 4;
bd = [0 10 20 32 50];
for t = 11:20
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
    Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%%
C = 5;
bd = [0 12 15 21 30 50];
for t = 21:30
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc == 2 || cc==4
    Adj(dex, dex, t) = rand(numel(dex)) < 0.02;
else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end
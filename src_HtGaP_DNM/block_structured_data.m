figure(1);
N = 60;
T = 6;
Adj = zeros(N, N, T);
C = 3;
bd = [0 15 35 80];
t = 1;
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
    Adj(dex, dex, t) = 1;
   Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
end
subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
for t = 2:T
    bd = bd + [0 6 7 0];
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        Adj(dex, dex, t) = 1;
        Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
    end
    subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end
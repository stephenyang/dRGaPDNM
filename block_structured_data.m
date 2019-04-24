figure(1);
N = 60;
T = 5;
B = zeros(N, N, T);
C = 3;
bd = [0 12 32 80];
t = 1;
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
    B(dex, dex, t) = 1;
    B(:, :, t) = B(:, :, t) - diag(diag(B(:, :, t)));
end
subplot(T, 4, (t-1)* 4 + 1);imagesc(B(:,:,t));
for t = 2:T
    bd = bd + [0 5 5 0];
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        B(dex, dex, t) = 1;
        B(:, :, t) = B(:, :, t) - diag(diag(B(:, :, t)));
    end
    subplot(T, 4, (t-1)* 4 + 1);imagesc(B(:,:,t));
end


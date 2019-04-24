% figure(1);
N = 80;
T = 35;
Adj = zeros(N, N, T);


%%
C = 4;
bd = [0 12 43 54 80];
for t = 1:7
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.035);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));

end

%%
C = 5;
bd = [0 8 16 24 48 80];
for t = 8:15
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.035);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));

end

%%
C = 4;
bd = [0 12 43 54 80];
for t = 16:20
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
    else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.035);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%%
bd = [0 8 16 36 60 80];C = 5;
for t = 21:25
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
        
    end
    Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.035);
    Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
  %  subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%%
C = 4;
bd = [0 12 43 54 80];
for t = 25:30
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
    else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.035);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end


%%
bd = [0 8 16 36 60 80];C = 5;
for t = 31:35
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        Adj(dex, dex, t) = rand(numel(dex)) < 0.65;
        % Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
    end
    Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
    Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%    subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end


Adj = double(Adj > 0);
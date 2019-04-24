% figure(1);
N = 80;
T = 23;
Adj = zeros(N, N, T);
C = 4;
bd = [0 12 43 54 80];

for t = 1:4
    
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = 1;
else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));

end

bd = [0 12 43 54 80];
for t = 5:8
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = 1;
    else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

bd = [0 8 16 36 60 80];C = 5;
for t = 9:12
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        Adj(dex, dex, t) = 1;
        
    end
    Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
    Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
  %  subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

C = 4;
bd = [0 12 43 54 80];
for t = 13:16
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = 1;
    else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end


bd = [0 8 16 36 60 80];C = 5;
for t = 17:19
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        Adj(dex, dex, t) = 1;
        % Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
    end
    Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
    Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%    subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

C = 4;
bd = [0 12 43 54 80];
for t = 20:23
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
if cc ~= 3
    Adj(dex, dex, t) = 1;
    else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.025;
end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.025);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
%subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end
Adj = double(Adj > 0);

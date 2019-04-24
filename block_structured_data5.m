
N = 50;
T =  50;
Adj = zeros(N, N, T);

%% 1-10
C = 5;
bd = [0 12 20 32 40 50];
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

%% 11-20
% C = 5;bd = [0 15 21 31 38 50];
C = 6;bd = [0 15 24 34 41 48 50];
for t = 11:20
    for cc = 1:C
        dex = bd(cc)+1:bd(cc+1);
        if cc == 2 || cc==4 || cc==6
            Adj(dex, dex, t) = rand(numel(dex)) < 0.02;
        else
            Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
        end
    end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%% 21- 30
C = 5;
bd = [0 10 20 28 36 50];
for t = 21:30
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);

    Adj(dex, dex, t) = rand(numel(dex)) < 0.55;

end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
end

%% 31-40
C = 5;
bd = [0 12 20 32 40 50];
for t = 31:40
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

%% 41- 50
% C = 5;
% bd = [0 10 20 28 36 50];
% for t = 41:50
% for cc = 1:C
%     dex = bd(cc)+1:bd(cc+1);
% if cc == 2 % || cc==4
%     Adj(dex, dex, t) = rand(numel(dex)) < 0.02;
% else
%     Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
% end
% end
% Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
% Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% % subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
% end
C = 3;
bd = [0 16 34 50];
for t = 41:50
for cc = 1:C
    dex = bd(cc)+1:bd(cc+1);
% if cc == 2 % || cc==4
%     Adj(dex, dex, t) = rand(numel(dex)) < 0.02;
% else
    Adj(dex, dex, t) = rand(numel(dex)) < 0.55;
% end
end
Adj(:, :, t) = Adj(:, :, t)+double(rand(N) < 0.03);
Adj(:, :, t) = Adj(:, :, t) - diag(diag(Adj(:, :, t)));
% subplot(T, 4, (t-1)* 4 + 1);imagesc(Adj(:,:,t));
remove_id = [3 5 7 8 9 11 13 15 18 20 21 23 24];
for id = 1:numel(remove_id)
    Adj(remove_id(id), :, t) = 0;
    Adj(:, remove_id(id), t) = 0;
end
end


% 
Adj = double(Adj > 0);
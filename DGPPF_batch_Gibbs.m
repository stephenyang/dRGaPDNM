function [result] = DGPPF_batch_Gibbs(B, options)

%%
% dynamic gamma process Poisson factorization model
% inference algo: batch Gibbs sampling
%%

burnin = options.burnin;
mcmcsamps = options.mcmcsamps;
maxiter = burnin + mcmcsamps;

N = size(B{1}, 1);
T = numel(B);
K = min(N, 30);

Ktrace = zeros(1, maxiter);

if isfield(options,'TrainRatio')
    idx_train = options.idx_train;
    idx_test = options.idx_test;
    BTrain_Mask_triu = cell(T, 1);
    BTrain_Mask = cell(T, 1);
    BTrain = cell(T, 1);
    for t = 1:T
        BTrain_Mask_triu{t} = (zeros(N));
        BTrain_Mask_triu{t}(idx_train{t}) = 1;
        BTrain_Mask{t} = full((BTrain_Mask_triu{t} + BTrain_Mask_triu{t}'));
        BTrain{t} = triu(B{t}, 1);
        BTrain{t}(idx_test{t}) = 0;
    end
end

r_KT = (ones(K,1)./K) * ones(1, T);
gamma_K = ones(K, 1);
a_0 = 1;
b_0 = 1;
c_0 = 1;
d_0 = 1;
e_0 = 1;
f_0 = 1;
g_0 = 1;
c = 1;
c_n = ones(N, 1);
r_0 = ones(K, 1)./K;L_0 = zeros(K, 1);
phi = randg(ones(N, K));
Probsamps = zeros(N, N, T);
ProbAve = cell(T,1);
Btriu =cell(1,T);
obs_links = cell(1,T);
dex_lef = cell(1,T);
dex_rig = cell(1,T);
for t = 1:T
    Btriu{t} = triu(BTrain{t}, 1);
    [dex_lef{t}, dex_rig{t}, obs_links{t}] = find(Btriu{t});
    idx{t} = sub2ind([N, N], dex_lef{t}, dex_rig{t});
end

m_n_dot_k = zeros(K, N, T);
m_dot_dot_k = zeros(K, T);
tic;
for iter = 1:maxiter
    
    
    m_n_dot_k = zeros(K, N, T);
    for t = 1:T
        if iter >options.burnin
            Btriu{t} = triu(BTrain{t}, 1);
            [dex_lef{t}, dex_rig{t}, obs_links{t}] = find(Btriu{t});
        end
        rate = sum(phi(dex_lef{t}, :) * diag(r_KT(:, t)) .* phi(dex_rig{t}, :), 2);
        hidden_links{t} = truncated_poisson_rnd(rate);
        
%        [m_n_dot_k(:, :, t)] = Multrnd_mijk(sparse(dex_lef{t}, dex_rig{t}, hidden_links{t}, N, N), phi, r_KT(:, t));
 
            for ij = 1:numel(idx{t})
                pmf = phi(dex_lef{t}(ij), :)  .* r_KT(:, t)'.* phi(dex_rig{t}(ij), :);
                pmf = pmf(:)./sum(pmf);
                mij_kk = multrnd_histc(obs_links{t}(ij),pmf(:))';
                % mij_kk = mnrnd(hidden_links{t}(ij), pmf);
                m_n_dot_k(:,dex_lef{t}(ij), t) = m_n_dot_k(:,dex_lef{t}(ij), t) + mij_kk';
                m_n_dot_k(:,dex_rig{t}(ij), t) = m_n_dot_k(:,dex_rig{t}(ij), t) + mij_kk';
            end
        % m_dot_k_k_dot = zeros(K,K);
      
        m_dot_dot_k(:, t) = sum(m_n_dot_k(:, :, t), 2)./2;
    end
    
    Ktrace(iter) = nnz(sum(m_dot_dot_k, 2));
    
   
        s_k = sum(phi .* bsxfun(@minus, sum(phi, 1), phi), 1)'./2;
        
        P_T = zeros(K, T);
        P_T(:, T) = s_k./(c+s_k);
        for t = (T-1):-1:1
            P_T(:, t) = 1 - c./(c + s_k - log(max(1 - P_T(:, t+1), realmin)));
        end
        P_0 = 1 + c./(log(max(1 - P_T(:, 1),realmin)) - c);
        
        L_T = zeros(K, T);
        for k = 1:K
            L_T(k,T) = CRT_sum_mex_matrix(sparse(m_dot_dot_k(k, T)), r_KT(k, T));
        end
        for t = (T-1):-1:2
            for k = 1:K
                L_T(k, t) = CRT_sum_mex_matrix(sparse(m_dot_dot_k(k, t) + L_T(k, t+1)), r_KT(k, t-1));
            end
        end
        for k = 1:K
            L_T(k, 1) = CRT_sum_mex_matrix(sparse(m_dot_dot_k(k, 1) + L_T(k, 2)), r_0(k));
            L_0(k) = CRT_sum_mex_matrix(sparse(L_T(k, 1)), gamma_K(k));
        end
        %% forward
        r_0 = randg(gamma_K + L_T(:, 1))./(c - log(max(1 - P_T(:, 1), realmin)));
        r_KT(:, 1) = randg(r_0 + m_dot_dot_k(:, 1) + L_T(:, 2))./(c + s_k - log(max(1- P_T(:, 2), realmin)));
        for t = 2:(T-1)
            r_KT(:, t) = randg(r_KT(:, t-1) + m_dot_dot_k(:, t) + L_T(:, t+1))./(c + s_k - log(max(1- P_T(:, t+1), realmin)));
        end
        r_KT(:, T) = randg(r_KT(:, T-1) + m_dot_dot_k(:, T))./(c + s_k);
     
    gamma_K = randg(a_0 + L_0)./(b_0 - log(max(1 - P_0, realmin)));
    
    
    m_dot_k_dot_n = sum(m_n_dot_k, 3);
    sum_r_KT = sum(r_KT, 2)';
    temp = sum(phi, 1) .* sum_r_KT;
    for i = randperm(N)
        temp = temp - phi(i, :) .* sum_r_KT;
        phi(i, :) = randg(e_0 + m_dot_k_dot_n(:, i))'./(c_n(i) + temp);
        temp = temp + phi(i, :) .* sum_r_KT;
    end
    
    
    c_n = randg(f_0 + e_0 * K)./(g_0 + sum(phi, 2));
    
     
    c = randg(c_0 + sum(gamma_K) + sum(r_0) + sum(sum(r_KT(:, 1:(T-1)))))./(d_0 + sum(r_0) + sum(sum(r_KT(:, 1:T))));
    
     
    for t = 1:T
        Prob = phi * sparse(diag(r_KT(:, t))) * phi' + eps;
        Prob = 1 - exp(- Prob);
        if iter > burnin
            Probsamps(:, :, t) = Probsamps(:, :, t) + Prob;
            ProbAve{t} = Probsamps(:, :, t)/(iter - burnin);
        else
            ProbAve{t} = Prob;
        end
        
        if iter > options.burnin && options.TrainRatio ~=1
            BTrain{t}(idx_test{t}) = rand() < ProbAve{t}(idx_test{t});
        end
    end

    if mod(iter, 10) == 0 

        fprintf('DGPPF_iter: %03d.\n', [iter]);
        if options.display==1
        figure(123);
        [~, rdex] = sort(m_dot_dot_k, 'descend');
        [~, z] = max(m_dot_k_dot_n(rdex, :),[],1);
        [~, Rankdex] = sort(z);
        frame = min(T,8);
            if T> 50
                t_ind = [1:10:70 70];
            else
                 t_ind = [1:T];
            end
        for t = 1:frame
%                         Prob = phi * sparse(diag(r_KT(:, t_ind(t)))) * phi' + eps;
%                         Prob = 1 - exp(- Prob);
%                         if iter > burnin
%                             Probsamps(:, :, t_ind(t)) = Probsamps(:, :, t_ind(t)) + Prob;
%                             ProbAve = Probsamps(:, :, t_ind(t))/(iter - burnin);
%                         else
%                             ProbAve = Prob;
%                         end
            subplot(frame, 4, (t-1)* 4 + 1);
            imagesc(BTrain{t_ind(t)} + BTrain{t_ind(t)}');
            subplot(frame, 4, (t-1)* 4 + 2);
            imagesc(ProbAve{t_ind(t)});
            % imagesc(B(Rankdex,Rankdex,t));
            subplot(frame, 4, (t-1)* 4 + 3);
            imagesc(log(phi(:, rdex(:,t_ind(t))) * diag(r_KT(rdex(:,t_ind(t)), t_ind(t))) + 1e-2));
            subplot(frame, 4, (t-1)* 4 + 4);
            bar(r_KT(rdex(:, t_ind(t)), t_ind(t))./sum(r_KT(rdex(:, t_ind(t)), t_ind(t))));
            
        end
        drawnow;
        end
    end
end
result.timecost =  toc
result.r_KT = r_KT;
% result.rdex = rdex;
result.phi = phi;
result.ProbAve = ProbAve;

%% quantify results
if isfield(options,'TrainRatio')
    if options.TrainRatio == 1
        idx_test = idx_train;        
    end
        rate_test = cell(T,1);
        links = cell(T, 1);
        coll_rate = [];
        coll_links = [];
        for t = 1:T
            B{t} = full(B{t});
            rate_test{t} = ProbAve{t}(idx_test{t});
            links{t} = B{t}(idx_test{t});
            coll_rate = [coll_rate;rate_test{t}];
            coll_links = [coll_links;links{t}];
        end
        % coll_rate = double(coll_rate > 1e-3);
        figure(234);
        subplot(1,2,1);
        [XX, YY, TT, AUCroc] = perfcurve(coll_links, coll_rate, 1);
        plot(XX,YY);
        axis([0 1 0 1]), grid on, xlabel('FPR'), ylabel('TPR'), hold on;
        x = [0:0.1:1];plot(x,x,'b--'), hold off; title(['AUCroc = ', num2str(AUCroc)])
        
        subplot(1,2,2);
        [prec, tpr, fpr, thresh] = prec_rec(coll_rate, coll_links,  'numThresh',3000);
        plot([0; tpr], [1 ; prec]); % add pseudo point to complete curve
        xlabel('recall');
        ylabel('precision');
        title('precision-recall graph');
        AUCpr = trapz([0;tpr],[1;prec]);
        F1= max(2*tpr.*prec./(tpr+prec));
        title(['AUCpr = ', num2str(AUCpr), '; F1 = ', num2str(F1)])
        fprintf('aucroc: %f, aucprec: %f.\n', [AUCroc, AUCpr]);
        result.aucroc = AUCroc;
        result.aucprec = AUCpr;
end

end
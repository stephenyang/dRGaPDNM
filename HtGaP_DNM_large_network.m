function [result] = HtGaP_DNM_large_network(Adj, options)

%%
% Hierarchial thinned gamma process dynamic relational model
% inference algo: batch Gibbs sampling
%%

%% simulate data
re_order = 0;

%% HtGP-PFA
burnin = options.burnin;
mcmcsamps = options.mcmcsamps;
maxiter = burnin + mcmcsamps;

T = numel(Adj);
N = size(Adj{1}, 1);
K = min(N, 50);
K_active = zeros(1, maxiter);

diagdex = sparse(1:K, 1:K, true);
triu1dex = triu(true(K), 1);

if isfield(options,'TrainRatio')
    idx_train = options.idx_train;
    idx_test = options.idx_test;
    ATrain_Mask_triu = cell(T, 1);
    ATrain_Mask = cell(T, 1);
    ATrain = cell(T, 1);
    for t = 1:T
        ATrain_Mask_triu{t} = (zeros(N));
        ATrain_Mask_triu{t}(idx_train{t}) = 1;
        ATrain_Mask{t} = full((ATrain_Mask_triu{t} + ATrain_Mask_triu{t}'));
        ATrain{t} = triu(Adj{t}, 1);
        ATrain{t}(idx_test{t}) = 0;
    end
end


    r_K = ones(K, 1)./K;
    gamma_0 = 1;  % gamma_0/K shape param of r_k
    xi = 1;
    c_n = ones(1,N);
    alpha_0 = 1;     % rate param of r_k
    
     
    jitter = 0.1;
    psi_dict = [1:1:20, 0.1,0.5,0.01,0.05,0.001,0.005];
    len_psi_dict = numel(psi_dict);
    U = 1./ones(1, len_psi_dict); % uniform prior for kernel width
    psi_inds = ones(K,1);
    
    Tvec = ones(T, 1);
    mus = (1:T)';
    kernel_dist = pdist2((1:T)', mus).^2;
    Km = compute_kernel_mat((1:T)', struct('K', K, 'mus', mus, 'psi_inds', psi_inds, 'psi_dict', psi_dict,'kernel_dist',kernel_dist));
    L = T;
    W = zeros(L+1, K); % kernel weights

    c0 = 1e-4;
    d0 = 1e-4;
    invSigma = (c0/d0)*ones(L+1,K); % variance of omega (kernel width)

    b_KT = zeros(K, T);
    for k = 1:K
        nu_k = betarnd(1,1);
        inds = [true; rand(L, 1)< nu_k];
        W(inds, k) = randn(sum(inds, 1), 1) .* (1./(invSigma(inds, k)));
        W(~inds, k) = 0;
        tmp = Km{k} * W(:, k);
        tmin = min(tmp);
        if tmin <0
            tmp = tmp + abs(tmin);
        end
        tmp = tmp./max(tmp);
        g_k = tmp;  % thining probs
        g_k(g_k==0) = 1e-16;
        g_k(g_k==1) = 1 - 1e-16;
        b_KT(k, :) = rand(1, T) < g_k';
    end
    b_KT = ones(K, T);
    zcode = [0 0 1];
    lp = zeros(1, 3);
    
    for t = 1:T
        if sum(b_KT(:, t)) == 0
            b_KT(randi([1 K]), t) = 1;  % at least one group per time
        end
    end

    b_tilde_KT = b_KT - 0.5; % auxiliary vars
    
    b_K_dot_r_K = (r_K * ones(1, T)) .* b_KT;
    
    beta = 1; % rate of lambda_kk
    betavec = beta .* ones(1, K);
    
    lambda_KKT = zeros(K, K, T);
    for t = 1:T
        lambda_KKT(:, :, t) = b_K_dot_r_K(:, t) * b_K_dot_r_K(:, t)';
        lambda_KKT(:, :, t) = triu(lambda_KKT(:, :, t), 1) + triu(lambda_KKT(:, :, t))';
        lambda_KKT(sparse(1:K, 1:K, true)) = xi * b_K_dot_r_K(:, t);
    end

    phi_0 = randg(ones(N, K));
    phi = zeros(N, K, T);
    for t = 1:T
        phi(:, :, t) = gamrnd(1,1, [N, K]);
    end
    
    % m_hat_nkt = zeros(N, K, T);
    P_n_k_t = zeros(N,K,T);
    P_n_k_0 = zeros(N, K);
    g_0= 1;
    f_0 = 1;
    u_0 = 1;
    v_0 = 1;
    a_0 = 1;
    b_0 = 1;
    % varbeta = 1;
    Y_n_k_t = zeros(N, K, T);
    
    
    Probsamps = zeros(N, N, T);
    Prob = zeros(N, N, T);
    ProbAve = cell(1,T);
    

    for t = 1:T
        Atriu{t} = triu(ATrain{t}, 1);
        [dex_left{t}, dex_rig{t}, ~] = find(Atriu{t});
    end
    M = cell(1,T);

    m_n_k_dot_dot = zeros(K, N, T);
    m_dot_k_k_dot = zeros(K, K, T);

for iter = 1:maxiter
    

    for t = 1:T
%         if iter >options.burnin
%             Atriu{t} = triu(ATrain{t}, 1);
%             [dex_left{t}, dex_rig{t}, ~] = find(Atriu{t});
%         end
        rate = sum(phi(dex_left{t}, :, t) * lambda_KKT(:, :, t) .* phi(dex_rig{t}, :, t), 2);
        M{t} = truncated_poisson_rnd(rate);
        [m_n_k_dot_dot(:, :, t), diag_kk] = Multrnd_mik1k2j(sparse(dex_left{t}, dex_rig{t}, M{t}, N, N), phi(:,:,t), lambda_KKT(:, :, t));
        diag_kk(sparse(1:K, 1:K, true)) = diag_kk(sparse(1:K, 1:K, true))./2;
        m_dot_k_k_dot(:, :, t) = diag_kk;
    end
    K_active(iter) = nnz(sum(sum(m_n_k_dot_dot, 3), 2));

 
    theta_KK = zeros(K);
    for t = 1:T
        theta_KK_t = phi(:, :, t)' * bsxfun(@minus, sum(phi(:, :, t), 1), phi(:, :, t));
        theta_KK_t = theta_KK_t .* (b_KT(:, t) * b_KT(:, t)');
        theta_KK = theta_KK + theta_KK_t;
    end
    theta_KK(diagdex) = theta_KK(diagdex)./2;
    
    m_dot_k_k_dot_dot = sum(m_dot_k_k_dot, 3);
    L_KK = zeros(K);
    p_kk_tilde_one_minus = zeros(K);
    temp = zeros(K, 1);
    betavec = beta .* ones(1, K);
    for k = randperm(K)
        R_K = r_K';
        R_K(k) = xi;
        p_kk_tilde_one_minus(k, :) = betavec./(betavec + theta_KK(k, :));
        L_KK(k, :) = CRT_sum_mex_matrix(sparse(m_dot_k_k_dot_dot(k,:)), r_K(k) * R_K);
        temp(k) = sum(R_K .* log(max(p_kk_tilde_one_minus(k, :), realmin)));
        r_K(k) = randg(gamma_0/K + sum(L_KK(k, :)))./(alpha_0 - temp(k));
    end

    ell = sum(CRT_sum_mex_matrix(sparse(m_dot_k_k_dot_dot(diagdex))', xi * r_K'));
    xi = randg(1e-2 + ell)./(1e-2 - sum(r_K .* log(max(diag(p_kk_tilde_one_minus), realmin))));

  
    R_KK = r_K * r_K';
    R_KK(diagdex) = xi * r_K;
    lambda_KK = zeros(K);
    lambda_KK(diagdex) = randg(m_dot_k_k_dot_dot(diagdex) + R_KK(diagdex))./(beta + theta_KK(diagdex));
    lambda_KK(triu1dex) = randg(m_dot_k_k_dot_dot(triu1dex) + R_KK(triu1dex))./(beta + theta_KK(triu1dex));
    lambda_KK = lambda_KK + triu(lambda_KK, 1)';
    for t = 1:T
        lambda_KKT(:, :,t) = (b_KT(:, t) * b_KT(:, t)') .* lambda_KK;
    end
    

    beta = randg(1 + sum(R_KK(diagdex)) + sum(R_KK(triu1dex)))./(1+ sum(lambda_KK(diagdex)) + sum(lambda_KK(triu1dex)) );
    

    L_k_tilde = CRT_sum_mex_matrix(sparse(sum(L_KK, 2)'), gamma_0/K);
    p_ttilde_k_one_minus = alpha_0./(alpha_0 - temp);
    gamma_0 = randg(1e-3 + sum(L_k_tilde))./(1e-3 - 1/K*sum(log(max(p_ttilde_k_one_minus, realmin))));
       

    alpha_0 = randg(1e-2 + gamma_0)./(1e-2 + sum(r_K));
    

    for rep = 1:5  
    lp = zeros(1, 3);
    for t = 1:T
        theta_KK_t = phi(:, :, t)' * bsxfun(@minus, sum(phi(:, :, t), 1), phi(:, :, t));
        %%
        m_n_k_dot_dot_t = m_n_k_dot_dot(:, :, t);
        pos_ind = sum(m_n_k_dot_dot_t, 2) > 0;
        %%
        for k = 1:K
            phi_k = gaussiancdf(Km{k} * W(:, k));
            phi_k(phi_k==0) = 1e-16;
            phi_k(phi_k==1)  = 1 - 1e-16;
            
            lPhix = log(phi_k(t));
            lPhix_c = log(1 - phi_k(t));
            
            if pos_ind(k)==1
                b_KT(k, t) = 1;
            else
                rate = r_K(k) * xi * theta_KK_t(k, k);
                %%% change to NB pdf
                denom = log(exp(-rate) * phi_k(t) + (1-phi_k(t)) + eps);
                %% 3
                lp(1) = -rate + lPhix_c - denom;
                %% 2
                lp(2) = lPhix_c + log(1 - exp(-rate) + eps) - denom;
                %% 1
                lp(3) = -rate + lPhix - denom;
                
                p = exp(lp);
                r = sample_discrete(p,1);
                b_KT(k, t) = zcode(r);
            end
        end
    end
    end


    for k = randperm(K)
        kinds = Km{k};
        inv_sigma = diag(invSigma(:, k));
        Sigma = (inv_sigma + kinds' * kinds + jitter * eye(L + 1));
        Mu = Sigma\(kinds'*b_tilde_KT(k, :)');
        W(:, k) = mvnrnd(Mu', inv(Sigma))';
        a = -inf(T, 1);
        b = zeros(T, 1);
        idx = b_KT(k, :)==1;
        a(idx)=0; b(idx)=1e4;
        b_tilde_KT(k, :) = randTN(kinds * W(:, k), 1, a, b);
    end
    

    for k = 1:K
        b_KT_k = b_KT(k,:)';
        lp = zeros(1, len_psi_dict);
        for p = 1:len_psi_dict
            Km_temp = [Tvec exp(-kernel_dist/psi_dict(p))];
            Phi_k = gaussiancdf(Km_temp * W(:, k));
            lp(p) = sum(b_KT_k .* log(Phi_k + eps) + (1-b_KT_k).*log(1 - Phi_k + eps));
            lp(p) = lp(p) + log(U(p));
        end
        
        lp = lp - max(lp);
        pr = exp(lp);
        pr = pr./sum(pr);
        psi_inds(k) = sum(rand > cumsum(pr)) + 1;
    end
    Km = compute_kernel_mat(Tvec, struct('K', K, 'mus', mus, 'psi_inds', psi_inds, 'psi_dict', psi_dict, 'kernel_dist', kernel_dist));
    

    c = c0 + 0.5;
    d = d0 + 0.5*W.^2;
    invSigma = gamrnd(c, 1./d, [L+1 K]);


    %% P_n_k_t :: backward
            phi_times_lambda_KK_T = phi(:, :, T) * lambda_KKT(:, :, T);
            temp = sum(phi_times_lambda_KK_T, 1);
            for i = randperm(N)
                temp_i = temp - phi_times_lambda_KK_T(i, :);
                P_n_k_t(i, :, T) = temp_i ./ (c_n(i) + temp_i);
            end
            for t = (T-1):-1:1
                phi_times_lambda_KK_t = phi(:, :, t) * lambda_KKT(:, :, t);
                temp = sum(phi_times_lambda_KK_t, 1);
                for i = 1:N
                    temp_i = temp - phi_times_lambda_KK_t(i,:) - log(max(1-P_n_k_t(i, :, t+1), realmin));
                    P_n_k_t(i, :, t) = temp_i./(c_n(i) + temp_i);
                end
            end
            
            %% Y_n_k_t :: backward
            Y_n_k_t(:,:,T) = CRT_matrix(m_n_k_dot_dot(:, :, T)', phi(:, :, T-1));
            for t = (T-1):-1:2
                Y_n_k_t(:,:,t) = CRT_matrix(m_n_k_dot_dot(:, :, t)' + Y_n_k_t(:, :, t+1), phi(:, :, t-1));
            end
            Y_n_k_t(:, :, 1) = CRT_matrix(m_n_k_dot_dot(:, :, 1)' + Y_n_k_t(:, :, 2),phi_0);
            Y_n_k_0 = CRT_matrix(Y_n_k_t(:, :, 1), g_0(ones(N, K)));
                        
            %% phi :: forward
            phi_0 = randg(g_0 + Y_n_k_t(:, :, 1))./(f_0 - log(max(1 - P_n_k_t(:, :, 1), realmin)));
            
            %%
            phi_lambda = phi(:, :, 1) * lambda_KKT(:, :, 1);
            temp = sum(phi_lambda, 1);
            for i = randperm(N)
                temp = temp - phi_lambda(i, :);
                phi(i, :, 1) = randg(phi_0(i, :) + Y_n_k_t(i, :, 2) + m_n_k_dot_dot(:, i, 1)')./(c_n(i) + temp - log(max(1 - P_n_k_t(i, :, 2), realmin)));
                temp = temp + phi(i, :, 1) * lambda_KKT(:, :, 1);
            end
            
            %%            
            for t = 2:(T-1)
                    phi_lambda = phi(:, :, t) * lambda_KKT(:, :, t);
                    temp = sum(phi_lambda, 1);

                    for i = randperm(N)
                        temp = temp - phi_lambda(i, :);
                        phi(i, :, t) = randg(phi(i, :, t-1) + Y_n_k_t(i, :, t+1) + m_n_k_dot_dot(:, i, t)')./(c_n(i) + temp - log(max(1 - P_n_k_t(i, :, t+1), realmin)));
                        temp = temp + phi(i, :, t) * lambda_KKT(:, :, t);
                    end  
            end
            
            %%
            phi_lambda = phi(:, :, T) * lambda_KKT(:, :, T);
            temp = sum(phi_lambda, 1);
            for i = randperm(N)
                  temp = temp - phi_lambda(i, :);
                  phi(i, :, T) = randg(phi(i, :, T-1) + m_n_k_dot_dot(:, i, T)')./(c_n(i) + temp);
                  temp = temp + phi(i, :, T) * lambda_KKT(:, :, T);
            end
            
            %% 
            g_0 = randg(u_0 + sum(sum(Y_n_k_0)))./(v_0 - sum(sum(log(max(1 - P_n_k_0, realmin)))));
            
            %%
            f_0 = randg(g_0 * N * K + a_0)./(b_0 + sum(sum(phi_0)));
            

    
    m_n_k_dot_dot_dot = sum(m_n_k_dot_dot, 3);
    [~, rdex] = sort(sum(m_n_k_dot_dot_dot, 2), 'descend');
    [~, z] = max(m_n_k_dot_dot_dot(rdex, :), [], 1);
    [~, Rankdex] = sort(z);
    
    for t = 1:T
        Prob(:,:,t) = phi(:,:,t) * lambda_KKT(:,:,t) * phi(:,:,t)' + eps;
        Prob(:,:,t) = 1 - exp(-Prob(:,:,t));
        if iter > burnin
            Probsamps(:,:,t) = Probsamps(:,:,t) + Prob(:,:,t);
            ProbAve{t} = Probsamps(:,:,t)/(iter - burnin);
        else
            ProbAve{t} = Prob(:,:,t);
        end
        
%         if iter > options.burnin && options.TrainRatio ~=1
%             ATrain{t}(idx_test{t}) = rand() < ProbAve{t}(idx_test{t});
%         end
    end
        
    if mod(iter, 100) == 0
        % toc
        fprintf('iter: %04d.\n', [iter]);
        if options.display
            slices = min(8, T);
            for t = 1:slices
                figure(123);
                if re_order== 1
                    subplot(slices, 4, (t-1)* 4 + 1); imagesc(ATrain{t}(Rankdex, Rankdex));
                    subplot(slices, 4, (t-1)* 4 + 2); imagesc(ProbAve{t}(Rankdex, Rankdex));
                    subplot(slices, 4, (t-1)* 4 + 3); imagesc(log(phi(Rankdex, rdex, t) * lambda_KKT(rdex, rdex, t) + 1e-2));
                    subplot(slices, 4, (t-1)* 4 + 4); imagesc(log(lambda_KKT(rdex, rdex,t) + 1e-2));
                else
                    subplot(slices, 4, (t-1)* 4 + 1); imagesc(ATrain{t});
                    subplot(slices, 4, (t-1)* 4 + 2); imagesc(ProbAve{t});
                    subplot(slices, 4, (t-1)* 4 + 3); imagesc(log(phi(:,:, t) * lambda_KKT(:,:, t) + 1e-2));
                    subplot(slices, 4, (t-1)* 4 + 4); imagesc(log(lambda_KKT(:,:,t) + 1e-2));
                end
                
                figure(234);colormap('gray')
                % subplot(1,5, mod(iter/50, 5)+1);
                imagesc(b_KT);title(['iter: ' num2str(iter)]);
                drawnow;
            end
        end
    end
end

%% collect mcmc samps & results
result.ProbAve = ProbAve;
result.r_K = r_K;
result.K = sum(r_K./sum(r_K) > 0.05);
result.phi = phi;
result.lambda_KK = lambda_KKT(:, :, end);
result.lambda_KKT = lambda_KKT;
result.Rankdex = Rankdex;
result.rdex = rdex;
result.W = W;
result.Km = Km;
result.L = L;
result.T = T;
result.psi_inds = psi_inds;
%%

end
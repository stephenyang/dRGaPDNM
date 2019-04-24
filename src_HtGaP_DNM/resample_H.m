%% pi_KT
pi_KT = zeros(K, T);
zcode = [0 0 1];

for t = 1:T
    m_dot_k_k_dot_t = m_dot_k_k_dot(:, :, t);
    pos_ind = m_dot_k_k_dot_t(diagdex) > 0;
    for k = 1:K

        phi_k = normcdf(Km{k} * W(:, k));
        phi_k(phi_k==0) = 1e-16;
        phi_k(phi_k==1)  = 1 - 1e-16;
        
        lPhix = log(phi_k(t));
        lPhix_c = log(1 - phi_k(t));
    
        if pos_ind(k)==1
            pi_KT(k, t) = 1;
        else
            rate = r_K(k) * xi *  phi(:, :, t)' * bsxfun(@minus, sum(phi(:, :, t), 1), phi(:, :, t));
            denom = log(exp(-rate) * phi_k(t) + (1-phi_k(t)) + eps);
            %% 3
            lp(1) = -rate + lPhix_c - denom;
            %% 2
            lp(2) = lPhix_c + log(1 - exp(-rate) + eps) - denom;
            %% 1
            lp(3) = -rate + lPhix - denom;
            
            p = exp(lp);
            cump(1) = p(1);
            cump(2) = sum(p(1:2));
            cump(3) =sum(p);
            
            r = find(rand * (1-eps) < cump, 1);
            pi_KT(k, t) = zcode(r);
        end
    end
end

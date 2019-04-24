function [Km] = compute_kernel_mat(X, params)

K = params.K;
% mus = params.mus;
kernel_dist = params.kernel_dist;
psi_inds = params.psi_inds;
psi_dict = params.psi_dict;

N = size(X, 1);
Km = cell(1, K);
Nones = ones(N, 1);
for k = 1:K
    temp = exp(-kernel_dist/psi_dict(psi_inds(k)));
    Km{k} = [Nones temp];
end

end
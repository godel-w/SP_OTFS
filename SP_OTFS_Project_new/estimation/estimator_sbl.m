function [supp, h_hat, state, diag_info] = estimator_sbl(y, Omega, cfg, noise_var, warm_state)
    [MN, Q] = size(Omega);
    % 列归一化
    col_norms = sqrt(sum(abs(Omega).^2, 1)).';
    col_norms(col_norms < 1e-9) = 1;
    Omega_n = Omega ./ col_norms.';
    
    sigma2 = noise_var;
    if ~isempty(warm_state), gamma = warm_state.gamma;
    else, gamma = (abs(Omega_n'*y).^2)/MN; end
    gamma(gamma < 1e-12) = 1e-12;
    
    for i = 1:cfg.est.sbl.max_iter
        g_old = gamma;
        % E-Step: Robust Update
        inv_Gamma = diag(1./gamma);
        A = (Omega_n' * Omega_n)/sigma2 + inv_Gamma;
        
        try
            L = chol((A+A')/2, 'lower');
            mu = L' \ (L \ ((Omega_n'*y)/sigma2));
            Sigma_diag = sum(abs(inv(L)).^2, 2);
        catch
            mu = A \ ((Omega_n'*y)/sigma2);
            Sigma_diag = real(diag(inv(A)));
        end
        
        % M-Step
        gamma = abs(mu).^2 + Sigma_diag;
        gamma(gamma < 1e-12) = 1e-12;
        
        if norm(gamma-g_old)/norm(g_old) < cfg.est.sbl.tol, break; end
    end
    
    th = max(gamma) * cfg.est.sbl.prune;
    supp = find(gamma > th);
    
    if length(supp) > cfg.est.sbl.max_active
        [~, s_idx] = sort(gamma(supp), 'descend');
        supp = supp(s_idx(1:cfg.est.sbl.max_active));
    end
    
    h_hat = Omega(:, supp) \ y; % LS Refine
    
    state.gamma = gamma;
    state.sigma2 = sigma2;
    diag_info.iters = i;
end
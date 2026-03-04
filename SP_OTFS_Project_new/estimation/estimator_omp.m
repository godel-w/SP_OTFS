function [supp, h_hat, info] = estimator_omp(y, Omega, cfg, noise_var)
    [MN, Q] = size(Omega);
    % 列归一化
    col_norms = sqrt(sum(abs(Omega).^2, 1)).';
    col_norms(col_norms < 1e-9) = 1;
    Omega_n = Omega ./ col_norms.';
    
    residual = y;
    supp = [];
    % 能量门限
    tol_energy = MN * noise_var * cfg.est.omp.resid_factor;
    
    iter = 0;
    for k = 1:cfg.est.omp.max_sparsity
        iter = k;
        proj = Omega_n' * residual;
        if ~isempty(supp), proj(supp) = 0; end
        [~, idx] = max(abs(proj));
        
        supp = [supp; idx];
        Omega_S = Omega(:, supp);
        h_hat = Omega_S \ y;
        
        residual = y - Omega_S * h_hat;
        if norm(residual)^2 < tol_energy, break; end
    end
    info.iters = iter;
end
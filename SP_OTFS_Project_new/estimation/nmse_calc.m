function nmse = nmse_calc(h_est, supp_est, grid, chan_true)
    % 映射回全空间计算 NMSE
    h_full_true = zeros(grid.n_cand, 1);
    for p = 1:chan_true.n_paths
        l_t = round(chan_true.l(p));
        k_t = round(chan_true.k(p));
        idx = find(grid.l == l_t & grid.k == k_t);
        if ~isempty(idx), h_full_true(idx(1)) = h_full_true(idx(1)) + chan_true.h(p); end
    end
    
    h_full_est = zeros(grid.n_cand, 1);
    h_full_est(supp_est) = h_est;
    
    err = norm(h_full_est - h_full_true)^2;
    ref = norm(h_full_true)^2;
    if ref < 1e-12, nmse = 0; else, nmse = err/ref; end
end
function grid = grid_generate(cfg)
    l_vec = 0 : cfg.search.l_max;
    k_vec = -cfg.search.k_max : cfg.search.k_max;
    [L, K] = meshgrid(l_vec, k_vec);
    grid.l = L(:);
    grid.k = K(:);
    grid.n_cand = length(grid.l);
end
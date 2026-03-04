function Heff = build_heff_from_support(supp, h_hat, grid, bases, cfg)
    M = cfg.grid.M; N = cfg.grid.N; MN = M*N;
    l_vec = grid.l(supp); k_vec = grid.k(supp);
    
    H_t = sparse(MN, MN);
    dop = exp(1j * 2 * pi * (0:MN-1)' / MN);
    idx = (1:MN)';
    
    for i = 1:length(supp)
        r = mod(idx - 1 + l_vec(i), MN) + 1;
        Pi = sparse(r, idx, ones(MN,1), MN, MN);
        Delta = spdiags(dop.^k_vec(i), 0, MN, MN);
        H_t = H_t + h_hat(i) * Pi * Delta;
    end
    Heff = bases.Brx * H_t * bases.Btx;
end
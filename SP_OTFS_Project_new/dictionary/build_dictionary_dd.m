function Omega = build_dictionary_dd(x_p, grid, bases, cfg)
    % 构造 DD 域字典: Col_i = Brx * Pi^l * Delta^k * Btx * x_p
    MN = cfg.grid.M * cfg.grid.N;
    Q = grid.n_cand;
    Omega = zeros(MN, Q);
    
    % 1. Pilot -> Time
    s_p = bases.Btx * x_p;
    dop_base = exp(1j * 2 * pi * (0:MN-1)' / MN);
    
    for i = 1:Q
        l = grid.l(i); k = grid.k(i);
        % Delta^k * s
        v = s_p .* (dop_base .^ k);
        % Pi^l * v
        v = circshift(v, l);
        % Time -> DD
        Omega(:, i) = bases.Brx * v;
    end
end
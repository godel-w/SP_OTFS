function pass = sanity_check_channel_alignment()
    cfg = cfg_default();
    cfg.grid.M = 8; cfg.grid.N = 4; MN=32; cfg.search.l_max=2;
    
    chan.l=2; chan.k=1; chan.h=1; chan.n_paths=1;
    x_p = randn(MN, 1);
    
    [Btx, Brx] = build_Btx_Brx(cfg.grid.M, cfg.grid.N);
    [s_cp, ~, cp_len] = otfs_mod(reshape(x_p,8,4), cfg);
    
    % Phys
    r_cp = apply_channel_phys(s_cp, chan, 0, cfg);
    [y_dd_phys, ~] = otfs_demod(r_cp, cfg, cp_len);
    
    % Matrix: Brx * Pi * Delta * Btx * x_p
    dop = exp(1j*2*pi*(0:MN-1)'/MN); idx=(1:MN)';
    r_idx = mod(idx-1+chan.l, MN)+1;
    Pi = sparse(r_idx, idx, ones(MN,1), MN, MN);
    Delta = spdiags(dop.^chan.k, 0, MN, MN);
    y_dd_mat = Brx * Pi * Delta * Btx * x_p;
    
    err = norm(y_dd_phys - y_dd_mat);
    pass = err < 1e-8;
    if pass, fprintf('[PASS] Channel Alignment (Err=%.2e)\n', err);
    else, fprintf('[FAIL] Channel Alignment (Err=%.2e)\n', err); end
end
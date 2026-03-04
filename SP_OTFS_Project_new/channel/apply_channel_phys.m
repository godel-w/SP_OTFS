function r_cp = apply_channel_phys(s_cp, chan, sigma_w2, cfg)
    % 物理层信道: r[n] = sum h * s[n-l] * exp(j*2pi*k*(n-l)/MN)
    % 关键: 对齐相位参考点 n=0 为去CP后的第一个样本
    
    MN = cfg.grid.M * cfg.grid.N;
    len_tot = length(s_cp);
    cp_len = len_tot - MN;
    
    % t_vec: CP部分为负索引，Body从0开始
    t_vec = (0:len_tot-1).' - cp_len;
    
    r_clean = zeros(len_tot, 1);
    
    for p = 1:chan.n_paths
        l = chan.l(p);
        k = chan.k(p);
        h = chan.h(p);
        
        l_int = round(l); % 整数时延近似
        s_del = circshift(s_cp, l_int);
        
        % 相位项 (t-l) 对应矩阵顺序 Pi * Delta
        phase = exp(1j * 2 * pi * k * (t_vec - l) / MN);
        
        r_clean = r_clean + h * s_del .* phase;
    end
    
    noise = (randn(len_tot, 1) + 1j*randn(len_tot, 1)) * sqrt(sigma_w2/2);
    r_cp = r_clean + noise;
end
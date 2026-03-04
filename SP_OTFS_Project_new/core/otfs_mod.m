function [s_cp, s_body, cp_len] = otfs_mod(X, cfg)
    % OTFS_MOD Reduced OTFS 调制
    % 输出:
    %   s_cp: 加CP后的信号
    %   s_body: 时域本体信号
    %   cp_len: 实际使用的CP长度

    N = cfg.grid.N;
    MN = numel(X);
    
    % Reduced OTFS: ISFFT along Doppler
    F_N_H = dftmtx(N)' / sqrt(N);
    S = X * F_N_H;
    s_body = S(:);
    
    % CP 长度逻辑
    if ~isempty(cfg.cp.len)
        cp_len = cfg.cp.len;
    else
        % 自动回退: max_search + guard
        cp_len = cfg.search.l_max + cfg.cp.guard;
    end
    
    % 加 CP
    s_cp = [s_body(MN-cp_len+1:end); s_body];
end
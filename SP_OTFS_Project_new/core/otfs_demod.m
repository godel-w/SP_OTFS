function [y_dd, y_time] = otfs_demod(r_cp, cfg, cp_len)
    % OTFS_DEMOD Reduced OTFS 解调
    % 输入: cp_len 必须由发送端传入，确保一致

    M = cfg.grid.M; 
    N = cfg.grid.N;
    MN = M * N;
    
    % 去 CP
    % 假设 r_cp 恰好是完整的一帧 (CP + Body)
    % 实际接收需做同步，此处假设已同步
    r_body = r_cp(cp_len+1 : cp_len+MN);
    y_time = r_body;
    
    % Reduced OTFS: SFFT along Doppler
    % Y = R * F_N
    R_mat = reshape(r_body, M, N);
    F_N = dftmtx(N) / sqrt(N);
    
    Y_mat = R_mat * F_N;
    y_dd = Y_mat(:);
end
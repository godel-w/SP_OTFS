function [bits_hat, stats] = spi_controller_lmmse(y_dd, x_p, cfg, bases, dict, est_type)
    sw2 = 10^(-cfg.sim.current_snr/10);
    sd2 = cfg.pwr.sigma_d2;
    MN = cfg.grid.M * cfg.grid.N;
    
    x_curr = zeros(MN, 1);
    sbl_state = [];
    Heff_curr = sparse(MN, MN);
    stats.sparsity_trace = [];
    
    for t = 1:cfg.alg.T_iter
        if t == 1
            y_in = y_dd;
            nv = sw2 + sd2;
        else
            % QAM 硬判决反馈
            % 1. 去除功率缩放，恢复到单位功率
            x_norm = x_curr / sqrt(sd2);
            
            % 2. 硬判决 (Slicing): Demod 得到整数 -> Mod 得到标准星座点
            %    这是找到最近星座点的最快方法
            sym_int = qamdemod(x_norm, cfg.mod.order, 'UnitAveragePower', true);
            sym_hard = qammod(sym_int, cfg.mod.order, 'UnitAveragePower', true);
            
            % 3. 恢复功率缩放
            x_hard = sym_hard * sqrt(sd2);
            
            % 4. 阻尼更新 (保持不变)
            x_fb = cfg.alg.beta * x_hard + (1-cfg.alg.beta) * x_curr;
            
            y_in = y_dd - Heff_curr * x_fb;
            nv = sw2 + cfg.alg.noise_floor_iter_factor * sd2;
        end
        
        if strcmp(est_type, 'SBL')
            [supp, h, sbl_state] = estimator_sbl(y_in, dict.Omega, cfg, nv, sbl_state);
        else
            [supp, h] = estimator_omp(y_in, dict.Omega, cfg, nv);
        end
        
        Heff_curr = build_heff_from_support(supp, h, dict.grid, bases, cfg);
        % [~, x_curr] = lmmse_detect(y_dd, x_p, Heff_curr, sw2, sd2);
        [~, x_curr] = lmmse_detect(y_dd, x_p, Heff_curr, sw2, sd2, cfg);
        stats.sparsity_trace(end+1) = length(supp);
    end
    
    % bits_hat = real(x_curr) < 0;
    % stats.h_est = h; stats.supp = supp;
    [bits_hat, ~] = lmmse_detect(y_dd, x_p, Heff_curr, sw2, sd2, cfg);
    stats.h_est = h; stats.supp = supp;
end
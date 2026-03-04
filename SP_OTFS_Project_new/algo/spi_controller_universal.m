function [bits_hat, x_soft, stats] = spi_controller_universal(y_dd, x_p, cfg, bases, dict, est_type, det_type)
    % det_type: 'MP' or 'LMMSE'
    sw2 = 10^(-cfg.sim.current_snr/10);
    sd2 = cfg.pwr.sigma_d2;
    MN = cfg.grid.M * cfg.grid.N;
    
    x_curr = zeros(MN, 1);
    sbl_state = [];
    mp_state = []; 
    Heff_curr = sparse(MN, MN);
    
    tic; % 开始计时该检测流程（包含信道估计+检测）
    for t = 1:cfg.alg.T_iter
        % --- 1. SIC 干扰消除 ---
        if t == 1
            y_in = y_dd; nv = sw2 + sd2;
        else
            x_norm = x_curr / sqrt(sd2);
            sym_int = qamdemod(x_norm, cfg.mod.order, 'UnitAveragePower', true);
            x_hard = qammod(sym_int, cfg.mod.order, 'UnitAveragePower', true) * sqrt(sd2);
            x_fb = cfg.alg.beta * x_hard + (1-cfg.alg.beta) * x_curr;
            y_in = y_dd - Heff_curr * x_fb;
            nv = sw2 + cfg.alg.noise_floor_iter_factor * sd2;
        end
        
        % --- 2. 稀疏信道估计 ---
        if strcmp(est_type, 'SBL')
            [supp, h, sbl_state] = estimator_sbl(y_in, dict.Omega, cfg, nv, sbl_state);
        else
            [supp, h] = estimator_omp(y_in, dict.Omega, cfg, nv);
        end
        
        % --- 3. 重构 Heff ---
        Heff_curr = build_heff_from_support(supp, h, dict.grid, bases, cfg);
        
        % --- 4. 符号检测 ---
        if strcmp(det_type, 'MP')
            [~, x_curr, mp_state] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);
        else
            [~, x_curr] = lmmse_detect(y_dd, x_p, Heff_curr, sw2, sd2, cfg);
        end
    end
    stats.time = toc; % 记录总执行时间
    
    % 最终判决
    if strcmp(det_type, 'MP')
        [bits_hat, x_soft, ~] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);
    else
        [bits_hat, x_soft] = lmmse_detect(y_dd, x_p, Heff_curr, sw2, sd2, cfg);
    end
    stats.h_est = h; stats.supp = supp;
end
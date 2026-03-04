function [bits_hat, stats] = spi_controller_mp(y_dd, x_p, cfg, bases, dict, est_type)
    sw2 = 10^(-cfg.sim.current_snr/10);
    sd2 = cfg.pwr.sigma_d2;
    MN = cfg.grid.M * cfg.grid.N;
    
    x_curr = zeros(MN, 1);
    
    % States for warm start
    sbl_state = [];
    mp_state = []; % [新增] MP 检测器状态
    
    Heff_curr = sparse(MN, MN);
    stats.sparsity_trace = [];
    
    for t = 1:cfg.alg.T_iter
        % 1. 干扰消除 (SIC) / 输入准备
        if t == 1
            y_in = y_dd;
            nv = sw2 + sd2;
        else
            % QAM 硬判决反馈 (解调再调制以获得干净的星座点)
            x_norm = x_curr / sqrt(sd2);
            sym_int = qamdemod(x_norm, cfg.mod.order, 'UnitAveragePower', true);
            sym_hard = qammod(sym_int, cfg.mod.order, 'UnitAveragePower', true);
            x_hard = sym_hard * sqrt(sd2);
            
            x_fb = cfg.alg.beta * x_hard + (1-cfg.alg.beta) * x_curr;
            y_in = y_dd - Heff_curr * x_fb;
            nv = sw2 + cfg.alg.noise_floor_iter_factor * sd2;
        end
        
        % 2. 稀疏估计
        if strcmp(est_type, 'SBL')
            [supp, h, sbl_state] = estimator_sbl(y_in, dict.Omega, cfg, nv, sbl_state);
        else
            [supp, h] = estimator_omp(y_in, dict.Omega, cfg, nv);
        end
        
        % 3. 重构 Heff
        Heff_curr = build_heff_from_support(supp, h, dict.grid, bases, cfg);
        
        % 4. MP 检测 (替换 LMMSE)
        % [注意] 仅支持 4QAM，若 cfg.mod.order!=4，mp_detect_qam 内部会报错
        [~, x_curr, mp_state] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);
        
        % Log
        stats.sparsity_trace(end+1) = length(supp);
    end
    
    % 最终硬判决输出
    [bits_hat, ~, ~] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);
    
    stats.h_est = h; 
    stats.supp = supp;
end
% SP-OTFS Main Execution Script (Updated for LMMSE vs MP Comparison)
addpath(genpath(pwd));
clc; clear; close all;

%% 1. Sanity Checks
if sanity_check_kron() && sanity_check_channel_alignment()
    disp('>>> All System Checks Passed.');
else
    error('>>> System Checks Failed!');
end

%% 2. Config & Init
% Select Case:
cfg = cfg_case_uwa_like_M32N16(); 

seed_rng(cfg.sim.rng_seed);
[bases.Btx, bases.Brx] = build_Btx_Brx(cfg.grid.M, cfg.grid.N);
MN = cfg.grid.M * cfg.grid.N;

%% 3. Simulation Loop
res.snr = cfg.sim.snr_list;

% 初始化结果存储
res.ber_lmmse = zeros(size(res.snr));
res.ber_mp    = zeros(size(res.snr));
res.nmse_sbl  = zeros(size(res.snr)); % 记录信道估计精度
res.time_lmmse = zeros(size(res.snr));
res.time_mp    = zeros(size(res.snr));

fprintf('\nRunning: %s (%dx%d) | Comparing LMMSE vs MP\n', cfg.channel.model, cfg.grid.M, cfg.grid.N);

for i = 1:length(res.snr)
    cfg.sim.current_snr = res.snr(i);
    sw2 = 10^(-res.snr(i)/10);
    
    % 累加器结构体
    % el: error lmmse, em: error mp, t_l: time lmmse, t_m: time mp
    acc = struct('el',0, 'em',0, 'ms',0, 'tl',0, 'tm',0, 'b',0);
    
    for f = 1:cfg.sim.n_frames
        % --- Source Generation ---
        n_bits = MN * cfg.mod.bps;
        bits_tx = randi([0, 1], n_bits, 1);
        syms_int = bit2int(bits_tx, cfg.mod.bps);
        syms_d = qammod(syms_int, cfg.mod.order, 'UnitAveragePower', true, 'InputType', 'integer');
        
        x_d = syms_d * sqrt(cfg.pwr.sigma_d2);
        Xd = reshape(x_d, cfg.grid.M, cfg.grid.N);
        
        % 导频生成
        p_bits = randi([0,3], MN, 1);
        xp = ((2*bitget(p_bits,1)-1)+1j*(2*bitget(p_bits,2)-1))/sqrt(2)*sqrt(cfg.pwr.sigma_p2);
        Xp = reshape(xp, cfg.grid.M, cfg.grid.N);
        
        % --- Channel & Noise ---
        [s_cp, ~, cp_len] = otfs_mod(Xd+Xp, cfg);
        chan = channel_factory(cfg);
        r_cp = apply_channel_phys(s_cp, chan, sw2, cfg);
        
        % --- Receiver ---
        [y_dd, ~] = otfs_demod(r_cp, cfg, cp_len);
        
        % Dictionary
        grid = grid_generate(cfg);
        Omega = build_dictionary_dd(xp, grid, bases, cfg);
        dict.Omega = Omega; dict.grid = grid;
        
        % --- 方案 A: SBL + LMMSE ---
        % 我们使用之前定义的控制器，增加计时
        tic;
        [bl, sts_l] = spi_controller_lmmse(y_dd, xp, cfg, bases, dict, 'SBL');
        acc.tl = acc.tl + toc;
        
        % --- 方案 B: SBL + MP (AMP) ---
        tic;
        [bm, sts_m] = spi_controller_mp(y_dd, xp, cfg, bases, dict, 'SBL');
        acc.tm = acc.tm + toc;

        % --- Stats Accumulation ---
        acc.el = acc.el + sum(bl ~= bits_tx);
        acc.em = acc.em + sum(bm ~= bits_tx);
        
        % 记录 SBL 的信道估计 MSE (两者使用相同的信道估计器，统计一个即可)
        acc.ms = acc.ms + nmse_calc(sts_m.h_est, sts_m.supp, grid, chan);
        
        acc.b = acc.b + n_bits;
    end
    
    % --- Store Averages ---
    res.ber_lmmse(i) = acc.el / acc.b;
    res.ber_mp(i)    = acc.em / acc.b;
    res.nmse_sbl(i)  = acc.ms / cfg.sim.n_frames;
    res.time_lmmse(i) = acc.tl / cfg.sim.n_frames;
    res.time_mp(i)    = acc.tm / cfg.sim.n_frames;
    
    fprintf('SNR %2d | BER: LMMSE=%.4f, MP=%.4f | MSE: %.3f | Time: L=%.3fs, M=%.3fs\n', ...
        res.snr(i), res.ber_lmmse(i), res.ber_mp(i), res.nmse_sbl(i), res.time_lmmse(i), res.time_mp(i));
end

%% 4. Plotting Results
figure('Name', 'OTFS Detector Comparison');
subplot(2,1,1);
semilogy(res.snr, res.ber_lmmse, 'r-o', 'LineWidth', 2, 'DisplayName', 'SBL + LMMSE'); hold on;
semilogy(res.snr, res.ber_mp, 'b-s', 'LineWidth', 2, 'DisplayName', 'SBL + MP (AMP)');
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('Bit Error Rate Performance');
legend show;

subplot(2,1,2);
plot(res.snr, res.time_lmmse, 'r--o', 'LineWidth', 1.5, 'DisplayName', 'LMMSE Time'); hold on;
plot(res.snr, res.time_mp, 'b--s', 'LineWidth', 1.5, 'DisplayName', 'MP Time');
grid on; xlabel('SNR (dB)'); ylabel('Average Time per Frame (s)');
title('Computational Complexity (Execution Time)');
legend show;

% 注意：如果你的 plot_results 函数不支持新字段，请使用上面自定义的绘图代码
% plot_results(res, cfg);
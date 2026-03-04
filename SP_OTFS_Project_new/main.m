% SP-OTFS Main Execution Script
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
% cfg = cfg_case_inttaps_M16N16();
cfg = cfg_case_uwa_like_M32N16();

seed_rng(cfg.sim.rng_seed);
[bases.Btx, bases.Brx] = build_Btx_Brx(cfg.grid.M, cfg.grid.N);
MN = cfg.grid.M * cfg.grid.N;

%% 3. Simulation Loop
res.snr = cfg.sim.snr_list;
res.ber_omp = zeros(size(res.snr)); res.ber_sbl = res.ber_omp;
res.nmse_omp = zeros(size(res.snr)); res.nmse_sbl = res.nmse_omp;
res.spar_omp = zeros(size(res.snr)); res.spar_sbl = res.spar_omp;

fprintf('\nRunning: %s (%dx%d)\n', cfg.channel.model, cfg.grid.M, cfg.grid.N);

for i = 1:length(res.snr)
    cfg.sim.current_snr = res.snr(i);
    sw2 = 10^(-res.snr(i)/10);
    
    acc = struct('eo',0,'es',0,'mo',0,'ms',0,'so',0,'ss',0,'b',0);
    
    for f = 1:cfg.sim.n_frames
        % Source
        % 1. 生成比特: 总比特数 = 网格点数 * 每个符号的比特数
        n_bits = MN * cfg.mod.bps;
        bits_tx = randi([0, 1], n_bits, 1);

        % 2. 映射为 QAM 符号 (DD域)
        % bit2int 将比特流转为整数索引
        syms_int = bit2int(bits_tx, cfg.mod.bps);
        % qammod 生成符号，并强制归一化功率为 1 (UnitAveragePower=true)
        syms_d = qammod(syms_int, cfg.mod.order, 'UnitAveragePower', true, 'InputType', 'integer');

        % 3. 功率缩放 (乘以数据功率 sigma_d2)
        x_d = syms_d * sqrt(cfg.pwr.sigma_d2);
        Xd = reshape(x_d, cfg.grid.M, cfg.grid.N);
        p_bits = randi([0,3], MN, 1);
        xp = ((2*bitget(p_bits,1)-1)+1j*(2*bitget(p_bits,2)-1))/sqrt(2)*sqrt(cfg.pwr.sigma_p2);
        Xp = reshape(xp, cfg.grid.M, cfg.grid.N);
        
        % Mod & Channel
        [s_cp, ~, cp_len] = otfs_mod(Xd+Xp, cfg);
        chan = channel_factory(cfg);
        r_cp = apply_channel_phys(s_cp, chan, sw2, cfg);
        
        % Demod
        [y_dd, ~] = otfs_demod(r_cp, cfg, cp_len);
        
        % Dict
        grid = grid_generate(cfg);
        Omega = build_dictionary_dd(xp, grid, bases, cfg);
        dict.Omega = Omega; dict.grid = grid;
        
        % Run
        % tic;
        [bo, sto] = spi_controller_mp(y_dd, xp, cfg, bases, dict, 'OMP');
        % time_omp = toc;
        
        % tic;
        [bs, sts] = spi_controller_mp(y_dd, xp, cfg, bases, dict, 'SBL');
        % time_sbl = toc;

        % Stats
        acc.eo = acc.eo + sum(bo~=bits_tx); acc.es = acc.es + sum(bs~=bits_tx);
        acc.mo = acc.mo + nmse_calc(sto.h_est, sto.supp, grid, chan);
        acc.ms = acc.ms + nmse_calc(sts.h_est, sts.supp, grid, chan);
        acc.so = acc.so + length(sto.supp); acc.ss = acc.ss + length(sts.supp);
        % acc.b = acc.b + MN;
        acc.b = acc.b + MN * cfg.mod.bps;
    end
    
    % Store
    res.ber_omp(i) = acc.eo/acc.b; res.ber_sbl(i) = acc.es/acc.b;
    res.nmse_omp(i) = acc.mo/cfg.sim.n_frames; res.nmse_sbl(i) = acc.ms/cfg.sim.n_frames;
    res.spar_omp(i) = acc.so/cfg.sim.n_frames; res.spar_sbl(i) = acc.ss/cfg.sim.n_frames;
    
    fprintf('SNR %2d | BER: O=%.4f S=%.4f | MSE: O=%.3f S=%.3f\n', ...
        res.snr(i), res.ber_omp(i), res.ber_sbl(i), res.nmse_omp(i), res.nmse_sbl(i));
end

plot_results(res, cfg);
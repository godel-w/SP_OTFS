% SP-OTFS Main Execution Script (Enhanced for 6 Figures)
addpath(genpath(pwd));
clc; clear; close all;

%% 1. Sanity Checks
if sanity_check_kron() && sanity_check_channel_alignment()
    disp('>>> All System Checks Passed.');
else
    error('>>> System Checks Failed!');
end

%% 2. Config & Init
cfg = cfg_case_uwa_like_M32N16();   % 你现在用的配置
seed_rng(cfg.sim.rng_seed);

[bases.Btx, bases.Brx] = build_Btx_Brx(cfg.grid.M, cfg.grid.N);
MN = cfg.grid.M * cfg.grid.N;

%% 2.1 额外：选择一组“抓图用”的参考点（用于图3/4/5）
% 选一个中等 SNR 和第 1 帧做快照（可自行改）
snr_ref = cfg.sim.snr_list(round(numel(cfg.sim.snr_list)/2));
frame_ref = 1;

%% 3. Simulation Loop
res.snr = cfg.sim.snr_list;

res.ber_omp  = zeros(size(res.snr));
res.ber_sbl  = zeros(size(res.snr));
res.nmse_omp = zeros(size(res.snr));
res.nmse_sbl = zeros(size(res.snr));
res.spar_omp = zeros(size(res.snr));
res.spar_sbl = zeros(size(res.snr));

% 新增：耗时（图6）
res.time_omp = zeros(size(res.snr));   % 平均每帧耗时
res.time_sbl = zeros(size(res.snr));

% 新增：保存一条代表性的迭代 trace（图3/图5）
res.trace_omp = cell(size(res.snr));
res.trace_sbl = cell(size(res.snr));

% 新增：保存一帧 taps 快照（图4）
res.snap = struct(); % res.snap.truth / res.snap.omp / res.snap.sbl / res.snap.grid

fprintf('\nRunning: %s (%dx%d)\n', cfg.channel.model, cfg.grid.M, cfg.grid.N);

for i = 1:length(res.snr)

    cfg.sim.current_snr = res.snr(i);
    sw2 = 10^(-res.snr(i)/10);

    acc = struct('eo',0,'es',0,'mo',0,'ms',0,'so',0,'ss',0,'b',0, ...
                 'to',0,'ts',0);  % 新增 to/ts 做耗时累计

    for f = 1:cfg.sim.n_frames

        % ====== Source ======
        n_bits = MN * cfg.mod.bps;
        bits_tx = randi([0, 1], n_bits, 1);

        syms_int = bit2int(bits_tx, cfg.mod.bps);
        syms_d = qammod(syms_int, cfg.mod.order, 'UnitAveragePower', true, 'InputType', 'integer');

        x_d = syms_d * sqrt(cfg.pwr.sigma_d2);
        Xd  = reshape(x_d, cfg.grid.M, cfg.grid.N);

        p_bits = randi([0,3], MN, 1);
        xp = ((2*bitget(p_bits,1)-1)+1j*(2*bitget(p_bits,2)-1))/sqrt(2)*sqrt(cfg.pwr.sigma_p2);
        Xp = reshape(xp, cfg.grid.M, cfg.grid.N);

        % ====== Mod & Channel ======
        [s_cp, ~, cp_len] = otfs_mod(Xd+Xp, cfg);
        chan = channel_factory(cfg);
        r_cp = apply_channel_phys(s_cp, chan, sw2, cfg);

        % ====== Demod ======
        [y_dd, ~] = otfs_demod(r_cp, cfg, cp_len);

        % ====== Dict ======
        grid = grid_generate(cfg);
        Omega = build_dictionary_dd(xp, grid, bases, cfg);
        dict.Omega = Omega;
        dict.grid = grid;

        % ====== Run OMP ======
        t0 = tic;
        [bo, sto] = spi_controller_mp_new(y_dd, xp, cfg, bases, dict, 'OMP');
        acc.to = acc.to + toc(t0);

        % ====== Run SBL ======
        t1 = tic;
        [bs, sts] = spi_controller_mp_new(y_dd, xp, cfg, bases, dict, 'SBL');
        acc.ts = acc.ts + toc(t1);

        % ====== Stats ======
        acc.eo = acc.eo + sum(bo~=bits_tx);
        acc.es = acc.es + sum(bs~=bits_tx);

        acc.mo = acc.mo + nmse_calc(sto.h_est, sto.supp, grid, chan);
        acc.ms = acc.ms + nmse_calc(sts.h_est, sts.supp, grid, chan);

        acc.so = acc.so + length(sto.supp);
        acc.ss = acc.ss + length(sts.supp);

        acc.b = acc.b + MN * cfg.mod.bps;

        % ====== 抓一帧：用于图3/图4/图5 ======
        if (res.snr(i) == snr_ref) && (f == frame_ref)
            % 1) 迭代 trace（图3、图5）
            if isfield(sto,'trace'); res.trace_omp{i} = sto.trace; end
            if isfield(sts,'trace'); res.trace_sbl{i} = sts.trace; end

            % 2) taps 快照（图4）
            % 真值 taps 从 chan 中取（你的 channel_factory/apply_channel_phys 里通常有 delays/dopplers/gains）
            res.snap.truth = chan;     % 直接存 chan（后面绘图函数再解析）
            res.snap.omp   = sto;
            res.snap.sbl   = sts;
            res.snap.grid  = grid;
        end
    end

    % ====== Store ======
    res.ber_omp(i)  = acc.eo/acc.b;
    res.ber_sbl(i)  = acc.es/acc.b;

    res.nmse_omp(i) = acc.mo/cfg.sim.n_frames;
    res.nmse_sbl(i) = acc.ms/cfg.sim.n_frames;

    res.spar_omp(i) = acc.so/cfg.sim.n_frames;
    res.spar_sbl(i) = acc.ss/cfg.sim.n_frames;

    res.time_omp(i) = acc.to/cfg.sim.n_frames;
    res.time_sbl(i) = acc.ts/cfg.sim.n_frames;

    fprintf('SNR %2d | BER: O=%.4f S=%.4f | NMSE: O=%.3f S=%.3f | t(s): O=%.3f S=%.3f\n', ...
        res.snr(i), res.ber_omp(i), res.ber_sbl(i), res.nmse_omp(i), res.nmse_sbl(i), ...
        res.time_omp(i), res.time_sbl(i));
end

plot_results6(res, cfg);  % 你原有的函数：下面我给你一个增强版实现

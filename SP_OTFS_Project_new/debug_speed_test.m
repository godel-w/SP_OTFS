%% SP-OTFS Speed Test Script (Single Frame, Single SNR)
% 专门用于对比 OMP 和 SBL 的运行速度
addpath(genpath(pwd));
clc; clear; close all;

%% 1. 配置与初始化
fprintf('>>> 正在初始化配置...\n');

% 加载指定案例配置
cfg = cfg_case_uwa_like_M32N16(); 

% [强制修改]：只跑 1 帧，只跑 20dB (高SNR下SBL迭代次数较稳定)
target_snr = 10;
cfg.sim.snr_list = target_snr;
cfg.sim.n_frames = 1;
cfg.sim.current_snr = target_snr; % 关键：spi_controller 内部依赖这个变量计算噪声

% [强制修改]：确保调制参数存在 (防止配置文件没更新)
if ~isfield(cfg, 'mod')
    cfg.mod.order = 4; 
    cfg.mod.bps = 2; 
end

% 初始化随机种和基矩阵
seed_rng(cfg.sim.rng_seed);
[bases.Btx, bases.Brx] = build_Btx_Brx(cfg.grid.M, cfg.grid.N);
MN = cfg.grid.M * cfg.grid.N;

%% 2. 准备数据 (这部分不计入算法耗时)
fprintf('>>> 正在生成信号与构建字典 (公共开销)...\n');

% --- A. 发送端 (Source) ---
% 1. 生成比特
n_bits = MN * cfg.mod.bps;
bits_tx = randi([0, 1], n_bits, 1);

% 2. QAM 映射
syms_int = bit2int(bits_tx, cfg.mod.bps);
syms_d = qammod(syms_int, cfg.mod.order, 'UnitAveragePower', true, 'InputType', 'integer');
x_d = syms_d * sqrt(cfg.pwr.sigma_d2);
Xd = reshape(x_d, cfg.grid.M, cfg.grid.N);

% 3. 导频
p_bits = randi([0,3], MN, 1);
xp = ((2*bitget(p_bits,1)-1)+1j*(2*bitget(p_bits,2)-1))/sqrt(2)*sqrt(cfg.pwr.sigma_p2);
Xp = reshape(xp, cfg.grid.M, cfg.grid.N);

% --- B. 调制与信道 (Mod & Channel) ---
[s_cp, ~, cp_len] = otfs_mod(Xd+Xp, cfg);
chan = channel_factory(cfg);
sw2 = 10^(-target_snr/10);
r_cp = apply_channel_phys(s_cp, chan, sw2, cfg);

% --- C. 解调 (Demod) ---
[y_dd, ~] = otfs_demod(r_cp, cfg, cp_len);

% --- D. 构建字典 (Dictionary) ---
% 字典构建通常只做一次，不算在迭代检测算法的耗时里
grid = grid_generate(cfg);
Omega = build_dictionary_dd(xp, grid, bases, cfg);
dict.Omega = Omega; dict.grid = grid;

fprintf('>>> 数据准备完毕。网格: %dx%d, 字典列数: %d\n', ...
    cfg.grid.M, cfg.grid.N, size(Omega, 2));

%% 3. JIT 预热 (Warm-up)
% MATLAB 首次运行函数会进行编译，这会影响计时准确性。
% 我们先“空跑”一次，不输出结果，也不计时。
fprintf('>>> 正在预热 JIT 编译器 (Warm-up)... \n');
[~, ~] = spi_controller(y_dd, xp, cfg, bases, dict, 'OMP');
[~, ~] = spi_controller(y_dd, xp, cfg, bases, dict, 'SBL');
fprintf('>>> 预热完成。\n');

%% 4. 正式计时对比 (Speed Test)
fprintf('\n========================================\n');
fprintf('       开始速度测试 (T_iter = %d)\n', cfg.alg.T_iter);
fprintf('========================================\n');

% --- 测试 OMP ---
tic;
[bo, sto] = spi_controller(y_dd, xp, cfg, bases, dict, 'OMP');
time_omp = toc;

ber_omp = sum(bo ~= bits_tx) / n_bits;
fprintf('OMP 耗时: %.6f 秒 | 稀疏度: %d | BER: %.4f\n', ...
    time_omp, length(sto.supp), ber_omp);

% --- 测试 SBL ---
tic;
[bs, sts] = spi_controller(y_dd, xp, cfg, bases, dict, 'SBL');
time_sbl = toc;

ber_sbl = sum(bs ~= bits_tx) / n_bits;
fprintf('SBL 耗时: %.6f 秒 | 稀疏度: %d | BER: %.4f\n', ...
    time_sbl, length(sts.supp), ber_sbl);

%% 5. 结果分析
ratio = time_sbl / time_omp;
fprintf('\n>>> 结论: 在此参数下，SBL 比 OMP 慢 %.2f 倍\n', ratio);

% 可选：画出检测结果对比
% figure; 
% subplot(2,1,1); stem(bits_tx(1:100)); hold on; stem(bo(1:100), 'rx'); title('OMP bits (first 100)');
% subplot(2,1,2); stem(bits_tx(1:100)); hold on; stem(bs(1:100), 'rx'); title('SBL bits (first 100)');
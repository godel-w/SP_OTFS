function plot_results6(res, cfg)

snr = res.snr;

%% 图1：BER vs SNR
figure; 
semilogy(snr, res.ber_omp, '-o'); hold on;
semilogy(snr, res.ber_sbl, '-s');
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('OMP+MP','SBL+MP','Location','southwest');
title('BER vs SNR');

%% 图2：NMSE vs SNR
figure;
semilogy(snr, res.nmse_omp, '-o'); hold on;
semilogy(snr, res.nmse_sbl, '-s');
grid on; xlabel('SNR (dB)'); ylabel('NMSE');
legend('OMP','SBL','Location','southwest');
title('NMSE vs SNR');

%% 图6：耗时 vs SNR（每帧平均）
figure;
plot(snr, res.time_omp, '-o'); hold on;
plot(snr, res.time_sbl, '-s');
grid on; xlabel('SNR (dB)'); ylabel('Time per frame (s)');
legend('OMP pipeline','SBL pipeline','Location','northwest');
title('Runtime vs SNR');

%% 图3 & 图5：迭代收敛 + 导频污染缓解（取 snap 里存的那条 trace）
% 找到非空 trace 的一个索引
idx = find(~cellfun(@isempty, res.trace_sbl), 1, 'first');
if ~isempty(idx)
    trS = res.trace_sbl{idx};
    trO = res.trace_omp{idx};

    % 图3：收敛（xchg vs iter）
    figure;
    plot(trO.xchg, '-o'); hold on;
    plot(trS.xchg, '-s');
    grid on; xlabel('Outer iteration t'); ylabel('||x^{t}-x^{t-1}|| / ||x^{t-1}||');
    legend('OMP pipeline','SBL pipeline','Location','northeast');
    title('Convergence of outer loop (symbol change)');

    % 图5：导频残差能量下降（yin_res_pilot vs iter）
    figure;
    semilogy(trO.yin_res_pilot, '-o'); hold on;
    semilogy(trS.yin_res_pilot, '-s');
    grid on; xlabel('Outer iteration t'); ylabel('||y_{in}-H x_p||^2 / MN');
    legend('OMP pipeline','SBL pipeline','Location','northeast');
    title('Pilot contamination mitigation evidence (pilot residual)');
else
    warning('No trace found. Please ensure spi_controller_mp outputs stats.trace');
end

%% 图4：DD taps 可视化（真值 vs 估计）
if isfield(res,'snap') && isfield(res.snap,'grid')
    plot_dd_taps_snapshot(res.snap, cfg);
else
    warning('No snapshot found. Please ensure res.snap is captured in main loop.');
end

end

function plot_results(res, cfg)
    figure('Color','w','Position',[100,100,1200,400]);
    lbl = cfg.channel.model;
    
    subplot(1,3,1); semilogy(res.snr, res.ber_omp, '--o'); hold on;
    semilogy(res.snr, res.ber_sbl, '-s'); grid on; legend('OMP','SBL');
    xlabel('SNR'); ylabel('BER'); title(['BER: ' lbl]);
    
    subplot(1,3,2); semilogy(res.snr, res.nmse_omp, '--o'); hold on;
    semilogy(res.snr, res.nmse_sbl, '-s'); grid on;
    xlabel('SNR'); ylabel('NMSE'); title('Channel Est');
    
    subplot(1,3,3); plot(res.snr, res.spar_omp, '--o'); hold on;
    plot(res.snr, res.spar_sbl, '-s'); yline(cfg.channel.n_paths,'k:');
    xlabel('SNR'); ylabel('Sparsity'); grid on; title('Avg Taps');
end
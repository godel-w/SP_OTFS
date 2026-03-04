function plot_dd_taps_snapshot(snap, cfg)
%PLOT_DD_TAPS_SNAPSHOT  Visualize DD taps: ground truth vs OMP vs SBL
%
% snap.truth : channel struct (contains delays/dopplers/gains or similar)
% snap.omp   : stats from spi_controller_mp (fields: supp, h_est)
% snap.sbl   : stats from spi_controller_mp (fields: supp, h_est)
% snap.grid  : grid struct with mapping arrays (fields: l, k)
%
% cfg        : config struct (for title labels)

    ddgrid = snap.grid;   % <-- avoid name 'grid' (MATLAB has grid() command)

    % ---------- Parse ground-truth taps from channel struct ----------
    chan = snap.truth;

    % delays
    if isfield(chan,'delays')
        l_true = chan.delays(:);
    elseif isfield(chan,'path_delays')
        l_true = chan.path_delays(:);
    else
        l_true = [];
    end

    % dopplers
    if isfield(chan,'dopplers')
        k_true = chan.dopplers(:);
    elseif isfield(chan,'path_dopplers')
        k_true = chan.path_dopplers(:);
    else
        k_true = [];
    end

    % gains
    if isfield(chan,'gains')
        h_true = chan.gains(:);
    elseif isfield(chan,'path_gains')
        h_true = chan.path_gains(:);
    elseif isfield(chan,'gains_dB')
        h_true = 10.^(chan.gains_dB(:)/20);
    elseif isfield(chan,'path_gains_dB')
        h_true = 10.^(chan.path_gains_dB(:)/20);
    else
        h_true = [];
    end

    % ---------- Estimated taps: map supp -> (l,k) and magnitudes ----------
    [l_o, k_o, mag_o] = supp_to_lk_mag(snap.omp.supp, snap.omp.h_est, ddgrid);
    [l_s, k_s, mag_s] = supp_to_lk_mag(snap.sbl.supp, snap.sbl.h_est, ddgrid);

    % ---------- Plot ----------
    figure;

    subplot(1,3,1);
    if isempty(l_true)
        text(0.1,0.5,'No ground-truth taps found in chan struct','FontSize',10);
        axis off;
    else
        stem3(l_true, k_true, abs(h_true), 'filled');
        xlabel('delay index l'); ylabel('doppler index k'); zlabel('|h|');
        grid(gca,'on');
    end
    title('Ground truth taps');

    subplot(1,3,2);
    if isempty(l_o)
        text(0.1,0.5,'No OMP taps','FontSize',10); axis off;
    else
        stem3(l_o, k_o, mag_o, 'filled');
        xlabel('delay index l'); ylabel('doppler index k'); zlabel('|h|');
        grid(gca,'on');
    end
    title('Estimated taps (OMP)');

    subplot(1,3,3);
    if isempty(l_s)
        text(0.1,0.5,'No SBL taps','FontSize',10); axis off;
    else
        stem3(l_s, k_s, mag_s, 'filled');
        xlabel('delay index l'); ylabel('doppler index k'); zlabel('|h|');
        grid(gca,'on');
    end
    title('Estimated taps (SBL)');

    sgtitle(sprintf('DD taps visualization (%dx%d)', cfg.grid.M, cfg.grid.N));
end


function [l_vec, k_vec, mag_vec] = supp_to_lk_mag(supp, h_est, ddgrid)
% Map support indices to delay/doppler indices using ddgrid.l/ddgrid.k
    if isempty(supp)
        l_vec = []; k_vec = []; mag_vec = [];
        return;
    end

    % These should exist in your grid_generate() output
    l_vec  = ddgrid.l(supp);
    k_vec  = ddgrid.k(supp);
    mag_vec = abs(h_est(:));
end

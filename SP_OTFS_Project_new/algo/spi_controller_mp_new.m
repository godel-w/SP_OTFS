function [bits_hat, stats] = spi_controller_mp_new(y_dd, x_p, cfg, bases, dict, est_type)
% SPI_CONTROLLER_MP  SP-OTFS iterative receiver controller (SBL/OMP + MP)
%
% Inputs:
%   y_dd   : DD-domain received vector (MN x 1), contains pilot+data
%   x_p    : pilot vector (MN x 1)
%   cfg    : config struct
%   bases  : bases.Btx, bases.Brx
%   dict   : dict.Omega, dict.grid
%   est_type : 'SBL' or 'OMP'
%
% Outputs:
%   bits_hat : detected bits (MN*bps x 1)
%   stats    : struct with fields:
%              - sparsity_trace
%              - trace (for plotting convergence & pilot contamination mitigation)
%              - h_est, supp

    % --- Basic params ---
    sw2 = 10^(-cfg.sim.current_snr/10);
    sd2 = cfg.pwr.sigma_d2;
    MN  = cfg.grid.M * cfg.grid.N;

    % --- Init ---
    x_curr = zeros(MN, 1);

    % Warm start states
    sbl_state = [];
    mp_state  = [];

    Heff_curr = sparse(MN, MN);

    % --- Stats containers ---
    T = cfg.alg.T_iter;
    stats.sparsity_trace = zeros(T,1);

    % Enhanced trace for Figures (3) and (5)
    trace.supp_len       = zeros(T,1);   % |supp|
    trace.nv             = zeros(T,1);   % effective noise var used in estimator
    trace.xchg           = nan(T,1);     % ||x^t-x^{t-1}||/||x^{t-1}||
    trace.yin_res_pilot  = zeros(T,1);   % ||y_in - H*x_p||^2 / MN
    trace.yin_norm       = zeros(T,1);   % ||y_in||^2 / MN

    x_prev = [];

    % =========================
    % Outer Loop
    % =========================
    for t = 1:T

        % ---------------------------------------------------------
        % 1) SIC / Estimator input preparation (data cancellation)
        % ---------------------------------------------------------
        if t == 1
            y_in = y_dd;
            nv   = sw2 + sd2;   % conservative: treat data as extra noise
        else
            % QAM hard decision feedback: demod -> remod to get clean constellation
            x_norm  = x_curr / sqrt(sd2);
            sym_int = qamdemod(x_norm, cfg.mod.order, 'UnitAveragePower', true);
            sym_hard = qammod(sym_int, cfg.mod.order, 'UnitAveragePower', true);
            x_hard   = sym_hard * sqrt(sd2);

            % feedback damping
            x_fb = cfg.alg.beta * x_hard + (1 - cfg.alg.beta) * x_curr;

            % cancel estimated data component at estimator input
            y_in = y_dd - Heff_curr * x_fb;

            % effective noise floor (to avoid overconfidence / error propagation)
            nv = sw2 + cfg.alg.noise_floor_iter_factor * sd2;
        end

        % record nv and input energy
        trace.nv(t)       = nv;
        trace.yin_norm(t) = (norm(y_in).^2) / MN;

        % ----------------------------
        % 2) Sparse channel estimation
        % ----------------------------
        if strcmpi(est_type, 'SBL')
            [supp, h, sbl_state] = estimator_sbl(y_in, dict.Omega, cfg, nv, sbl_state);
        else
            [supp, h] = estimator_omp(y_in, dict.Omega, cfg, nv);
        end

        % ----------------------------
        % 3) Rebuild Heff from support
        % ----------------------------
        Heff_curr = build_heff_from_support(supp, h, dict.grid, bases, cfg);

        % ----------------------------
        % 4) MP/AMP detection (pilot cancellation inside)
        % ----------------------------
        % Note: mp_detect_qam supports QPSK/4QAM only; it will error for other orders.
        [~, x_curr, mp_state] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);

        % ----------------------------
        % Logging for plots
        % ----------------------------
        trace.supp_len(t)  = length(supp);
        trace.yin_res_pilot(t) = (norm(y_in - Heff_curr * x_p).^2) / MN;

        if ~isempty(x_prev)
            trace.xchg(t) = norm(x_curr - x_prev) / (norm(x_prev) + 1e-12);
        end
        x_prev = x_curr;

        stats.sparsity_trace(t) = length(supp);
    end

    % =========================
    % Final hard decision output
    % =========================
    [bits_hat, ~, ~] = mp_detect_qam(y_dd, x_p, Heff_curr, sw2, sd2, cfg, mp_state);

    % Final stats
    stats.h_est = h;
    stats.supp  = supp;
    stats.trace = trace;
end

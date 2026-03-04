function [bits_hat, x_soft, state] = mp_detect_qam(y_dd, x_p, Heff, sw2, sd2, cfg, state)
% MP_DETECT_QAM Message Passing Detector for OTFS (AMP Algorithm)
% Currently supports 4QAM (QPSK) only.
%
% Model: y_d = Heff * x_d + w, where y_d = y_dd - Heff * x_p
%
% Inputs:
%   y_dd : Received signal in DD domain (with pilot)
%   x_p  : Pilot vector
%   Heff : Effective channel matrix (Sparse)
%   sw2  : Noise variance (sigma_w^2)
%   sd2  : Data variance (sigma_d^2)
%   cfg  : Config struct (must contain cfg.mod.order, cfg.det.mp.*)
%   state: Previous state for warm start (optional)
%
% Outputs:
%   bits_hat: Hard decision bits [MN*bps x 1]
%   x_soft  : Soft symbol estimates
%   state   : Internal state for next iteration

    % 1. Sanity Check & Setup
    if cfg.mod.order ~= 4
        error('MP Detector currently only supports 4QAM (QPSK).');
    end
    
    MN = length(y_dd);
    
    % SIC: Subtract Pilot to get data observation
    y = y_dd - Heff * x_p;
    
    % Load Config or Defaults
    if isfield(cfg, 'det') && isfield(cfg.det, 'mp')
        MAX_ITER = cfg.det.mp.max_iter;
        DAMP     = cfg.det.mp.damp;
        TOL      = cfg.det.mp.tol;
        VMIN     = cfg.det.mp.vmin;
    else
        % Fallback defaults
        MAX_ITER = 30; DAMP = 0.6; TOL = 1e-4; VMIN = 1e-10;
    end
    
    % 2. Initialization / Warm Start
    if isempty(state)
        x_t = zeros(MN, 1);     % Initial estimate (Mean)
        z_t = y;                % Initial residual
        v_t = ones(MN, 1)*sd2;  % Initial variance (Prior variance)
    else
        x_t = state.x;
        z_t = state.z;
        v_t = state.v;
        % Re-calc residual to be safe with updated Heff
        % z_t = y - Heff * x_t; 
    end
    
    % Pre-compute Constellation (4QAM scaled by sqrt(sd2))
    % 4QAM: (+-1 +-1j) / sqrt(2) * sqrt(sd2)
    scale = sqrt(sd2);
    S = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2) * scale; 
    
    % Pre-compute H Hermitian (Sparse)
    Heff_H = Heff'; 
    
    % 3. AMP Iterations
    for iter = 1:MAX_ITER
        x_old = x_t;
        
        % --- Step A: Residual Update (with Onsager Correction) ---
        % b_t = (1/M) * sum( v_t / (tau_prev + sw2) )
        % Here we approximate effective noise variance from residual energy
        tau_z = mean(abs(z_t).^2); 
        if tau_z < VMIN, tau_z = VMIN; end
        
        % Onsager term (Standard complex AMP approximation)
        % b_t = (N/M) * mean(v_t) / (tau_z + sw2) ? 
        % Simplified Onsager for Unitary-like OTFS:
        b_t = mean(v_t) / (tau_z + sw2); 
        
        % z_new = y - H*x_t + z_old * b_t
        z_new = y - Heff * x_t + z_t * b_t;
        
        % Update Residual
        z_t = z_new;
        
        % --- Step B: Effective Observation (Matched Filter) ---
        % r = x_t + H' * z_t
        r = x_t + Heff_H * z_t;
        
        % Effective noise variance on r
        tau_r = mean(abs(z_t).^2); 
        if tau_r < VMIN, tau_r = VMIN; end
        
        % --- Step C: Discrete Denoising (4QAM) ---
        % Compute P(S_k | r_n) propto exp( -|r_n - S_k|^2 / tau_r )
        
        % Vectorized Distance Calculation: |r - s|^2
        % r is [MN x 1], S is [4 x 1] -> Result [MN x 4]
        dist_sq = abs(r - S.').^2; 
        
        % Exponents: -dist / tau_r
        exponents = -dist_sq ./ tau_r;
        
        % Log-Sum-Exp Trick for stability
        max_exp = max(exponents, [], 2);
        probs = exp(exponents - max_exp); % [MN x 4]
        sum_probs = sum(probs, 2);
        probs = probs ./ sum_probs;       % Normalize
        
        % Posterior Mean (Soft Symbol)
        x_post = probs * S;
        
        % Posterior Variance: E[|x|^2] - |E[x]|^2
        % E[|x|^2] = sum( |S|^2 * p ) = |S|^2 (since all |S| are same for 4QAM)
        % Actually |S_k|^2 = sd2 for all k in QPSK
        x_sq_mean = sum(probs .* (abs(S.').^2), 2); % Should be sd2
        v_post = x_sq_mean - abs(x_post).^2;
        
        % --- Step D: Damping & Update ---
        x_t = DAMP * x_post + (1-DAMP) * x_t;
        v_t = DAMP * v_post + (1-DAMP) * v_t;
        
        % --- Step E: Check Convergence ---
        diff_norm = norm(x_t - x_old) / (norm(x_old) + 1e-9);
        if diff_norm < TOL
            break;
        end
    end
    
    % 4. Output Generation
    x_soft = x_t;
    
    % Save state for next SP-I iteration
    state.x = x_t;
    state.z = z_t;
    state.v = v_t;
    
    % Hard Decision (Bit output)
    % Normalize power back to 1 for qamdemod
    x_norm = x_soft / sqrt(sd2);
    
    bits_hat = qamdemod(x_norm, cfg.mod.order, ...
        'UnitAveragePower', true, ...
        'OutputType', 'bit');
end
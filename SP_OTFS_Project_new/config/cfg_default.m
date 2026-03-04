function cfg = cfg_default()
    % Grid
    cfg.grid.M = 32;
    cfg.grid.N = 16;
    
    % CP
    cfg.cp.len = [];      % 为空则自动计算
    cfg.cp.guard = 2;     % 自动计算时的保护间隔
    
    % Power
    cfg.pwr.sigma_p2 = 0.3;
    cfg.pwr.sigma_d2 = 1 - 0.3;
    
    % Simulation
    cfg.sim.snr_list = 0:4:24;
    cfg.sim.n_frames = 10;
    cfg.sim.rng_seed = 42;
    
    % Modulation Parameters
    cfg.mod.order = 4; % 4 = QPSK, 16 = 16QAM, 64 = 64QAM
    cfg.mod.bps = log2(cfg.mod.order); % bit per symbol
    
    % Channel Defaults
    cfg.channel.model = 'INT_TAPS';
    cfg.channel.n_paths = 4;
    cfg.channel.max_delay = 5;
    cfg.channel.max_doppler = 2;
    
    % Search Defaults
    cfg.search.l_max = 6;
    cfg.search.k_max = 3;
    
    % Algo
    cfg.alg.T_iter = 3;
    cfg.alg.beta = 0.6;
    cfg.alg.noise_floor_iter_factor = 0.05;
    
    % Estimators
    cfg.est.omp.max_sparsity = 10;
    cfg.est.omp.resid_factor = 1.1;
    
    cfg.est.sbl.max_iter = 15;
    cfg.est.sbl.tol = 1e-4;
    cfg.est.sbl.prune = 1e-3;
    cfg.est.sbl.max_active = 8;

    % [新增] MP Detector Config
    cfg.det.mp.max_iter = 30;
    cfg.det.mp.damp = 0.6;
    cfg.det.mp.tol = 1e-4;
    cfg.det.mp.vmin = 1e-10;
    
    
end
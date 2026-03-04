function chan = channel_factory(cfg)
    model = cfg.channel.model;
    M = cfg.grid.M; N = cfg.grid.N; P = cfg.channel.n_paths;
    
    l_bound = min(cfg.channel.max_delay, M-1);
    k_bound = min(cfg.channel.max_doppler, floor(N/2)-1);
    
    chan.l = []; chan.k = []; chan.h = [];
    
    switch model
        case 'INT_TAPS'
            chan.l = randi([0, l_bound], P, 1); chan.l(1)=0;
            chan.k = randi([-k_bound, k_bound], P, 1);
            pwr = exp(-0.3*(0:P-1)');
            h_raw = sqrt(pwr) .* exp(1j*2*pi*rand(P,1));
            
        case 'RAYLEIGH_INT_TAPS'
            chan.l = randi([0, l_bound], P, 1); chan.l(1)=0;
            chan.k = randi([-k_bound, k_bound], P, 1);
            pwr = exp(-0.3*(0:P-1)');
            h_raw = (randn(P,1)+1j*randn(P,1)).*sqrt(pwr/2);
            
        case 'UWA_LIKE'
            % 稀疏长时延
            p_idx = sort(randperm(l_bound, P))'; p_idx(1)=0;
            chan.l = p_idx;
            % chan.k = (rand(P,1)-0.5)*2.5; % 小分数多普勒
            chan.k = randi([-k_bound, k_bound], P, 1); % 整数数多普勒
            pwr = exp(-0.15*chan.l);
            h_raw = (randn(P,1)+1j*randn(P,1)).*sqrt(pwr/2);
            
        case 'FRACTIONAL_DOPPLER'
            chan.l = randi([0, l_bound], P, 1); chan.l(1)=0;
            chan.k = randi([-k_bound, k_bound], P, 1) + (rand(P,1)-0.5);
            pwr = exp(-0.3*(0:P-1)');
            h_raw = (randn(P,1)+1j*randn(P,1)).*sqrt(pwr/2);
            
        otherwise
            error('Unknown Channel Model');
    end
    chan.h = h_raw / norm(h_raw);
    chan.n_paths = P;
end
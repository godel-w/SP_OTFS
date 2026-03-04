function cfg = cfg_case_uwa_like_M32N16()
    cfg = cfg_default();
    cfg.grid.M = 32;
    cfg.grid.N = 16;
    
    % UWA 特性: 稀疏但时延跨度大
    cfg.channel.model = 'UWA_LIKE';
    cfg.channel.n_paths = 5;
    cfg.channel.max_delay = 32; 
    cfg.channel.max_doppler = 7;
    
    % 搜索范围需覆盖
    cfg.search.l_max = 32; 
    cfg.search.k_max = 8;
    
    % CP 设置
    cfg.cp.len = 32; 
end
function cfg = cfg_case_inttaps_M16N16()
    cfg = cfg_default();
    cfg.grid.M = 16;
    cfg.grid.N = 16;
    cfg.channel.model = 'INT_TAPS';
    cfg.search.l_max = 6;
    cfg.search.k_max = 3;
end
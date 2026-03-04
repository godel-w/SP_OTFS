function [bits_hat, x_soft] = lmmse_detect(y_dd, x_p, Heff, sw2, sd2, cfg)
    % [输入增加 cfg 以获取调制阶数]
    MN = length(y_dd);
    
    % 1. 减去导频 (SIC) - 不变
    y_d = y_dd - Heff * x_p;
    
    % 2. LMMSE 均衡 - 不变
    reg = (sw2 / sd2) * eye(MN);
    x_soft = (Heff' * Heff + reg) \ (Heff' * y_d);
    
    % 3. QAM 解调 (修改部分)
    if nargout > 0
        % 归一化去功率
        x_norm = x_soft / sqrt(sd2);
        
        % 直接解调为比特流 ('OutputType', 'bit')
        % 注意：这里假设使用的是默认的 Binary 映射，如果生成时用了 Gray 需要对应
        bits_hat = qamdemod(x_norm, cfg.mod.order, ...
            'UnitAveragePower', true, ...
            'OutputType', 'bit');
    end
end
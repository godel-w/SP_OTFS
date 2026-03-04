function [Btx, Brx] = build_Btx_Brx(M, N)
    % BUILD_BTX_BRX 构建 Reduced OTFS 的基变换矩阵
    %
    % 系统定义 (Reduced OTFS):
    %   Tx: S = X * F_N'  (F_N' 是共轭转置)
    %   Rx: Y = R * F_N
    %
    % 数学恒等式: vec(AXB) = (B.' ⊗ A) vec(X)
    %
    % 推导 Btx (Tx):
    %   S = I_M * X * F_N'
    %   s = vec(S) = ((F_N').' ⊗ I_M) * vec(X)
    %              = (conj(F_N) ⊗ I_M) * vec(X)
    %
    % 推导 Brx (Rx):
    %   Y = I_M * R * F_N
    %   y_dd = vec(Y) = (F_N.' ⊗ I_M) * vec(R)
    %   注意: 这里必须用 F_N.' (非共轭转置)，虽然 DFT 矩阵对称，但在复数域下
    %   保持公式严谨性非常重要。

    F_N = dftmtx(N) / sqrt(N);
    I_M = eye(M);

    % Btx: DD -> Time
    Btx = kron(conj(F_N), I_M);

    % Brx: Time -> DD
    Brx = kron(F_N.', I_M);
end
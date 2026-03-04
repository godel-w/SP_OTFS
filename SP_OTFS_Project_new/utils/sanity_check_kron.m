function pass = sanity_check_kron()
    M = 4; N = 4;
    [Btx, Brx] = build_Btx_Brx(M, N);
    
    % Tx Test: vec(X * F') ?= Btx * vec(X)
    X = randn(M, N) + 1j*randn(M,N);
    F_H = dftmtx(N)'/sqrt(N);
    S = X * F_H;
    err1 = norm(vec(S) - Btx * vec(X));
    
    % Rx Test: vec(R * F) ?= Brx * vec(R)
    R = randn(M, N) + 1j*randn(M,N);
    F = dftmtx(N)/sqrt(N);
    Y = R * F;
    err2 = norm(vec(Y) - Brx * vec(R));
    
    pass = (err1 < 1e-10) && (err2 < 1e-10);
    if pass, fprintf('[PASS] Kronecker Check\n');
    else, fprintf('[FAIL] Kronecker Check: Tx=%.2e, Rx=%.2e\n', err1, err2); end
end
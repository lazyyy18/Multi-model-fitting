function [optimizedF] = OptimizeFundamentalLM(src, dst, initialF, weights)
    % src, dst: Nx2 對應點集
    % initialF: 初始 Fundamental Matrix (3x3)
    % weights: 每點對應的權重 (Nx1)

    f_vec = initialF(:);  % 展平為 9x1 向量

    % 設定優化目標函數
    objectiveFunction = @(f_vec) FundamentalEstimator(f_vec, src, dst, weights);

    options = optimset('Display', 'none', ...
                       'TolX', 1e-8, ...
                       'TolFun', 1e-8, ...
                       'MaxIter', 200);

    % 執行優化
    optimized_f_vec = fminunc(objectiveFunction, f_vec, options);

    % 還原為 3x3 矩陣並強制 rank-2
    F = reshape(optimized_f_vec, 3, 3);
    [U, S, V] = svd(F);
    S(3,3) = 0;
    optimizedF = U * S * V';
end

function cost = FundamentalEstimator(f_vec, src, dst, weights)
    % 將展平的向量還原為 3x3 F matrix
    F = reshape(f_vec, 3, 3);

    % % 轉為齊次座標
    % src_h = [src, ones(size(src,1), 1)]';  % 3xN
    % dst_h = [dst, ones(size(dst,1), 1)]';  % 3xN
    
    % 計算 Sampson 距離 (近似重投影誤差)
    Fx1 = F * src';       % 3xN
    Ftx2 = F' * dst';     % 3xN
    x2tFx1 = sum(dst' .* (F * src'), 1);  % 1xN

    % 分母項
    denom = Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2;

    % Sampson error（平方）
    sampson_error = (x2tFx1.^2) ./ denom;

    % 加權 cost
    cost = sum(weights(:) .* sampson_error(:));
end
function [optimizedH] = OptimizeHomographyLM(src, dst, initialH, weights)
    % src: 原始點集 (Nx2 矩陣)
    % dst: 目標點集 (Nx2 矩陣)
    % initialH: 初始的同源矩陣 (3x3 矩陣)
    % weights: 權重向量 (Nx1)

    % 將 initialH 展平為 9x1 向量
    h_vec = initialH(:);

    % fminunc 設定（現代參數名稱）
    opts = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...   % 無梯度 → 用 quasi-newton
        'Display','none', ...
        'MaxIterations',10000, ...
        'FunctionTolerance',1e-5, ...
        'StepTolerance',1e-3);

    % 最小化加權重投影誤差（Estimator 內部用 projective2d/transformPointsForward）
    h_opt = fminunc(@(h) HomographyEstimator(h, src, dst, weights), h_vec, opts);

    % 還原為 3×3 並規模正規化
    H = reshape(h_opt, 3, 3);
    optimizedH = H / H(3,3);
end

function cost = HomographyEstimator(h_vec, src, dst, weights)
    % HomographyEstimator: 計算同源矩陣的加權平方誤差
    % h_vec: 展平的 9x1 同源矩陣參數向量
    % src: 原始點集 (Nx2 矩陣)
    % dst: 目標點集 (Nx2 矩陣)
    % weights: 權重向量 (Nx1)
    %src = src(:,1:3);
    %dst = dst(:,1:3);
    % 將展平的 h_vec 轉換為 3x3 矩陣
    % 將 9x1 參數還原為 3x3 H
    H = reshape(h_vec, 3, 3);

    projected_points =  (H * src')';
    projected_points(projected_points(:,3)<1e-8,3) = 1e-8;
    projected_pts = projected_points(:, 1:2) ./ projected_points(:,3);    
    
    % 加權平方重投影誤差
    e = dst(:,1:2) - projected_pts(:,1:2);                      % N×2
    cost = sum(weights(:) .* sum(e.^2, 2));                 % 1×1
end
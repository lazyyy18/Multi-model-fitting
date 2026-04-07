function residuals = sampson_residual(f, pts1, pts2)
    F = reshape(f, 3, 3);

    N = size(pts1, 1);
    pts1_h = pts1';  % 3×N
    pts2_h = pts2';  % 3×N

    Fx1 = F * pts1_h;             % 3×N
    Ftx2 = F' * pts2_h;           % 3×N

    x2tFx1 = sum(pts2_h .* (F * pts1_h), 1);  % 1×N

    % 分母為 gradient 向量的平方和
    denom = Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2;

    residuals = (x2tFx1.^2)  ./ denom;  % → N×1, always ≥ 0
    residuals = sqrt(residuals(:));  % N×1
end
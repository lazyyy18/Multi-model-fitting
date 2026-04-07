function [error] = res_error(X,Y,model)
    % 使用单应性矩阵 H 对所有点进行投影变换
    projected_pts_homog = (model * X')';    
    % 归一化投影点坐标
    projected_pts = projected_pts_homog(:, 1:2) ./ projected_pts_homog(:,3);    
    % 计算每个点的投影误差
    error = vecnorm(projected_pts - Y(:,1:2), 2, 2);
end
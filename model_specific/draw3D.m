function [counts, length_bin, angle_bin, xbin, ybin] = draw3D(source_pt, target_pt)
    
    % 计算 motion 向量
    motion = target_pt(:,1:2) - source_pt(:,1:2);
    
    % 计算 motion 向量的长度
    motion_lengths = hypot(motion(:,1), motion(:,2));
    % 计算 motion 的角度
    motion_angle = atan2d(motion(:, 2), motion(:, 1));

    % 將 labels_gmm 對齊 gth
    N = length(motion_lengths);
    sigma = std(motion_angle, 1);
    bin_width_fd_angle = 3.49 *sigma /(N^(1/3));
    % 找到最小與最大值
    min_input = round(min(motion_angle) * 10) / 10;
    max_input = round(max(motion_angle) * 10) / 10;
    
    % 直接生成範圍
    edges_angle = min_input - bin_width_fd_angle / 2 : bin_width_fd_angle : max_input + bin_width_fd_angle / 2;
    if(edges_angle(end) <= max_input)
        edges_angle(end) = max_input+1;
    end
    sigma = std(motion_lengths, 1);  % 計算標準差 (σ)
    bin_width_fd_length = 3.49 * sigma /(N^(1/3));


    
    % 找到最小與最大值
    min_input = round(min(motion_lengths) * 10) / 10;
    max_input = round(max(motion_lengths) * 10) / 10;
    
    % 直接生成範圍
    edges_length = min_input - bin_width_fd_length / 2 : bin_width_fd_length : max_input + bin_width_fd_length / 2;
    if(edges_length(end) <= max_input)
        edges_length(end) = max_input+1;
    end

    [counts, length_bin, angle_bin, xbin, ybin] = histcounts2(motion_angle,motion_lengths, edges_angle, edges_length);
end


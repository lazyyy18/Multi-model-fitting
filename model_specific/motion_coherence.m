function [cluster_pts] = motion_coherence(source_pt, target_pt,numModels, peak_mun, gth)
    % 计算 motion 向量
    motion = target_pt(:,1:2) - source_pt(:,1:2);

    % 更快更穩：長度用 hypot、角度用 atan2（度或弧度都可）
    motion_lengths = hypot(motion(:,1), motion(:,2));
    motion_angles  = atan2d(motion(:,2), motion(:,1));  % 想再快可用 atan2 改成弧度    
    % [cluster_pts] = motion_histogram_2d(motion_lengths, motion_angles,numModels, peak_mun, 0);
    [cluster_pts] = motion_histogram_fast15(motion_lengths, motion_angles,peak_mun);
    % cluster_pts = motion_histogram_2d_m(motion_lengths, motion_angles, numModels, peak_mun, 0);
end


% function cluster_pts = motion_histogram_2d_m(lengths, angles, numModels, peak_mun, witch_first)
%     cluster_pts = motion_histogram_2d_mex(lengths, angles, numModels, peak_mun, witch_first);
%     % N = numel(angles);
%     % 
%     % %% =======================
%     % %  ANGLE HISTOGRAM (FASTER)
%     % %% =======================
%     % 
%     % minA = round(min(angles) * 10) / 10;
%     % maxA = round(max(angles) * 10) / 10;
%     % 
%     % meanA = mean(angles);
%     % sigmaA = sqrt(mean(angles.^2) - meanA^2);
%     % bwA = 3.49 * sigmaA / (N^(1/3));
%     % 
%     % edgesA = minA - bwA/2 : bwA : maxA + bwA/2;
%     % if edgesA(end) < maxA
%     %     edgesA(end) = maxA;
%     % end
%     % 
%     % % -------- O(N) quantized binning (替代 histcounts) --------
%     % inv_bwA = 1 / bwA;
%     % binA = floor((angles - edgesA(1)) * inv_bwA) + 1;
%     % binA(binA < 1) = 1;
%     % binA(binA >= numel(edgesA)) = numel(edgesA)-1;
%     % 
%     % anglecount = accumarray(binA, 1, [numel(edgesA)-1, 1])';
%     % 
%     % [~, peak1] = max(anglecount);
%     % leftA  = max(1, peak1 - peak_mun);
%     % rightA = min(numel(anglecount), peak1 + peak_mun);
%     % 
%     % angle_low  = edgesA(leftA);
%     % angle_high = edgesA(rightA+1);
%     % 
%     % compact_idx = (angles >= angle_low & angles < angle_high);
%     % 
%     % 
%     % %% ========================
%     % %  LENGTH HISTOGRAM (FASTER)
%     % %% ========================
%     % 
%     % second_info = lengths;
%     % lengths_f = lengths(compact_idx);
%     % 
%     % if isempty(lengths_f)
%     %     cluster_pts = [];
%     %     return;
%     % end
%     % 
%     % minL = round(min(lengths_f) * 10) / 10;
%     % maxL = round(max(lengths_f) * 10) / 10;
%     % 
%     % meanL = mean(lengths_f);
%     % sigmaL = sqrt(mean(lengths_f.^2) - meanL^2);
%     % bwL = 3.49 * sigmaL / (N^(1/3));
%     % 
%     % edgesL = minL - bwL/2 : bwL : maxL + bwL/2;
%     % if edgesL(end) < maxL
%     %     edgesL(end) = maxL;
%     % end
%     % 
%     % % -------- O(N) quantized binning (替代 histcounts) --------
%     % inv_bwL = 1 / bwL;
%     % binL_full = floor((lengths_f - edgesL(1)) * inv_bwL) + 1;
%     % binL_full(binL_full < 1) = 1;
%     % binL_full(binL_full >= numel(edgesL)) = numel(edgesL)-1;
%     % 
%     % lengthcount = accumarray(binL_full, 1, [numel(edgesL)-1, 1])';
%     % 
%     % [~, peak2] = max(lengthcount);
%     % 
%     % leftL  = max(1, peak2 - peak_mun);
%     % rightL = min(numel(lengthcount), peak2 + peak_mun);
%     % 
%     % len_low  = edgesL(leftL);
%     % len_high = edgesL(rightL+1);
%     % 
%     % % 完全等價 replicate
%     % if peak2 ~= numel(lengthcount)
%     %     validL = (second_info >= len_low & second_info <  len_high);
%     % else
%     %     validL = (second_info >= len_low & second_info <= len_high);
%     % end
%     % 
%     % cluster_pts = find(compact_idx & validL);
% end



%F:使用% 使用 scott's 方法确定 bin 数

function [cluster_pts]  = motion_histogram_2d(motion_lengths, motion_angle,numModels, peak_mun,witch_first)
% % 
    first = motion_angle;
    second = motion_lengths;
    % first_N = length(first);
    % first = motion_lengths;
    % second = motion_angle;
    first_N = length(second);
    compact_idx = [];

    %angle
    bins_angle = 0; bins_length = 0;

    % 使用 Freedman-Diaconis 方法确定 bin 数
    % bin_width_fd_angle = 2 * iqr(motion_angle) / (first_N^(1/3));
    % 使用 scott's 方法确定 bin 数
    mean_val = mean(first);  % 計算均值
    sigma = sqrt(mean(first.^2) - mean_val^2);  % 計算標準差 (σ)
    bin_width_fd_angle = 3.49 *sigma /(first_N^(1/3));
    %使用(Sturges' Rule)
    %in_width_fd_angle = log2(first_N) + 1;

    % bins_angle = round((max(motion_angle) - min(motion_angle)) * 10.0) / 10.0 / bin_width_fd_angle;
    % 找到最小與最大值
    min_input = round(min(first) * 10) / 10;
    max_input = round(max(first) * 10) / 10;

    % 直接生成範圍
   % edges_angle = [min_input - bin_width_fd_angle / 2, min_input : bin_width_fd_angle : (max_input + bin_width_fd_angle)];
    edges_angle = min_input - bin_width_fd_angle / 2 : bin_width_fd_angle : max_input + bin_width_fd_angle / 2;
    if(edges_angle(end) < max_input)
        edges_angle(end) = max_input;
    end
    % 创建 bin
    % edges_angle = linspace(min(motion_angle), max(motion_angle),bins_angle+ 1);
%     [anglecount,edges_angle] = histcounts(motion_angle);
    anglecount = histcounts(first, edges_angle);
    [~, peak1] = max(anglecount);
    %取angle的最大角度
    corres_angle_bin = discretize(first, edges_angle);
%     angle_bin(isnan(length_bin) & second_info < edges_length(1)) = 1; % 小於最小邊界
%     angle_bin(isnan(length_bin) & second_info > edges_length(-1)) = length(edges_length); % 大於最大邊界

    leftpeak = max(1, peak1 - peak_mun);       % 防止索引低於 1
    rightpeak = min(length(anglecount), peak1 + peak_mun); % 防止索引超過範圍

    compact_idx = ismember(corres_angle_bin, leftpeak:rightpeak); 
    % % 畫圖
    % figure
    % bar(anglecount)
    % hold on;
    % set(gca, 'XTickLabel', round(edges_angle(1:end-1), 2));
    % xlabel('Motion angle');
    % ylabel('frequency');
    % title('2D Histogram of Motion angle');
    % hold off
    % % length


    second_info = second;

    second = second(compact_idx);
    %iqr_value = iqr(motion_lengths);
    %motion_lengths = motion_lengths(:);
    %iqr_length = iqr(motion_lengths);
    % 使用 Freedman-Diaconis 方法确定 bin 数
    % bin_width_fd_length = 2 * iqr(motion_lengths) / (first_N ^(1/3));
    % 使用 scott's 方法确定 bin 数
    mean_val = mean(second);  % 計算均值
    sigma = sqrt(mean(second.^2) - mean_val^2);  % 計算標準差 (σ)
    bin_width_fd_length = 3.49 * sigma /(first_N^(1/3));


    % 使用 (Sturges' Rule) 方法确定 bin 数
    %bin_width_fd_length = log2(second_N) + 1;
    %自訂數量
    %bin_width_fd_length = 70;
    % bins_length = ceil((max(motion_lengths) - min(motion_lengths)) / bin_width_fd_length);
    % 创建 bin
    % edges_length = linspace(min(motion_lengths), max(motion_lengths), bins_length + 1);
%     [lengthcount,edges_length] = histcounts(motion_lengths);

    % 找到最小與最大值
    min_input = round(min(second) * 10) / 10;
    max_input = round(max(second) * 10) / 10;
    % 直接生成範圍
 %   edges_length = [min_input - bin_width_fd_length / 2, min_input : bin_width_fd_length : (max_input + bin_width_fd_length)];
    edges_length = min_input - bin_width_fd_length / 2 : bin_width_fd_length : max_input + bin_width_fd_length / 2;
    if(edges_length(end) < max_input)
        edges_length(end) = max_input;
    end
    [lengthcount] = histcounts(second, edges_length);
    corres_length_bin = discretize(second_info, edges_length);
    corres_length_bin(isnan(corres_length_bin) & second_info < edges_length(1)) = 1; % 小於最小邊界
    corres_length_bin(isnan(corres_length_bin) & second_info > edges_length(end)) = length(edges_length); % 大於最大邊界

    % figure 
    % bar(lengthcount)
    % hold on;
    % set(gca, 'XTickLabel', round(edges_length(1:end-1), 2));
    % xlabel('Motion length');
    % ylabel('frequency');
    % title('2D Histogram of Motion Lengths');
    % hold off

    [~,peak2] = max(lengthcount);
    leftpeak = max(1, peak2 - peak_mun);       % 防止索引低於 1
    rightpeak = min(length(lengthcount), peak2 + peak_mun); % 防止索引超過範圍
    width = [edges_length(leftpeak) ; edges_length(rightpeak+1)];
    if peak2 ~= length(lengthcount)
        cluster_pts = find(compact_idx == 1 & width(1) <= second_info & second_info < width(2));
    else
        cluster_pts = find(compact_idx == 1 & width(1) <= second_info & second_info <= width(2));
    end
end


function cluster_pts = motion_histogram_fast15(lengths, angles, peak_mun)

    N = numel(angles);

    %% =======================
    %  ANGLE HISTOGRAM (FASTER)
    %% =======================

    minA = round(min(angles) * 10) / 10;
    maxA = round(max(angles) * 10) / 10;

    % meanA = mean(angles);
    % sigmaA = sqrt(mean(angles.^2) - meanA^2);
    sigmaA = std(angles,1);
    bwA = 3.49 * sigmaA / (N^(1/3));
    % bwA = 2 * iqr(angles) / (N^(1/3));

    edgesA = minA - bwA/2 : bwA : maxA + bwA/2;
    if edgesA(end) < maxA
        edgesA(end) = maxA;
    end


    % -------- O(N) quantized binning (替代 histcounts) --------
    inv_bwA = 1 / bwA;
    binA = floor((angles - edgesA(1)) * inv_bwA) + 1;
    binA(binA < 1) = 1;
    binA(binA >= numel(edgesA)) = numel(edgesA)-1;

    anglecount = accumarray(binA, 1, [numel(edgesA)-1, 1])';

    [~, peak1] = max(anglecount);
    leftA  = max(1, peak1 - peak_mun);
    rightA = min(numel(anglecount), peak1 + peak_mun);

    angle_low  = edgesA(leftA);
    angle_high = edgesA(rightA+1);

    compact_idx = find(angles >= angle_low & angles < angle_high);

    %---------plot angle bin ------------
    binCenterA = (edgesA(1:end-1) + edgesA(2:end)) / 2;
    figure; hold on;
    
    histogram(angles, ...
        'BinEdges', edgesA, ...
        'FaceColor', '#C0C0C0', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1);
    
    % 高亮區間
    % yl = ylim;
    histogram(angles(compact_idx), ...
        'BinEdges', edgesA(leftA:rightA+1), ...
        'FaceColor', '#6495ED', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1);
    xlabel('Angle');
    ylabel('Frequency');
    title('Angle Histogram');


    %% ========================
    %  LENGTH HISTOGRAM (FASTER)
    %% ========================

    second_info = lengths;
    lengths_f = lengths(compact_idx);
    if isempty(lengths_f)
        cluster_pts = [];
        return;
    end

    minL = round(min(lengths_f) * 10) / 10;
    maxL = round(max(lengths_f) * 10) / 10;

    % meanL = mean(lengths_f);
    sigmaL = std(lengths_f,1);
    bwL = 3.49 * sigmaL / (N^(1/3));
    % bwL = 2 * iqr(lengths_f) / (N^(1/3));

    edgesL = minL - bwL/2 : bwL : maxL + bwL/2;
    if edgesL(end) < maxL
        edgesL(end) = maxL;
    end

    % -------- O(N) quantized binning (替代 histcounts) --------
    inv_bwL = 1 / bwL;
    binL_full = floor((lengths_f - edgesL(1)) * inv_bwL) + 1;
    binL_full(binL_full < 1) = 1;
    binL_full(binL_full >= numel(edgesL)) = numel(edgesL)-1;

    lengthcount = accumarray(binL_full, 1, [numel(edgesL)-1, 1])';

    [~, peak2] = max(lengthcount);

    leftL  = max(1, peak2 - peak_mun);
    rightL = min(numel(lengthcount), peak2 + peak_mun);

    len_low  = edgesL(leftL);
    len_high = edgesL(rightL+1);

    % 完全等價 replicate
    if peak2 ~= numel(lengthcount)
        validL = (lengths_f >= len_low & lengths_f <  len_high);
    else
        validL = (lengths_f >= len_low & lengths_f <= len_high);
    end
    

    %---------plot length bin----------
    figure; hold on;
    
    % ---- 全部先畫成黑色（不要的）----
    histogram(lengths_f, ...
        'BinEdges', edgesL, ...
        'FaceColor', '#C0C0C0', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1);
    
    % ---- 畫選中的 bins（藍色）----
    histogram(lengths_f(validL), ...
        'BinEdges', edgesL(leftL:rightL+1), ...
        'FaceColor', '#6495ED', ...
        'EdgeColor', 'none',...
        'FaceAlpha', 1);
    
    xlabel('Length');
    ylabel('Frequency');
    title('Length Histogram (After Angle Filtering)');
    % grid on;


    cluster_pts = compact_idx(validL);
end


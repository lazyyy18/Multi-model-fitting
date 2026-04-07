function [time,label,H] = GH_motion_coosac(X, Y, numModels, img1, img2)
    rng(42);
    threshold = 0.1;
    [nor_x, ~] = normalise2dpts(X');
    [nor_y, ~] = normalise2dpts(Y');
    nor_x = nor_x';
    nor_y = nor_y';
    tic;
    imgy = Y + [size(img1,2) , 0 , 0]; 
    N = size(X,1);
    [hiscount, ~, ~, xbin, ybin] = draw3D(X, imgy); 
    inlier_rate = sum(hiscount(hiscount > 30)) / N;
    if (inlier_rate) > 0.75
        R = 10;
        threshold = threshold * 1.0;
    elseif (inlier_rate) > 0.5
        R = 9;
        threshold = threshold * 1.83;
    else
        R = 8;
        threshold = threshold * 2.45;
    end
    H = cell(1,numModels); 
    model_threshold = zeros(numModels,1);
    model_residual = ones(N,numModels) * 1e7; 
    neighbor_inlier = false(N,numModels); 

    best_inlier = false(N,1); 
    active = true(N,1);
    point_idx = 1:N;
    motion_field = imgy - X;
    motion_field = motion_field(:,1:2);
    pca_X = [nor_x(:, 1:2), nor_y(:, 1:2), motion_field];

    for i = 1:numModels
        
        [~, idx] = max(hiscount, [], 'all');
        [row, col] = ind2sub(size(hiscount), idx);
        cluster_pts = find(xbin == row & ybin == col);
        
        if (length(cluster_pts) < 10)
            continue
        end
        before = rng;
        BIC = zeros(1, numModels - i + 1);
        GMModels = cell(1, numModels - i + 1);
        stream = RandStream('mt19937ar','Seed',22);
        options = statset('MaxIter',300, 'UseParallel',false, ...
                          'UseSubstreams',false, 'Streams',stream);
        X_cluster = pca_X(point_idx(cluster_pts), :);
        parfor gmm_model_idx = 1: numModels - i + 1
            GMModels{gmm_model_idx} = fitgmdist(X_cluster, gmm_model_idx, 'RegularizationValue', 1e-5, ...
                'Options', options, 'Replicates', 5);
            BIC(gmm_model_idx) = GMModels{gmm_model_idx}.BIC;
        end
        [~, bestK] = min(BIC);
        bestGMM = GMModels{bestK};
        GM_idx = cluster(bestGMM, X_cluster);
        major_label = mode(GM_idx);
        GM_pts = find(GM_idx == major_label);
        if length(GM_pts) >= 10
            cluster_pts = cluster_pts(GM_pts);
        end
        rng(before);
        [H{i}, model_threshold(i), ~] = coosac_motion(nor_x(point_idx,:), nor_y(point_idx,:), cluster_pts, threshold, R); 
        model_residual(:,i) = res_error(nor_x, nor_y, H{i});
        neighbor_inlier(:,i) = (model_residual(:,i) < model_threshold(i));
        active = active & ~neighbor_inlier(:,i);
        point_idx = find(active);
        [hiscount, ~, ~, xbin, ybin] = draw3D(X(point_idx,:), imgy(point_idx,:));
    end
    eigen_inlier = false(N,numModels);
    neighbor_count = sum(neighbor_inlier, 2);
    onesN = ones(N,1);
    
    parfor i = 1:numModels
        a = neighbor_inlier(:,i) & (neighbor_count == 1);
        initial_H = model(nor_x(a,:), nor_y(a,:), onesN(a));
        optimizedH = OptimizeHomographyLM( ...
        nor_x(a,:), nor_y(a,:), initial_H, onesN(a));
        model_residual(:,i) = res_error(nor_x, nor_y, optimizedH);
        eigen_inlier(:,i) = model_residual(:,i) < model_threshold(i);
    end
    
    [~, best_idx] = min(model_residual, [], 2);
    label = zeros(N,1);
    
    for i = 1:numModels
        a = eigen_inlier(:,i) & (best_idx == i);
        label(a) = i;
    end
    
    time = toc;
end


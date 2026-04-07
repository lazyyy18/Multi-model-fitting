function [best_model, iteration,best_inliers,jump] = mix_plane_tinycoosac(x,inlier_threshold,weight,cel)
    best_j = 0;
    jump = false;
    p = 6;
    iteration = 0;
    best_inliers  = [];
    best_model = eye(3);
    tinyconfidence = 1000;
    Xs = x(:,1:3);
    Ys = x(:,4:6);
    Nx = size(x,1);
    Nx_p = Nx;
    best_inlier_cnt = 0;
    [psub, spidx] = datasample(x(cel,:), p, Replace= false);
     adjusted_weights = adjustWeights(weight(cel(spidx)));       
    while iteration <= 40 && tinyconfidence > p+1 && jump == 0
        H = model(psub(:,1:3), psub(:,4:6), adjusted_weights);
        error = res_error(Xs, Ys, H);
        % 只算數量，不建 logical 向量（省很多）
        inlier_cnt = sum(error < inlier_threshold);
    
        t = 1 - error / inlier_threshold;
        model_quality = inlier_cnt + sum(t(t>0));
    
        if model_quality > best_j
            best_model = H;
            best_j = model_quality;
            best_inlier_cnt = inlier_cnt;   % <- 存數量就好
        end
    
        tinyconfidence = cal_confidence(Nx_p, best_inlier_cnt, p, iteration, 0.995);
        [~, spidx] = mink(error, 15);
        psub = x(spidx(15-p+1:15),:);
        adjusted_weights = adjustWeights(weight(spidx));
        jump = norm(error(spidx)) < 1e-4;  
        iteration = iteration + 1;
    end
end


function [confidence] = cal_confidence(s,hd,p,k,max)
    a = (1-max)^(1/k);
    b = (1-a)^(1/p);
    confidence = (s-hd) *b;
end
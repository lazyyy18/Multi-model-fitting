function [best_model, iteration,best_inliers,jump] = mix_motion_tinycoosac(x,inlier_threshold,weight,c_rel)
    best_j = 0;
    jump = false;
    p =10;
    tinyconfidence = 100;
    iteration = 0;
    best_inliers  = [];
    best_model = eye(3); 
    Nx_p = size(x,1);
    [psub,pidx] = datasample(x(c_rel,:), p, Replace= false);
    adjusted_weights = adjustWeights(weight(c_rel(pidx)));
    while  (iteration <= 40 && tinyconfidence > p+1 && jump == 0 )
        H = model(psub(:,1:3), psub(:,4:6), adjusted_weights);
        
        error = res_error(x(:,1:3), x(:,4:6), H);
   
        % === 只取最小 15 個 residual ===
        [~, idx15] = mink(error, 15);
        
        psub = x(idx15(15-p+1:15), :);
        adjusted_weights = adjustWeights(weight(idx15));
        inliers = error <  inlier_threshold;
        num_inliers = sum(inliers);
        t = 1 - error / inlier_threshold;
        model_quality = num_inliers + sum(t(t>0));
        if model_quality > best_j
            best_model = H;
            best_j = model_quality;
        end      
        tinyconfidence = cal_confidence(Nx_p, num_inliers, p, iteration, 0.995);
        jump = norm(error(idx15)) < 1e-3;       
        iteration = iteration + 1;
    end

end


function [confidence] = cal_confidence(s,hd,p,k,max)
    a = (1-max)^(1/k);
    b = (1-a)^(1/p);
    confidence = (s-hd) *b;
end
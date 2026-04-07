function [H,best_inlie_res, c_rel] = coosac_motion(norsrc,nortar,c_rel,threshold, R)

   %init data
    N = length(norsrc);
    best_j = 0;
    best_inlie_res = 1e7;
    p = 10;
    maxIters = 200000;
    final_iteration = 0;
	QR_times = 1;
    max_confidence =  0.995;
    mask = false(N,1);
    mask(c_rel) = 1; 
    %%%
    best_r = 0;
    history_inlier = zeros(N,1);
    H = eye(3);
    weight = zeros(N,1); 
    while (final_iteration < maxIters && QR_times <= R && length(c_rel)>=10)
        x = [norsrc, nortar];
        [best_model, iteration, ~, jump] = mix_motion_tinycoosac(x,threshold,weight,c_rel);
        if(jump)
    		final_iteration = final_iteration *2 + iteration;
        end
        final_iteration = final_iteration + iteration;
        glo_error = res_error(norsrc, nortar, best_model);
        Cinl = glo_error <  threshold;
        model_quality = sum(Cinl) + sum(max(0, 1 - glo_error / threshold));
        part = 7;
        max_sigma = 0;
        max_sigma = max(max_sigma,max(glo_error(Cinl)));
        prob = computeInlierProbability(glo_error, threshold);
        weight = weight + max(min(part, part - floor(glo_error ./ (threshold / 3))), 0);
        for i = 1:N
            if Cinl(i)
                if ~mask(i) && ( history_inlier(i) || (~history_inlier(i) && rand() <= prob(i)))
                    mask(i) = true;
                end
            else
                if history_inlier(i) && rand() > prob(i)
                    mask(i) = false;
                end
            end
        end
        history_inlier = Cinl;
        c_rel = find(mask);
        if ~isempty(glo_error(Cinl))
            max_sigma = max(max_sigma, max(glo_error(Cinl)));
            threshold = (max_sigma * 0.5 + threshold * 0.5);
        end
        if model_quality > best_j
            H = best_model;
            best_inlie_res = threshold;
            best_j = model_quality;
        end
        final_confidence = cal_confidence(N,sum(Cinl),p,final_iteration+1,max_confidence);
		final_iteration = final_iteration +1;
		QR_times = QR_times +1;
        if (final_confidence <= p+1)
            break;
        end
    end
    initial_H = model(norsrc(c_rel,:),nortar(c_rel,:),weight(c_rel));
    eigen_H = OptimizeHomographyLM(norsrc(c_rel,:), nortar(c_rel,:), initial_H, weight(c_rel));
    eigen_error = res_error(norsrc,nortar,eigen_H);
    eigen_inlier = eigen_error(:) <  threshold;
    eigen_score = sum(eigen_inlier) + sum(eigen_error <=  threshold/2);
    if eigen_score > best_j 
		H = eigen_H;
        best_inlie_res = best_inlie_res;
    end
end

function [confidence] = cal_confidence(s,hd,p,k,max)
    a = (1-max)^(1/k);
    b = (1-a)^(1/p);
    confidence = (s-hd) *b;
end


function [probability] = computeInlierProbability(r, sigma)
    a = 2 / 2.0; 
    z = (r.^2) / (2 * (sigma.^2)); % 計算 z
    % gammainc(z, a) 返回的是 P(a, z) = gamma(a, z) / Gamma(a)
    cdf = gammainc(z, a, 'lower'); % 計算歸一化的下不完全伽馬函數
    probability = 1.0 - cdf; % 計算最終的概率
end
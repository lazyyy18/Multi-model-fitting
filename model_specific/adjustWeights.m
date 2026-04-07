function [adjusted_weights]  = adjustWeights(original_weights)
    adjusted_weights = original_weights;
    if sum(original_weights) == 0
			adjusted_weights(:) = 1;
    else  
			adjusted_weights(adjusted_weights == 0) =  1e-6;		
    end
end
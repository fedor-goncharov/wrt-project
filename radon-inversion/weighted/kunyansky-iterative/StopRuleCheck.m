% Stopping a cycle of iterative inversions in 'iterative_wradon_inversion.m'
% Stopping rule is based on stopping_params in the input arguments of 'iterative_wradon_inversion.m'

function stop = StopRuleCheck(stopping_params, cycle_params)
  stop = false; 
  
  stopping_flag = stopping_params{1};
  stopping_iterations = stopping_params{2}; 
  stopping_precision = stopping_params{3};
  
  iteration = cycle_params{1};
  current_precision = cycle_params{2};
  
  if (stopping_flag == 'FLAG_NITERATIONS') 
    if (iteration >= stopping_iterations) 
      stop = true;
    end
  end
  
  if (stopping_flag == 'FLAG_PRECISION')
    if (current_precision < stopping_precision)
      stop = true;
    end 
  end 
  
  if (stopping_flag == 'FLAG_COMPOSITE')
    if (iteration >= stopping_iteration || current_precision < stopping_precision)
      stop = true;
    end
  end 
end
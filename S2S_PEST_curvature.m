if PEST_ON==1
    stepsize_min          = 0.02;
    curvature_limits      = [-1 curvature_fixed]; % for concave
    %curvature_limits      = [curvature_fixed 1]; % for convex
    
    if exist('dyn','var') && dyn.trialNumber > 5 && sum([trial.completed])>=5
        %selected_target_t=trial(dyn.trialNumber-1).target_selected(PEST_effector);
        completed_trials= trial([trial.completed]);
        success_rate= sum(completed_trials(end-4:end).success)/5;
        
        preferred_behaviour=success_rate==desired_success_rate;
        if trial(dyn.trialNumber-1).completed
            [curvature_adjusted, stepsize] = PEST(curvature_adjusted,curvature_limits,trial(dyn.trialNumber-1).success,preferred_behaviour,stepsize,stepsize_min);
        end
    else
        stepsize              = 1;
        curvature_adjusted    = curvature_start;
    end
    PEST_list(dyn.trialNumber)=fix_eye_x;
             %stepsize
             %PEST_list
    if stepsize<=stepsize_min
            dyn.state = STATE.CLOSE; 
            PEST_list
            return
    end
else
    curvature_adjusted             = curvature_start;
end


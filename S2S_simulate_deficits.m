function out=S2S_simulate_deficits(impairment,sampling_rate)
% generation of data needed to run the S2S_session_plot
% -> production of fake data as from the human match-to-sample task
% out needs to be generated
% it is a structure consisting of several stuctures:
% keys, selected, task, timing, states, saccades, reaches, binary, counts,
% statistic, correlation, rundescriptions, emptyflag

% how to:
% 1) data_generation
% 2) M2S_output=S2S_session_analysis(out,0,6);
% 3) S2S_session_plot(M2S_output)


% for the purpose of data-generation, only the following stuctures need to
% be generated:
% binary & saccades


%%

%impairment='sensory_neglect_L';
%impairment='motivational_neglect_L';
%impairment='object_related_neglect_L';
%impairment='ideal';

clear saccades;
saccade_duration=1.5;
inspection_time_distribution_factor=0.1;

% sequence generation
sphi1=sin(pi/3)*20;
sphi2=sin(pi/3)*10;
cphi1=cos(pi/3)*20;
cphi2=cos(pi/3)*10;
possible_tar_pos = [ -cphi1 + sphi1*1i,  -cphi2  + sphi2*1i,  cphi2 + sphi2*1i,  cphi1 + sphi1*1i,...
    -20.0 + 0.0i,  -10.0 + 0.0i,  10.0 + 0.0i, 20.0 + 0.0i, ...
    -cphi1 - sphi1*1i,  -cphi2  - sphi2*1i, cphi2 - sphi2*1i,  cphi1 - sphi1*1i];

%     possible_tar_pos = [ -22.0 + 5.0i,  -14.0 + 5.0i,  -6.0 + 5.0i,  6.0 + 5.0i,  14.0 + 5.0i, 22.0 + 5.0i,...
%                          -22.0 + 0.0i,  -14.0 + 0.0i,  -6.0 + 0.0i,  6.0 + 0.0i,  14.0 + 0.0i, 22.0 + 0.0i, ...
%                          -22.0 - 5.0i,  -14.0 - 5.0i,  -6.0 - 5.0i,  6.0 - 5.0i,  14.0 - 5.0i, 22.0 - 5.0i];

%target distributions
all_target_distributions=possible_tar_pos([1,8,9;2,7,10;3,6,11;4,5,12]);

n_targets = 3;
%
%     possible_convex_sides = {'LR', 'LR', 'LR', 'LR'};
%     possible_convexities = [-0.5,-0.3,0.3,0.5];


possible_convex_sides = {'L', 'LR', 'R', 'LR'};
possible_convexities = [-0.7,-0.7,-0.7,-0.3];


%     possible_tar_pos = [13.0 + 5.0i, 13.0 - 5.0i, -13.0 - 5.0i, -13.0 + 5.0i];
%     n_targets = 4;

%     possible_convex_sides = {'LR', 'TB', 'L', 'T'};
%     possible_convexities = [-1,-1,1,1,];



N_repetitions=10;

possible_shapes = [1,2,3,4];
shape_to_repeat = [3,4,1,2];
match_position_conditions = [-1,1];
distractor_position_conditions = [0];
L_or_R_cue=[15,-15];

target_distribution_conditions=1:size(all_target_distributions,1);
sequence=repmat(combvec(L_or_R_cue,possible_shapes,match_position_conditions,distractor_position_conditions,target_distribution_conditions),1,N_repetitions);
N_trials=size(sequence,2);


%% impairment parameters !!!!
switch impairment
    case 'ideal'
%         target_preferences =[1,  1,  1,   1,   1,  1,...
%             1,  1,  1,   1,   1,  1,...
%             1,  1,  1,   1,   1,  1];

        target_preferences =[1,  1,  1,  1,...
                             1,  1,  1,  1,...
                             1,  1,  1,  1];
        Confusion_matrix.L =[1 ,0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0;...
            0, 0, 0, 1];
        Confusion_matrix.R =[1 ,0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0;...
            0, 0, 0, 1];
        Inspection_time_bias_L_R    =[0,0];
        RT_time_bias_L_R.L            =[0,0];
        RT_time_bias_L_R.R            =[0,0];
     case 'sensory_neglect_L'
%         target_preferences =[1,  1,  2,   2,   1,  1,...
%             1,  2,  4,   4,   2,  1,...
%             1,  4,  6,   6,   4,  1];

        target_preferences =[1,  1,  1,  1,...
            1,  2,  2,  1,...
            1,  4,  4,  1];
        Confusion_matrix.L =[0.3, 0.3,  0.3,    0.3;...
            0.3, 0.3,  0.3,    0.3;...
            0.3,   0.3,  0.3,  0.3;...
            0.3,   0.3,  0.3,  0.3];
        Confusion_matrix.R =[0.9, 0.1,  0,   0;...
            0.05, 0.9,  0.05, 0;...
            0,   0.05,  0.9, 0.05;...
            0,     0,  0.1, 0.9];
        Inspection_time_bias_L_R=[0,0];
        RT_time_bias_L_R.L            =[0.3,0.3];
        RT_time_bias_L_R.R            =[0,0];
    case 'motivational_neglect_L'
        
        target_preferences =[0.7,  1,   30,  20,...
            0.7,  1,   40,  20,...
            0.7,  1,   50,  20];
        
        Confusion_matrix.L = [0.6, 0.38,   0.015,0.005;...
            0.195, 0.6,   0.195, 0.01;...
            0.01,0.195,   0.6, 0.195;...
            0.005,0.015,  0.38, 0.6];
        Confusion_matrix.R =[0.6, 0.38,   0.015,0.005;...
            0.195, 0.6,   0.195, 0.01;...
            0.01,0.195,   0.6, 0.195;...
            0.005,0.015,  0.38, 0.6];
        Inspection_time_bias_L_R    =[0,0.3];
        RT_time_bias_L_R.L            =[0.3,0];
        RT_time_bias_L_R.R            =[0.3,0];
        
    case 'object_related_neglect_L'
%         target_preferences =[1,  1,  2,   2,   1,  1,...
%             1,  2,  4,   4,   2,  1,...
%             1,  4,  6,   6,   4,  1];
        target_preferences =[1,  1,  1,  1,...
            1,  2,  2,  1,...
            1,  4,  4,  1];
        Confusion_matrix.L =[ 0.48, 0.02,  0.02, 0.48;...
            0.02, 0.48,  0.48, 0.02;...
            0.02, 0.48,  0.48, 0.02;...
            0.48, 0.02,  0.02, 0.48];
        Confusion_matrix.R =[ 0.48, 0.02,  0.02, 0.48;...
            0.02, 0.48,  0.48, 0.02;...
            0.02, 0.48,  0.48, 0.02;...
            0.48, 0.02,  0.02, 0.48];
        Inspection_time_bias_L_R    =[0,0];
        RT_time_bias_L_R.L            =[0,0];
        RT_time_bias_L_R.R            =[0,0];
end


CM_fieldnames=fieldnames(Confusion_matrix);
for FN_idx=1:numel(CM_fieldnames)
    FN=CM_fieldnames{FN_idx};
    Confusion_RT_matrix.(FN)=Confusion_matrix.(FN);
    Confusion_RT_matrix.(FN)=Confusion_RT_matrix.(FN).*(1-repmat(Confusion_matrix.(FN)(eye(size(Confusion_matrix.(FN),1))==1),1,size(Confusion_matrix.(FN),2)))*1;
end;

%run task generation
for trial=1:N_trials
    task(trial).type     = 6;                                              % masked
    task(trial).effector = 0;                                              % eye
    binary(trial).choice = 0;                                              % eye
    states(trial).state_abo=-1;
    states(trial).start_obs=0;
    
    binary(trial).completed             = true;
    selected(trial).session             = 20200101;
    selected(trial).run                 = 1;
%     raw(trial).x_eye                    = NaN;
%     raw(trial).y_eye                    = NaN;

    saccades(trial).tar_rad             = 5;
    saccades(trial).fix_pos             = 0+1i;
    saccades(trial).n_targets           = n_targets;
    saccades(trial).cue_pos             = sequence(1,trial);
    saccades(trial).col_dim             = [0 0 0];
    %     saccades(trial).all_convex_sides    = [possible_convex_sides(sequence(2,trial)) possible_convex_sides(possible_shapes~=sequence(2,trial)) repmat(possible_convex_sides(shape_to_repeat(sequence(2,trial))),1,n_targets-numel(possible_shapes))];
    %     saccades(trial).all_convexities     = [possible_convexities(sequence(2,trial)) possible_convexities(possible_shapes~=sequence(2,trial)) repmat(possible_convexities(shape_to_repeat(sequence(2,trial))),1,n_targets-numel(possible_shapes))];
    
    second_shape=mod(sequence(2,trial)+2,4);
    if second_shape==0; second_shape=4; end;
    random_shape=Shuffle(find(~ismember([1,2,3,4],[sequence(2,trial) second_shape])));
    saccades(trial).all_convex_sides    = [possible_convex_sides(sequence(2,trial)) possible_convex_sides(second_shape) possible_convex_sides(random_shape(1))];
    saccades(trial).all_convexities     = [possible_convexities(sequence(2,trial))  possible_convexities(second_shape) possible_convexities(random_shape(1))];
    
    target_distribution     = all_target_distributions(sequence(5,trial),:);
    target_distribution     = [target_distribution(1) Shuffle(target_distribution(2:end))];
    match_position          = Shuffle(target_distribution(sign(real(target_distribution))==sign(sequence(3,trial))));
    distractor_position     = Shuffle(target_distribution(sign(real(target_distribution))~=sign(sequence(4,trial))));
    distractor_position     = distractor_position(distractor_position~=match_position(1));
    
    %saccades(trial).all_tar_pos         = [sequence(3,trial) Shuffle(possible_tar_pos(possible_tar_pos~=sequence(3,trial)))] ;
    saccades(trial).all_tar_pos         = [match_position(1) distractor_position(1) Shuffle(target_distribution(~ismember(target_distribution,match_position(1)) & ~ismember(target_distribution,distractor_position(1))))] ;
    saccades(trial).all_tar_pos         = saccades(trial).all_tar_pos(1:n_targets);
    
end


%behaviour generation

for trial=1:N_trials
    % resets
    target_inspected         =0;
    n=0;
    raw(trial).x_eye=[];
    raw(trial).y_eye=[];
    raw(trial).states=[];
    % cue assignment for confusion matrix
    if sign(saccades(trial).cue_pos)<0
        cue_letter='L';
    else
        cue_letter='R';
    end
    
    % current trial's target positions
    [A,B]=ismember(possible_tar_pos,saccades(trial).all_tar_pos);
    ordered_potential_targets= B(B~=0);
    current_target_preferences=target_preferences(A);
    
    
    % Exploration algorithm according to target_preferences (for inspecting next target) and Confusion
    % matrix (for selecting shapes)
    while true 
        random_hitrate=rand; % Importantly modified each iteration (= each new target inspection)
        n=n+1;
        pool_of_targets     =   ordered_potential_targets(ordered_potential_targets~=target_inspected); % exclude the target that was inspected in the previous trial
        target_probability  =   current_target_preferences(ordered_potential_targets~=target_inspected); % define probabilities for inspecting each potential target
        temp_targets_n      =   randsample(pool_of_targets,1,true,target_probability); % randomly choose a target to inspect next according to the target_probabilities.
        target_inspected    =   temp_targets_n(1);
        
        current_Confusion_matrix            = Confusion_matrix.(cue_letter);
        current_Confusion_RT_matrix         = Confusion_RT_matrix.(cue_letter);
        
        
        
        % calculate selected shape index (for confusion matrix)
        tar_sid             = strcmp(possible_convex_sides, saccades(trial).all_convex_sides(target_inspected));
        tar_con             = possible_convexities == saccades(trial).all_convexities(target_inspected);
        tar_sid_corr        = strcmp(possible_convex_sides, saccades(trial).all_convex_sides(1));
        tar_con_corr        = possible_convexities == saccades(trial).all_convexities(1);
        tar_shape_selected  = find(tar_sid & tar_con);
        tar_shape_correct   = find(tar_sid_corr & tar_con_corr);
        
        
        saccades(trial).targets_inspected(1,n) = target_inspected;
        saccades(trial).all_inspection_durations(1,n) = 0.5 + current_Confusion_RT_matrix(tar_shape_correct,tar_shape_selected) + ...
            Inspection_time_bias_L_R(1+(sign(real(saccades(trial).all_tar_pos(target_inspected)))==1)) + ...
            randn*inspection_time_distribution_factor;
        % this could be improved !! and combined with raw data!
        saccades(trial).all_inspection_intervals(1,n) = 0.5 + current_Confusion_RT_matrix(tar_shape_correct,tar_shape_selected) + ...
            Inspection_time_bias_L_R(1+(sign(real(saccades(trial).all_tar_pos(target_inspected)))==1)) + ...
            randn*inspection_time_distribution_factor;
        
        raw(trial).x_eye=[raw(trial).x_eye repmat(real(saccades(trial).all_tar_pos(target_inspected)),1,round((saccades(trial).all_inspection_durations(1,n) + saccade_duration)*sampling_rate))];
        raw(trial).y_eye=[raw(trial).y_eye repmat(imag(saccades(trial).all_tar_pos(target_inspected)),1,round((saccades(trial).all_inspection_durations(1,n) + saccade_duration)*sampling_rate))];
        raw(trial).states=[raw(trial).states repmat(13,1,round((saccades(trial).all_inspection_durations(1,n)+ saccade_duration)*sampling_rate))];
        
        % first inspected target
        if n==1
            saccades(trial).endpos=saccades(trial).all_tar_pos(target_inspected);
            saccades(trial).lat = 0.27 + RT_time_bias_L_R.(cue_letter)(1+(sign(real(saccades(trial).all_tar_pos(target_inspected)))==1)) + randn*inspection_time_distribution_factor;
            raw(trial).x_eye=[zeros(1,round(saccades(trial).lat*sampling_rate)) raw(trial).x_eye ];
            raw(trial).y_eye=[zeros(1,round(saccades(trial).lat*sampling_rate))  raw(trial).y_eye ];
            raw(trial).states=[ones(1,round(saccades(trial).lat*sampling_rate))*13  raw(trial).states];
        end
        
        if random_hitrate<=current_Confusion_matrix(tar_shape_correct,tar_shape_selected) % target selection criterion
            saccades(trial).all_inspection_durations(1,n) = 1;
            break
        end
        
    end
    
    % final assignments once target is selected
    if target_inspected == 1;
        binary(trial).success=true;
    else
        binary(trial).success=false;
    end
    
    saccades(trial).exploration_time        = sum(saccades(trial).all_inspection_durations) + numel(saccades(trial).all_inspection_durations)*saccade_duration;
    saccades(trial).tar_pos                 = saccades(trial).all_tar_pos(target_inspected);
    saccades(trial).selected_convex_sides   = saccades(trial).all_convex_sides{target_inspected};
    saccades(trial).selected_convexity      = saccades(trial).all_convexities(target_inspected);
end

out{1}.saccades=saccades;
out{1}.reaches=saccades;
out{1}.task=task;
out{1}.binary=binary;
out{1}.selected=selected;
out{1}.states=states;
out{1}.raw=raw;
out{1}.keys.calc_val.i_sample_rate=sampling_rate;
end



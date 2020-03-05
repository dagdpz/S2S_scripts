function [revealed_target_position,revealed_nonnan,target_position_log_idx]=get_revealed_target_positions_LR(...
    revealed_target_positions_in_order,revealed_targets_in_order,a,b_before,previous_target_position_log_idx,...
    revealed_target_position,revealed_nonnan,target_position_log_idx)
%[revealed_target_position,revealed_target_position_correct,revealed_nonnan,target_position_log_idx]=...

max_tar_revealed=size(revealed_target_positions_in_order,2);
all_n_positions=[1,2];
position_log_idx=true(size(all_n_positions));
if nargin<3
    revealed_target_position={};
    revealed_nonnan={};
    target_position_log_idx={};
b_before={};
     a=1;
     previous_target_position_log_idx=true(size(revealed_target_positions_in_order,1),1);
end

for target_counter=[all_n_positions(position_log_idx) max(all_n_positions)+1]
    
    b=[b_before,{target_counter}];
    if target_counter==max(all_n_positions)+1
        target_position_log_idx{a}                  =isnan(revealed_target_positions_in_order(:,a));
        nonnan_norm_per_pos_index                   =isnan(revealed_target_positions_in_order(:,a));
    else
        if target_counter==1
            target_position_log_idx{a}                  =real(revealed_target_positions_in_order(:,a))<=0 ;
            nonnan_norm_per_pos_index                   =any(real(revealed_target_positions_in_order)<=0,2);
        elseif target_counter==2
            target_position_log_idx{a}                  =real(revealed_target_positions_in_order(:,a))>0 ;
            nonnan_norm_per_pos_index                   =any(real(revealed_target_positions_in_order)>0,2);
    end 
    end
    revealed_target_position{a}(b{:})           =sum(target_position_log_idx{a} & previous_target_position_log_idx);
    %revealed_target_position_correct{a}(b{:})   =sum(target_position_log_idx{a} & previous_target_position_log_idx & revealed_targets_in_order(:,a)==1);
    revealed_nonnan{a}(b{:})                    =sum(~isnan(revealed_target_positions_in_order(:,a)) & nonnan_norm_per_pos_index);
    
    if a<max_tar_revealed
        c=a+1;
        next_target_position_log_idx=target_position_log_idx{a} & previous_target_position_log_idx;
        [revealed_target_position,revealed_nonnan,target_position_log_idx]=get_revealed_target_positions_LR(...
            revealed_target_positions_in_order,revealed_targets_in_order,c,b,next_target_position_log_idx,...
            revealed_target_position,revealed_nonnan,target_position_log_idx);
    end
    
end
end
function [revealed_target_position,revealed_nonnan,target_position_log_idx]=get_revealed_target_positions(...
    revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,...
    all_n_positions,calculate_only_first,a,b_before,position_log_idx,previous_target_position_log_idx,...
    revealed_target_position,revealed_nonnan,target_position_log_idx)
%[revealed_target_position,revealed_target_position_correct,revealed_nonnan,target_position_log_idx]=...

max_tar_revealed=size(revealed_target_positions_in_order,2);
if nargin<7
    revealed_target_position={};
    revealed_nonnan={};
    target_position_log_idx={};
    b_before={};
    a=1;
    previous_target_position_log_idx=true(size(revealed_target_positions_in_order,1),1);
    position_log_idx=true(size(all_n_positions));
end
for idx=1:numel(b_before)
    position_log_idx=position_log_idx & all_n_positions~=b_before{idx};
end
for target_counter=[all_n_positions(position_log_idx) max(all_n_positions)+1]
    
    b=[b_before,{target_counter}];
    if target_counter==max(all_n_positions)+1
        target_position_log_idx{a}                  =isnan(revealed_target_positions_in_order(:,a));
        nonnan_norm_per_pos_index                   =isnan(revealed_target_positions_in_order(:,a));
         revealed_nonnan{a}(b{:})                    =0;
    else
        target_position_log_idx{a}                  =revealed_target_positions_in_order(:,a)==possible_positions(target_counter) ;
        nonnan_norm_per_pos_index                   =any(revealed_target_positions_in_order==possible_positions(target_counter),2);
         revealed_nonnan{a}(b{:})                    =sum(sum(allpos==possible_positions(target_counter)));
    end
    revealed_target_position{a}(b{:})           =sum(target_position_log_idx{a} & previous_target_position_log_idx);
%    revealed_target_position_correct{a}(b{:})   =sum(target_position_log_idx{a} & previous_target_position_log_idx & revealed_targets_in_order(:,a)==1);
    %revealed_nonnan{a}(b{:})                    =sum(~isnan(revealed_target_positions_in_order(:,a)) & nonnan_norm_per_pos_index);
    
    if a<max_tar_revealed && ~calculate_only_first
        c=a+1;
        next_target_position_log_idx=target_position_log_idx{a} & previous_target_position_log_idx;
        [revealed_target_position,revealed_nonnan,target_position_log_idx]=get_revealed_target_positions(...
            revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,all_n_positions,...
            calculate_only_first,c,b,position_log_idx,next_target_position_log_idx,revealed_target_position,revealed_nonnan,target_position_log_idx);
    end
    
end
end
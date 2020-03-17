function M2S_output=S2S_session_analysis(out,settings)

%% trial selection according to type and effector, and if they are completed, also invisible ones and choice?

sampling_rate=out.keys.calc_val.i_sample_rate;
switch settings.type
    case 5; Disp.masked_or_not='not masked';
    case 6; Disp.masked_or_not='masked';
end
if ~settings.placeholders_visible
    Disp.masked_or_not='invisible';
end

switch settings.choice
    case 0; Disp.n_cues='one cue';
    case 1; Disp.n_cues='two cues';
end
log_idx.fixation_breaks=[out.states.state_abo]<=7 & [out.states.state_abo]>0;
log_idx.fix=[out.states.state_abo]==7 | [out.states.state_abo]==6;
log_idx.trials=~log_idx.fixation_breaks;


Sacpos_fix_break=NaN;
Delay_fix_break=NaN;
samplepos_fix_break=NaN;
if settings.effector == 0;          observed_structure=[out.saccades];   Disp.S_or_R='saccades';
    Raw_fix=out.raw(log_idx.fix);
    for n=1:sum(log_idx.fix)
        Sacpos_fix_break(n)=Raw_fix(n).x_eye(end);
        Delay_fix_break(n)=sum(Raw_fix(n).states>=6);
    end
    samplepos_fix_break           = [out.saccades(log_idx.fix).cue_pos];
elseif settings.effector == 1;      observed_structure=[out.reaches];    Disp.S_or_R='reaches';
end

colors={observed_structure.col_dim};
placeholders_visible=cellfun(@(x) any(x(:)>0), colors); % ?
observed_effector=[out.task.effector]==settings.effector & [out.task.type]==settings.type & [out.binary.choice]==settings.choice & placeholders_visible==settings.placeholders_visible;
observed_structure=observed_structure(log_idx.trials & observed_effector);
saccades=[out.saccades(log_idx.trials & observed_effector)];
% reaches=[out.reaches(log_idx.trials & observed_effector)];

%% Display
n_targets               =unique([observed_structure.n_targets]);
sessions                =[out.selected([out.binary.completed] & observed_effector).session];
runs                    =[out.selected([out.binary.completed] & observed_effector).run];

Disp.settings           = settings;
Disp.sessions_runs      = [num2str(sessions(1)) '_' num2str(runs(1)) ' - ' num2str(sessions(end)) '_' num2str(runs(end))];
Disp.n_targets_string   = num2str(n_targets);
Disp.tar_radius         = max([observed_structure.tar_rad]);

Disp.plot_only_LR=any(n_targets>5); %%!
%Disp.plot_only_LR=numel(unique(revealed_target_positions_in_order(~isnan(revealed_target_positions_in_order))))>5; %%!

%% Restructuring Positions and shapes per trial
binary          = out.binary(log_idx.trials & observed_effector);
selected        = out.selected(log_idx.trials & observed_effector);
raw             = out.raw(log_idx.trials & observed_effector);
states          = out.states(log_idx.trials & observed_effector);
completed_index = [binary.completed];
log_idx.success = [binary.success]';
log_idx.error   = ~[binary.success]';

% fix_y           = imag([observed_structure(1).fix_pos]);
% fix_x           = real([observed_structure(1).fix_pos]);
fix_y           = 0;
fix_x           = 0;
Disp.fix_y      = nanmean(fix_y);

insp                = {observed_structure.targets_inspected}; % order of inspected targets per trial (the MAIN parameter )
insp_dur            = {observed_structure.all_inspection_durations};
insp_int            = {observed_structure.all_inspection_intervals};
samplepos           = [observed_structure.cue_pos] - [observed_structure.fix_pos]; %%
Disp.samplepos      = num2str(min(abs(samplepos)));


allpos              = cell2mat_expanded({observed_structure.all_tar_pos}).';
allpos              = allpos - repmat([observed_structure.fix_pos], size(allpos,1), 1);
allconvexities      = cell2mat_expanded({observed_structure.all_convexities})';
allconvexsides      = cell2mat_expanded({observed_structure.all_convex_sides})';

chosenpos           = [observed_structure.tar_pos]  - [observed_structure.fix_pos]; %%
chosenconvexity     = [observed_structure.selected_convexity]';
chosensides         = {observed_structure.selected_convex_sides}';

matchpos            = allpos(1,:);
matchconvexity      = allconvexities(1,:);
matchsides          = allconvexsides(1,:);

log_idx.distpos=allpos~=repmat(matchpos, size(allpos,1), 1);

%% Defining indexes per trial
log_idx.sample_left     =real(samplepos)'<=fix_x;
log_idx.sample_right    =real(samplepos)'>fix_x;
log_idx.match_left      =real(matchpos)' <=fix_x;
log_idx.match_right     =real(matchpos)' >fix_x;

log_idx.all_right       =real(allpos)>fix_x;
log_idx.all_left        =real(allpos)<fix_x;

log_idx.more_left   = sum(real(allpos)<0,1)'>=sum(real(allpos)>0,1)';
log_idx.more_right  = sum(real(allpos)>0,1)'>=sum(real(allpos)<0,1)';

log_idx.RT.LL    =log_idx.sample_left  & real([observed_structure.endpos].')<fix_x;
log_idx.RT.LR    =log_idx.sample_left  & real([observed_structure.endpos].')>fix_x;
log_idx.RT.RL    =log_idx.sample_right & real([observed_structure.endpos].')<fix_x;
log_idx.RT.RR    =log_idx.sample_right & real([observed_structure.endpos].')>fix_x;

log_idx.FB.LL    =real(samplepos_fix_break)<=fix_x & real(Sacpos_fix_break)<fix_x;
log_idx.FB.LR    =real(samplepos_fix_break)<=fix_x & real(Sacpos_fix_break)>fix_x;
log_idx.FB.RL    =real(samplepos_fix_break)> fix_x & real(Sacpos_fix_break)<fix_x;
log_idx.FB.RR    =real(samplepos_fix_break)> fix_x & real(Sacpos_fix_break)>fix_x;

log_idx.dist.LL    =log_idx.sample_left  & any(log_idx.all_left & log_idx.distpos)';
log_idx.dist.LR    =log_idx.sample_left  & any(log_idx.all_right & log_idx.distpos)';
log_idx.dist.RL    =log_idx.sample_right  & any(log_idx.all_left & log_idx.distpos)';
log_idx.dist.RR    =log_idx.sample_right  & any(log_idx.all_right & log_idx.distpos)';

% possible convexities and positions per run (!!)
unique_runs=uniqueRowsCA(num2cell([[selected.session];[selected.run]])');
for run=1:size(unique_runs,1)
    run_idx=[selected.session]==unique_runs{run,1}&[selected.run]==unique_runs{run,2};
    allconvexities_run      = cell2mat_expanded({observed_structure(run_idx).all_convexities})';
    allconvexsides_run      = cell2mat_expanded({observed_structure(run_idx).all_convex_sides})';
    allshapes_run= [num2cell(allconvexities_run(~isnan(allconvexities_run))) ,allconvexsides_run(cellfun(@(x) ~isempty(x) && ~all(isnan(x)),allconvexsides_run)) ];
    unique_shapes=uniqueRowsCA(allshapes_run); %% relevant function uniqueRowsCA!!
    
    possible_convexities_per_run{run}    =vertcat(unique_shapes{:,1});
    possible_sides_per_run{run}          =unique_shapes(:,2);
    possible_positions_per_run{run}      =unique(allpos(~isnan(allpos)));
    %
    %     possible_convexities(:,run)    =vertcat(unique_shapes{:,1});
    %     possible_sides(:,run)          =unique_shapes(:,2);
    %     possible_positions(:,run)      =unique(allpos(~isnan(allpos)));
end

% sides and positions are still taken from run 1 only, asusming that they
% do not change across runs... this might be crucial for the future
possible_sides          =possible_sides_per_run{1};
possible_positions      =possible_positions_per_run{1};
all_n_positions         =1:size(possible_positions,1);

%% reordering positions
unique_abs_positions=unique(abs(possible_positions).*sign(real(possible_positions)));
poscounter=0;
for abspos=1:numel(unique_abs_positions)
    same_xpos=possible_positions(abs(possible_positions)==abs(unique_abs_positions(abspos)) & sign(real(possible_positions))==sign(unique_abs_positions(abspos)));
    [unique_y_positions,idx]=sort(imag(same_xpos));
    unique_x_positions=real(same_xpos(idx));
    unique_x_positions=flipud(unique_x_positions);
    unique_y_positions=flipud(unique_y_positions);
    for ypos=1:numel(unique_y_positions)
        poscounter=poscounter+1;
        reordered_positions(poscounter)=unique_x_positions(ypos) + 1i*unique_y_positions(ypos);
    end
end
possible_positions=reordered_positions.';

%% combining shapes (convexities)
if settings.concatinate_option==1 %% for difficulty adjustment
    all_possible_convexities=vertcat(possible_convexities_per_run{:});
    all_possible_sides=vertcat(possible_sides_per_run{:});
    unique_possible_convexities=unique(all_possible_convexities);
    if numel(unique(sign(all_possible_convexities)))~=1
        limit_convexity_N=min(abs(all_possible_convexities(sign(all_possible_convexities)==-1))); %%!!
        limit_convexity_P=min(abs(all_possible_convexities(sign(all_possible_convexities)==1))); %%!!
    elseif isfield(settings,'limit_convexity')
        limit_convexity_N=settings.limit_convexity; %%!!
        limit_convexity_P=settings.limit_convexity; %%!!
    else
        limit_convexity_N=min(abs(all_possible_convexities)); %%!!
        limit_convexity_P=min(abs(all_possible_convexities)); %%!!
    end
    if numel(unique(all_possible_sides))==1 && numel(unique(sign(all_possible_convexities)))==1 %% new special case, not sure if it's working        
        possible_convexities_tmp{1}=all_possible_convexities(all_possible_convexities==unique_possible_convexities(1));
        possible_convexities_tmp{2}=all_possible_convexities(all_possible_convexities==unique_possible_convexities(2));
        possible_convexities_tmp{3}=all_possible_convexities(all_possible_convexities==unique_possible_convexities(3));
        possible_convexities_tmp{4}=all_possible_convexities(all_possible_convexities==unique_possible_convexities(4));
        possible_sides={'LR';'LR';'LR';'LR'};
    elseif numel(unique(all_possible_sides))==1 
        possible_convexities_tmp{1}=all_possible_convexities(sign(all_possible_convexities)==-1 & (all_possible_convexities<  abs(limit_convexity_N)*-1));
        possible_convexities_tmp{2}=all_possible_convexities(sign(all_possible_convexities)==-1 & (all_possible_convexities>= abs(limit_convexity_N)*-1));
        possible_convexities_tmp{3}=all_possible_convexities(sign(all_possible_convexities)==1  & (all_possible_convexities<= abs(limit_convexity_P)));
        possible_convexities_tmp{4}=all_possible_convexities(sign(all_possible_convexities)==1  & (all_possible_convexities>  abs(limit_convexity_P)));
        possible_sides={'LR';'LR';'LR';'LR'};
    else
        possible_convexities_tmp{1}=all_possible_convexities(strcmp(all_possible_sides,'L')  & (all_possible_convexities <-limit_convexity_N | all_possible_convexities >limit_convexity_P));
        possible_convexities_tmp{2}=all_possible_convexities(strcmp(all_possible_sides,'LR') & (all_possible_convexities <-limit_convexity_N | all_possible_convexities >limit_convexity_P));
        possible_convexities_tmp{3}=all_possible_convexities(strcmp(all_possible_sides,'R')  & (all_possible_convexities <-limit_convexity_N | all_possible_convexities >limit_convexity_P));
        possible_convexities_tmp{4}=all_possible_convexities(strcmp(all_possible_sides,'LR') & (all_possible_convexities >=-limit_convexity_N & all_possible_convexities <=limit_convexity_P));
        possible_sides={'L';'LR';'R';'LR'};
    end
    for new_convexity_index=1:4
        if isempty(possible_convexities_tmp{new_convexity_index})
            possible_convexities_tmp{new_convexity_index}=NaN;
        end
        current_convexity=round(max(abs(possible_convexities_tmp{new_convexity_index}))*sign(possible_convexities_tmp{new_convexity_index}(1))*100)/100;
        allconvexities(ismember(allconvexities,possible_convexities_tmp{new_convexity_index}))      = current_convexity;
        chosenconvexity(ismember(chosenconvexity,possible_convexities_tmp{new_convexity_index}))    = current_convexity;
        matchconvexity(ismember(matchconvexity,possible_convexities_tmp{new_convexity_index}))      = current_convexity;
        possible_convexities(new_convexity_index)                                                   = current_convexity;
    end
    possible_convexities=possible_convexities';
else
    possible_convexities=possible_convexities_per_run{1};
end

%% Inspected targets and ordering revealed target positions, shapes, and inspection durations
for trial_idx=1:numel(insp) % Looping through every trial again!
    targets_inspected_tr=insp{trial_idx}(~isnan(insp{trial_idx}));
    current_inspection_durations=insp_dur{trial_idx}(~isnan(insp_dur{trial_idx}));
    current_inspection_intervals=insp_int{trial_idx};
    
    %     % only relevant for debugging ?!
    %     if isempty(current_inspection_durations)
    %         current_inspection_durations=NaN(size(targets_inspected_tr));
    %         fprintf([' fail: ' num2str(trial_idx)]);
    %     elseif binary(trial_idx).completed==1
    %         current_inspection_durations(end)=NaN;
    %     end
    
    %inspected_targets(trial_idx)                          =numel(targets_inspected_tr);
    [revealed_targets_tr, revealed_targets_tr_indexes]    =unique(targets_inspected_tr,'first');
    targets_ordered_tr                                    =targets_inspected_tr(sort(revealed_targets_tr_indexes));
    revealed_targets(trial_idx)                           =numel(revealed_targets_tr);
    
    if ~isempty(targets_inspected_tr)
        current_positions_L                                             =real(allpos(sub2ind(size(allpos),targets_inspected_tr,repmat(trial_idx,size(targets_inspected_tr)))))<=0;
        returns_to_revealed_per_tar_pos.left(trial_idx)                 =max(hist(targets_inspected_tr([true, diff(targets_inspected_tr)~=0] &  current_positions_L),revealed_targets_tr));
        returns_to_revealed_per_tar_pos.right(trial_idx)                =max(hist(targets_inspected_tr([true, diff(targets_inspected_tr)~=0] & ~current_positions_L),revealed_targets_tr));
        current_inspection_durations(targets_inspected_tr==targets_inspected_tr(end))=NaN;
    else % trials with no inspected target (No exploration)
        returns_to_revealed_per_tar_pos.left(trial_idx)=NaN;
        returns_to_revealed_per_tar_pos.right(trial_idx)=NaN;
    end
    
    for rev_idx=1:max(n_targets) % rev_idx indexes each (potentially revealed) target only ONCE (NOT in order of revealing, therefore always used as an index for targets_ordered_tr)
        if revealed_targets(trial_idx)>=rev_idx % if
            revealed_targets_in_order(trial_idx,rev_idx)            =targets_ordered_tr(rev_idx);
            revealed_target_positions_in_order(trial_idx,rev_idx)   =allpos(sub2ind(size(allpos),targets_ordered_tr(rev_idx),trial_idx));
            revealed_convexities_in_order(trial_idx,rev_idx)        =allconvexities(sub2ind(size(allconvexities),targets_ordered_tr(rev_idx),trial_idx));
            revealed_convexsides_in_order(trial_idx,rev_idx)        =allconvexsides(sub2ind(size(allconvexsides),targets_ordered_tr(rev_idx),trial_idx));
            inspection_time_per_target_total(trial_idx,rev_idx)     =nansum(current_inspection_durations(targets_inspected_tr==targets_ordered_tr(rev_idx)));
            inspection_time_per_target_max(trial_idx,rev_idx)       =max(current_inspection_durations(targets_inspected_tr==targets_ordered_tr(rev_idx)))*1000;
            temp_int=current_inspection_intervals(targets_inspected_tr==targets_ordered_tr(rev_idx));
            search_time_per_target(trial_idx,rev_idx)=temp_int(1)*1000;
        else % rev_idx indexes each target only ONCE (in order of revealing)
            revealed_targets_in_order(trial_idx,rev_idx)            =NaN;
            revealed_target_positions_in_order(trial_idx,rev_idx)   =NaN+1i*NaN;
            revealed_convexities_in_order(trial_idx,rev_idx)        =NaN;
            revealed_convexsides_in_order(trial_idx,rev_idx)        ={''};
            inspection_time_per_target_total(trial_idx,rev_idx)     =NaN;
            inspection_time_per_target_max(trial_idx,rev_idx)       =NaN;
            search_time_per_target(trial_idx,rev_idx)               =NaN;
        end
    end
    search_time_per_target(trial_idx,1)=search_time_per_target(trial_idx,1)-observed_structure(trial_idx).lat*1000; % removing reaction time
end

%% Exploration Pie chart calculation
[revealed_tar_pos.selected,revealed_tar_pos.nonnan]=DAG_get_revealed_target_positions(revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,all_n_positions,Disp.plot_only_LR);

%if Disp.plot_only_LR==0
%     [revealed_tar_pos.more_right,revealed_tar_pos.more_right_nonnan]=...
%         DAG_get_revealed_target_positions(revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,all_n_positions,false);
%
%     revealed_tar_pos_LR.more_right              =revealed_tar_pos.more_right;
%     revealed_tar_pos_LR.more_right_nonnan       =revealed_tar_pos.more_right_nonnan;

%     [revealed_tar_pos.more_left,revealed_tar_pos.more_left_nonnan]=...
%         DAG_get_revealed_target_positions(revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,all_n_positions,false);

%     revealed_tar_pos_LR.more_left               =revealed_tar_pos.more_left;
%     revealed_tar_pos_LR.more_left_nonnan        =revealed_tar_pos.more_left_nonnan ;
% else

%     revealed_tar_pos_LR.more_right              =NaN;
%     revealed_tar_pos_LR.more_right_nonnan       =NaN;
%
%     revealed_tar_pos_LR.more_left               =NaN;
%     revealed_tar_pos_LR.more_left_nonnan        =NaN;
%     [revealed_tar_pos.more_left,revealed_tar_pos.more_left_nonnan]=...
%         DAG_get_revealed_target_positions(revealed_target_positions_in_order,revealed_targets_in_order,allpos,possible_positions,all_n_positions,true);
%       [revealed_tar_pos.more_right,revealed_tar_pos.more_right_nonnan]=...
%         DAG_get_revealed_target_positions(revealed_target_positions_in_order(more_targets_right,:),revealed_targets_in_order(more_targets_right,:),allpos(:,more_targets_right),possible_positions,all_n_positions,1);
%
% end



%% Exploration & Inspection  times
% Per revealed target & cue position

%Bins.exploration=0:0.2:5;
Bins.reaction_time=0:50:1000;
Bins.fixation_breaks=0:50:1000;
%Bins.inspection_total=0:50:1000;
Bins.inspection_mean=0:50:1000;
Bins.search_time=0:50:1000;
Bins.exploration_time=0:200:5000;
Bins.n_inspected=0:1:10;
Bins.n_return=0:1:5;
Bins.stepsize=0.5;
Bins.raw_2d{1}=(nanmean(fix_x)-30):Bins.stepsize:(nanmean(fix_x)+30);
Bins.raw_2d{2}=(nanmean(fix_y)-20):Bins.stepsize:(nanmean(fix_y)+20);


RT_subfieldnames={'LL','LR','RL','RR'};
for FN=RT_subfieldnames
    RTs.raw.(FN{:})             = [observed_structure(log_idx.RT.(FN{:})).lat, single(NaN)]*1000; % why is this a single??
    RTs.median.(FN{:})          = nanmedian(RTs.raw.(FN{:}));
    RTs.histograms.(FN{:})      = hist(RTs.raw.(FN{:}),Bins.reaction_time);
    RTs.histograms_fb.(FN{:})   = hist(Delay_fix_break(log_idx.FB.(FN{:})),Bins.fixation_breaks);
    RTs.median_fb.(FN{:})       = nanmedian(Delay_fix_break(log_idx.FB.(FN{:})));
    revealed_tar_pos.first.(FN{:}) = sum(log_idx.RT.(FN{:}));
end

for k=1:numel(possible_positions)
    inspection_time.sum_per_position{1}(k)=nansum(inspection_time_per_target_total(revealed_target_positions_in_order==possible_positions(k)));
    [~, inspection_time.spot_in_revealing_sequence{k}]=ind2sub(size(revealed_target_positions_in_order),find(revealed_target_positions_in_order==possible_positions(k)));
    inspection_time.mean_spot_in_revealing_sequence(k)=nanmean(inspection_time.spot_in_revealing_sequence{k});
    %revealed_tar_pos.success_rate_per_target{1}(k) = sum([binary(matchpos==possible_positions(k)).success])/sum(matchpos==possible_positions(k));
    revealed_tar_pos.success_rate_per_target{1}(k) = sum([binary(chosenpos==possible_positions(k)).success])/sum(chosenpos==possible_positions(k));
end

L_R_sign=[-1,1,-1,1];
FN={'LL','LR','RL','RR'};
SFN={'sample_left','sample_left','sample_right','sample_right'};
for k=1:numel(FN)
    position_k.(FN{k})                                  = repmat(log_idx.(SFN{k}) ,1,max(n_targets)) & sign(real(revealed_target_positions_in_order))==L_R_sign(k);
    inspection_time.max_per_position.(FN{k})            = inspection_time_per_target_max(position_k.(FN{k}));
    inspection_time.median_max_per_position.(FN{k})     = nanmedian(inspection_time.max_per_position.(FN{k}));
    inspection_time.hist_per_position.(FN{k})           = hist(inspection_time.max_per_position.(FN{k}),Bins.inspection_mean);
end


for k=1:size(possible_convexities,1)
    log_idx_match{k}=ismember(matchconvexity',possible_convexities(k,1)) & strcmp(matchsides',possible_sides(k)) & completed_index';
    for m=1:size(possible_convexities,1)
        shape{k,m}                                                          = repmat(log_idx_match{k},1,max(n_targets)) & ismember(revealed_convexities_in_order,possible_convexities(m,:)) & strcmp(revealed_convexsides_in_order,possible_sides(m));
        inspection_time.max_per_convexity{k,m}(1:sum(sum(shape{k,m})))     = inspection_time_per_target_max(shape{k,m});
        inspection_time.hist_per_convexity{k,m}                        = hist(inspection_time.max_per_convexity{k,m},Bins.inspection_mean);
        inspection_time.median_max_per_convexity{k,m}                      = nanmedian(inspection_time.max_per_convexity{k,m});
        
        log_idx_chosen{m}              =ismember(chosenconvexity,possible_convexities(m,:)) & strcmp(chosensides,possible_sides(m));
        log_idx_current_shape{m}       =ismember(allconvexities,possible_convexities(m,:)) & strcmp(allconvexsides,possible_sides(m));
        
        log_idx_shape_left      = any(log_idx_current_shape{m} & log_idx.all_left,1)';
        log_idx_shape_right     = any(log_idx_current_shape{m} & log_idx.all_right,1)';
        
        Sel.LL{k,m}=sum(log_idx_match{k} & log_idx.sample_left  & log_idx_shape_left  & real(chosenpos.')<0 & log_idx_chosen{m});
        Sel.LR{k,m}=sum(log_idx_match{k} & log_idx.sample_left  & log_idx_shape_right & real(chosenpos.')>0 & log_idx_chosen{m});
        Sel.RL{k,m}=sum(log_idx_match{k} & log_idx.sample_right & log_idx_shape_left  & real(chosenpos.')<0 & log_idx_chosen{m});
        Sel.RR{k,m}=sum(log_idx_match{k} & log_idx.sample_right & log_idx_shape_right & real(chosenpos.')>0 & log_idx_chosen{m});
        
        Tot.LL{k,m}=sum(log_idx_match{k} & log_idx.sample_left  & log_idx_shape_left);
        Tot.LR{k,m}=sum(log_idx_match{k} & log_idx.sample_left  & log_idx_shape_right);
        Tot.RL{k,m}=sum(log_idx_match{k} & log_idx.sample_right & log_idx_shape_left);
        Tot.RR{k,m}=sum(log_idx_match{k} & log_idx.sample_right & log_idx_shape_right);
        
        Tot.all{k,m}=sum(log_idx_match{k}  & (log_idx_shape_right | log_idx_shape_left));
        
        Summary.confusion.Sel{k,m}=Sel.LL{k,m}+Sel.LR{k,m}+Sel.RL{k,m}+Sel.RR{k,m};
        Summary.confusion.Tot{k,m}=Tot.all{k,m};
        
        Confusion_matrix.shapes{k,m}={possible_convexities(k,1),possible_sides{k,1};possible_convexities(m,1),possible_sides{m,1}}; % k is match, m is selected
    end
end



%% object_related_neglect indexes

OBJ_FN={'obj_L','obj_R','obj_LR'};
L_or_R={'R','LR','L';'L','LR','R';'L','R','LR'};
reference_convexity=min(min((abs(possible_convexities))))*sign(sum(sum(possible_convexities)));
for idx=1:numel(OBJ_FN)
    FN=['idx_' OBJ_FN{idx}];
    for conf_idx= 1:numel(Confusion_matrix.shapes)
        Confusion_matrix.(FN)(conf_idx)=  (Confusion_matrix.shapes{conf_idx}{1,1}==Confusion_matrix.shapes{conf_idx}{2,1} ...
            &&  any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,1})) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,2})))...
            || (any([Confusion_matrix.shapes{conf_idx}{:,1}]'==reference_convexity) && any([Confusion_matrix.shapes{conf_idx}{:,1}]'~=reference_convexity & strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,3})));
        %    Confusion_matrix.idx_obj_L(conf_idx)=  (Confusion_matrix.shapes{conf_idx}{1,1}==Confusion_matrix.shapes{conf_idx}{2,1} ...
        %                                              &&  any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'R')) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'LR')))...
        %                                              || (any([Confusion_matrix.shapes{conf_idx}{:,1}]==0) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'L')));
    end
    
    Confusion_matrix.(FN)=reshape(Confusion_matrix.(FN),size(possible_convexities,1),size(possible_convexities,1));
    
    Summary.confusion.(['Sel_' OBJ_FN{idx}])=sum([Summary.confusion.Sel{Confusion_matrix.(FN)}]);
    Summary.confusion.(['Tot_' OBJ_FN{idx}])=sum([Summary.confusion.Tot{Confusion_matrix.(FN)}]);
    Summary.confusion.(['dwell_time_' OBJ_FN{idx}])=[inspection_time.max_per_convexity{Confusion_matrix.(FN)}];
    Summary.confusion.(['mean_dwell_time_' OBJ_FN{idx}])=nanmean(Summary.confusion.(['dwell_time_' OBJ_FN{idx}]));
    
end

%% Remaining summaries

log_idx.match_revealed                  =any(revealed_targets_in_order==1,2);
for   rev_idx=1:max(n_targets)
    Summary.total_revealed(rev_idx)                         = sum(revealed_targets==rev_idx);
    Summary.completed_revealed(rev_idx)                     = sum(revealed_targets==rev_idx & completed_index);
    Summary.completed_revealed_including_match(rev_idx)     = sum(revealed_targets==rev_idx & completed_index & log_idx.match_revealed');
    Summary.successful_revealed(rev_idx)                    = sum(revealed_targets==rev_idx & [binary.success] & completed_index);
end

for tar=1:max(n_targets)
    inspection_time.mean_search_to_reveal(tar)=nanmean(search_time_per_target(:,tar));
    inspection_time.hist_search_to_reveal{tar}=hist(search_time_per_target(:,tar),Bins.search_time);
end

et_suc=[observed_structure([binary.success]).exploration_time].*1000;
et_err=[observed_structure(~[binary.success]).exploration_time].*1000;

Summary.successrate_over_time               = hist(et_suc,Bins.exploration_time)./(hist(et_suc,Bins.exploration_time) + hist(et_err,Bins.exploration_time));
%Bins.exploration_time

Summary.First                               = revealed_tar_pos.first;
Summary.RTs                                 = RTs.raw;
Summary.inspection_time_per_trial           = inspection_time.max_per_position;
Summary.inspection_time_median              = inspection_time.median_max_per_position;

Summary.TimeOut.Match_revealed              = sum(~completed_index' & log_idx.match_revealed);
Summary.TimeOut.None_revealed               = sum(all(isnan(revealed_targets_in_order),2));
Summary.TimeOut.All                         = sum(~completed_index');
Summary.All_match_revealed_completed        = sum(completed_index' & log_idx.match_revealed);
Summary.Match_revealed_but_not_selected     = sum(completed_index' & log_idx.match_revealed & ~[binary.success]');
Summary.All_success                         = sum([binary.success]==1);
Summary.cue_L_success                       = sum([binary.success]==1 & log_idx.sample_left') ;
Summary.cue_R_success                       = sum([binary.success]==1 & log_idx.sample_right') ;
Summary.All                                 = numel(binary);
Summary.Completed                           = sum([binary.completed]==1);
Summary.Fixation_breaks                     = sum(log_idx.fix);
Summary.cue_L                               = sum(log_idx.sample_left);
Summary.cue_R                               = sum(log_idx.sample_right);
Summary.Succes_rate                         = Summary.All_success/Summary.All;
Summary.Succes_rate_cue_L                   = Summary.cue_L_success/Summary.cue_L;
Summary.Succes_rate_cue_R                   = Summary.cue_R_success/Summary.cue_R;



SEL_FN      = {'RL','RR','LL','LR'};
samp_FN     = {'sample_right','sample_right','sample_left','sample_left'};
match_FN    = {'match_left','match_right','match_left','match_right'};
tar_FN      = {'left','right','left','right'};
for k=1:numel(SEL_FN)
    FN=SEL_FN(k);
    Summary.Fix.(FN{:})                     = sum(log_idx.FB.(FN{:}));
    %Summary.Fixdel.(FN{:})                  = nanmean(Delay_fix_break(log_idx.fix.(FN{:})));
    Summary.TimeOut.(FN{:})                 = sum(~completed_index' & log_idx.(samp_FN{k}) & log_idx.(match_FN{k}) );
    Summary.Hit.(FN{:})                     = sum([Sel.(FN{:}){eye(size(possible_convexities,1))==1}]);
    Summary.Tot.(FN{:})                     = sum([Tot.(FN{:}){eye(size(possible_convexities,1))==1}]);
    Summary.Err.(FN{:})                     = sum([Sel.(FN{:}){eye(size(possible_convexities,1))==0}]);
    %Summary.TWr.(FN{:})                     = sum([Tot.(FN{:}){eye(size(possible_convexities,1))==0}]); %% Normalization...?
    Summary.TWr.(FN{:})                     = sum(log_idx.dist.(FN{:})); %% Normalization...?
    
    returns_to_revealed_per_tar_pos.(FN{:}) = returns_to_revealed_per_tar_pos.(tar_FN{k})(log_idx.(samp_FN{k}));
    returns_to_revealed.hist.(FN{:})        = hist(returns_to_revealed_per_tar_pos.(FN{:}) ,Bins.n_return);
    returns_to_revealed.mean.(FN{:})        = nanmean(returns_to_revealed_per_tar_pos.(FN{:}));
    %Summary.Notexplored.(FN{:})             = sum(returns_to_revealed_per_tar_pos.(FN{:})==0);
    Summary.Returns.(FN{:})                 = sum(returns_to_revealed_per_tar_pos.(FN{:})>1);
    Summary.Total_inspections.(FN{:})       = numel(returns_to_revealed_per_tar_pos.(FN{:}));
end




%% raw eye traces for heat map calculation

observed_states=[11,12,13,14];
if ~all(isnan(raw(1).x_eye))
    n_rows=numel(raw);
    % during ITI
    n_samples=arrayfun(@(x) sum(ismember(x.states,[50,1])), raw);
    n_columns=max(n_samples);
    ITI_x_eye=NaN(n_rows,n_columns);ITI_y_eye=NaN(n_rows,n_columns);
    % during trial
    n_samples=arrayfun(@(x) sum(ismember(x.states,observed_states)), raw);
    n_columns=max(n_samples);
    temp_x_eye=NaN(n_rows,n_columns);temp_y_eye=NaN(n_rows,n_columns);
    %     quantile_factor=3/4;
    %     n_colors=ceil(quantile(n_samples,quantile_factor));
    %n_colors=ceil(n_columns/2);
    %     time_factor=jet(n_colors);
    %     time_factor=[flipud(time_factor);repmat([0 0 1],n_columns-n_colors,1)];
    
    %     for tf=1:n_columns
    %         time_factor(tf,:)=time_factor(tf,:)/sum(time_factor(tf,:));
    %     end
    
    for row=1:n_rows
        trace_index=ismember(raw(row).states,observed_states) & (1:numel(raw(row).states))/sampling_rate > states(row).start_obs+saccades(row).lat;
        %trace_index=ismember(raw(row).states,observed_states) & abs(raw(row).x_eye -fix_x + 1i*(raw(row).y_eye -fix_y) )>4;
        n_cell_elements=sum(trace_index);
        temp_x_eye(row,1:n_cell_elements)=raw(row).x_eye(trace_index) - real(observed_structure(row).fix_pos);
        temp_y_eye(row,1:n_cell_elements)=raw(row).y_eye(trace_index) - imag(observed_structure(row).fix_pos);
        
        
        trace_index=ismember(raw(row).states,[50,1]);
        %trace_index=ismember(raw(row).states,observed_states) & abs(raw(row).x_eye -fix_x + 1i*(raw(row).y_eye -fix_y) )>4;
        n_cell_elements=sum(trace_index);
        ITI_x_eye(row,1:n_cell_elements)=raw(row).x_eye(trace_index) - real(observed_structure(row).fix_pos);
        ITI_y_eye(row,1:n_cell_elements)=raw(row).y_eye(trace_index) - imag(observed_structure(row).fix_pos);
    end
end

%% raw eye traces for heat map
Summary.explore=struct('LL',{NaN},'LR',{NaN},'RL',{NaN},'RR',{NaN});
if ~all(isnan(raw(1).x_eye))
    n_SL=0; n_SR=0;
    for n=1:numel(raw)
        
        trace_index=ismember(raw(n).states,observed_states) & (1:numel(raw(n).states))/sampling_rate > states(n).start_obs+saccades(n).lat;
        ITI_index=ismember(raw(n).states,[50,1]);
        if samplepos(n)<=0
            n_SL=n_SL+1;
            Summary.explore.LL(n_SL) = round(1000/sampling_rate)*sum(trace_index & raw(n).x_eye<0)';
            Summary.explore.LR(n_SL) = round(1000/sampling_rate)*sum(trace_index & raw(n).x_eye>0)';
            Summary.explore.total(n) = Summary.explore.LL(n_SL)+Summary.explore.LR(n_SL);
            Summary.ITIexpl.LL(n_SL) = round(1000/sampling_rate)*sum(ITI_index & raw(n).x_eye<0)';
            Summary.ITIexpl.LR(n_SL) = round(1000/sampling_rate)*sum(ITI_index & raw(n).x_eye>0)';
            Summary.ITIexpl.total(n) = Summary.ITIexpl.LL(n_SL)+Summary.ITIexpl.LR(n_SL);
            Summary.Notexplored.LL(n_SL)             = ~any(trace_index & raw(n).x_eye<0);
            Summary.Notexplored.LR(n_SL)             = ~any(trace_index & raw(n).x_eye>0);
    
            
        else
            n_SR=n_SR+1;
            Summary.explore.RL(n_SR) = round(1000/sampling_rate)*sum(trace_index & raw(n).x_eye<0)';
            Summary.explore.RR(n_SR) = round(1000/sampling_rate)*sum(trace_index & raw(n).x_eye>0)';
            Summary.explore.total(n) = Summary.explore.RL(n_SR)+Summary.explore.RR(n_SR);
            Summary.ITIexpl.RL(n_SR) = round(1000/sampling_rate)*sum(ITI_index & raw(n).x_eye<0)';
            Summary.ITIexpl.RR(n_SR) = round(1000/sampling_rate)*sum(ITI_index & raw(n).x_eye>0)';
            Summary.ITIexpl.total(n) = Summary.ITIexpl.RL(n_SR)+Summary.ITIexpl.RR(n_SR);
            Summary.Notexplored.RL(n_SR)             = ~any(trace_index & raw(n).x_eye<0);
            Summary.Notexplored.RR(n_SR)             = ~any(trace_index & raw(n).x_eye>0);
        end
    end
    Summary.Notexplored.LL=sum(Summary.Notexplored.LL);
    Summary.Notexplored.LR=sum(Summary.Notexplored.LR);
    Summary.Notexplored.RL=sum(Summary.Notexplored.RL);
    Summary.Notexplored.RR=sum(Summary.Notexplored.RR);
    
    temp_x_success=temp_x_eye(log_idx.success,:);
    temp_y_success=temp_y_eye(log_idx.success,:);
    temp_x_error=temp_x_eye(log_idx.error,:);
    temp_y_error=temp_y_eye(log_idx.error,:);
    [inspection_time.raw,~] = hist3([temp_x_eye(~isnan(temp_x_eye)),temp_y_eye(~isnan(temp_y_eye))],'Edges',Bins.raw_2d);
    [inspection_time.heat_success,~] = hist3([temp_x_success(~isnan(temp_x_success)),temp_y_success(~isnan(temp_y_success))],'Edges',Bins.raw_2d);
    [inspection_time.heat_error,~] = hist3([temp_x_error(~isnan(temp_x_error)),temp_y_error(~isnan(temp_y_error))],'Edges',Bins.raw_2d);
    [inspection_time.ITI,~] = hist3([ITI_x_eye(~isnan(ITI_x_eye)),ITI_y_eye(~isnan(ITI_y_eye))],'Edges',Bins.raw_2d);
    
else
    Summary.explore.total=NaN;
end


%%Jackson Pollock plot

% log_idx.success = [binary.success]';
% log_idx.error   = ~[binary.success]';
pattern_labels={'sample_left','more_left','success';...
    'sample_left','more_right','success';...
    'sample_right','more_left','success';...
    'sample_right','more_right','success';...
    'sample_left','more_left','error';...
    'sample_left','more_right','error';...
    'sample_right','more_left','error';...
    'sample_right','more_right','error'};
Disp.pattern_labels={'SUCCESS: Sample left','More targets left';...
    'SUCCESS: Sample left','More targets right';...
    'SUCCESS: Sample right','More targets left';...
    'SUCCESS: Sample right','More targets right';...
    'ERROR: Sample left','More targets left';
    'ERROR: Sample left','More targets right';
    'ERROR: Sample right','More targets left';
    'ERROR: Sample right','More targets right'};

if ~all(isnan(raw(1).x_eye))
    quantile_factor=3/4;
    n_colors=ceil(quantile(n_samples,quantile_factor));
    time_factor=jet(n_colors);
    time_factor=[flipud(time_factor);repmat([0 0 1],n_columns-n_colors,1)];
    
    %     for tf=1:n_columns
    %         time_factor(tf,:)=time_factor(tf,:)/sum(time_factor(tf,:));
    %     end
    
    %     R=0;G=0;B=0;
    %     for c=1:n_columns
    %         R=R+hist3([temp_x_eye(:,c),temp_y_eye(:,c)],'Edges',Bins.raw_2d).*time_factor(c,1);
    %         G=G+hist3([temp_x_eye(:,c),temp_y_eye(:,c)],'Edges',Bins.raw_2d).*time_factor(c,2);
    %         B=B+hist3([temp_x_eye(:,c),temp_y_eye(:,c)],'Edges',Bins.raw_2d).*time_factor(c,3);
    %     end
    %     norm_M=NaN(size(R));
    %     for n=1:numel(R)
    %         norm_M(n)=max([R(n),G(n),B(n)]);
    %     end
    %     Rn=R./norm_M;
    %     Gn=G./norm_M;
    %     Bn=B./norm_M;
    for k=1:8
        T_hist=0;T_hist_weigthed=0;
        idx_pattern=log_idx.(pattern_labels{k,1}) & log_idx.(pattern_labels{k,2}) & log_idx.(pattern_labels{k,3});
        Disp.pattern_target_positions{k}=unique(allpos(:,idx_pattern));
        for c=1:n_columns % This is the part that takes time....
            T_hist=T_hist+hist3([temp_x_eye(idx_pattern,c),temp_y_eye(idx_pattern,c)],'Edges',Bins.raw_2d);
            T_hist_weigthed=T_hist_weigthed+hist3([temp_x_eye(idx_pattern,c),temp_y_eye(idx_pattern,c)],'Edges',Bins.raw_2d).*c;
        end
        
        T_hist_average=round(T_hist_weigthed./T_hist);
        
        Rn=NaN(size(T_hist_average));
        Gn=NaN(size(T_hist_average));
        Bn=NaN(size(T_hist_average));
        for n=1:numel(Rn)
            if isnan(T_hist_average(n))
                continue
            end
            Rn(n)=time_factor(T_hist_average(n),1);
            Gn(n)=time_factor(T_hist_average(n),2);
            Bn(n)=time_factor(T_hist_average(n),3);
        end
        
        Pattern{k}(:,:,1) =rot90(Rn);
        Pattern{k}(:,:,2) =rot90(Gn);
        Pattern{k}(:,:,3) =rot90(Bn);
    end
    Summary.max_pattern_samples=n_columns*quantile_factor;
else
    Pattern=NaN;
    Summary.max_pattern_samples=NaN;
end

M2S_output=struct('Summary',Summary,'Sel',Sel,'Tot',Tot,'inspection_time',inspection_time,'revealed_tar_pos',revealed_tar_pos,...
    'possible_convexities',possible_convexities(:,1),'possible_positions',possible_positions,'Disp',Disp,'RTs',RTs,'Bins',Bins,'returns_to_revealed',returns_to_revealed);
M2S_output.possible_sides=possible_sides;
M2S_output.Pattern=Pattern;
end


function [Sum_ses,Sum_norm,Sum_tot]=S2S_group_summary_UMG(M2S_output,plot_specs)
Summary                 =[M2S_output.Summary];

%% Normalization per subject (?) && Statistics && Mean + STD

FN_in={'LL','LR';'RL','RR';'LL','RL';'LR','RR'};
FN_sum={'SL','SR','TL','TR'};
FN_first={'LL','RR';'RL','LR';'LL','RL';'LR','RR'};

OBJ_FN={'obj_L','obj_R','obj_LR'};
% for idx=1:numel(OBJ_FN)
%     FN=['idx_' OBJ_FN{idx}];
% for conf_idx= 1:numel(Confusion_matrix.shapes)
%    Confusion_matrix.(FN)(conf_idx)=  (Confusion_matrix.shapes{conf_idx}{1,1}==Confusion_matrix.shapes{conf_idx}{2,1} ...
%                                              &&  any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,1})) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,2})))...
%                                              || (any([Confusion_matrix.shapes{conf_idx}{:,1}]'==0) && any([Confusion_matrix.shapes{conf_idx}{:,1}]'~=0 & strcmp(Confusion_matrix.shapes{conf_idx}(:,2),L_or_R{idx,3})));
% %    Confusion_matrix.idx_obj_L(conf_idx)=  (Confusion_matrix.shapes{conf_idx}{1,1}==Confusion_matrix.shapes{conf_idx}{2,1} ...
% %                                              &&  any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'R')) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'LR')))...
% %                                              || (any([Confusion_matrix.shapes{conf_idx}{:,1}]==0) && any(strcmp(Confusion_matrix.shapes{conf_idx}(:,2),'L')));
% end
%
% Confusion_matrix.(FN)=reshape(Confusion_matrix.(FN),size(possible_convexities,1),size(possible_convexities,1));
%
% Summary.confusion.(['Sel_' OBJ_FN{idx}])=sum([Summary.confusion.Sel{Confusion_matrix.(FN)}]);
% Summary.confusion.(['Tot_' OBJ_FN{idx}])=sum([Summary.confusion.Tot{Confusion_matrix.(FN)}]);
% Summary.confusion.(['dwell_time_' OBJ_FN{idx}])=[inspection_time.mean_per_convexity{Confusion_matrix.(FN)}];
% Summary.confusion.(['mean_dwell_time_' OBJ_FN{idx}])=nanmean(Summary.confusion.(['dwell_time_' OBJ_FN{idx}]));
%
% end

%  Sum_norm(p,g).RTs.(FN{:}) = cellfun(@(x) nanmean(x),[Sum_ses(p,g).RTs.(FN{:})]);
%             Sum_norm(p,g).DTs.(FN{:}) = cellfun(@(x) nanmean(double(x)),[Sum_ses(p,g).DTs.(FN{:})]); %why is this not double?
%             Sum_norm(p,g).Raw.(FN{:}) = cellfun(@(x) nanmean(x),[Sum_ses(p,g).Raw.(FN{:})]);

for n=1:numel(Summary) % n: subject
    % object related problems
    for k=1:numel(OBJ_FN)
        
        Sum_ses.(['Sel_' OBJ_FN{k}])(n)=Summary(n).confusion.(['Sel_' OBJ_FN{k}]);
        Sum_ses.(['Tot_' OBJ_FN{k}])(n)=Summary(n).confusion.(['Tot_' OBJ_FN{k}]);
        Sum_ses.(['dwell_time_' OBJ_FN{k}]){n}=Summary(n).confusion.(['dwell_time_' OBJ_FN{k}]);
        Sum_ses.(['mean_dwell_time_' OBJ_FN{k}])(n)=Summary(n).confusion.(['mean_dwell_time_' OBJ_FN{k}]);
        
        mean_DT=nanmean([Summary(n).confusion.mean_dwell_time_obj_L,Summary(n).confusion.mean_dwell_time_obj_R,Summary(n).confusion.mean_dwell_time_obj_LR]);
        
        Sum_norm.Hit.(OBJ_FN{k})(n)=Sum_ses.(['Sel_' OBJ_FN{k}])(n)/Sum_ses.(['Tot_' OBJ_FN{k}])(n);
        Sum_norm.DTs.(OBJ_FN{k})(n)=Sum_ses.(['mean_dwell_time_' OBJ_FN{k}])(n)-mean_DT;
        
        
        Sum_stat.Hit.(OBJ_FN{k})(n)=Sum_ses.(['Sel_' OBJ_FN{k}])(n)/Sum_ses.(['Tot_' OBJ_FN{k}])(n);
        Sum_stat.DTs.(OBJ_FN{k})(n)=Sum_ses.(['mean_dwell_time_' OBJ_FN{k}])(n)-mean_DT;
        %Sum_norm.DTs.(OBJ_FN{k})(n)=Sum_ses.(['mean_dwell_time_' OBJ_FN{k}])(n)-mean_DT;
        
    end
    %Sum_ref.DTs
    
    
    % Hitrates
    for k=1:numel(FN_sum)
        Sum_ses.Hit.(FN_sum{k})(n)=Summary(n).Hit.(FN_in{k,1}) + Summary(n).Hit.(FN_in{k,2});
        %Sum_ses.Tot.(FN_sum{k})(n)=Summary(n).Tot.(FN_in{k,1}) + Summary(n).Tot.(FN_in{k,2});
        Sum_tot.Hit.(FN_sum{k})(n)=Summary(n).Tot.(FN_in{k,1}) + Summary(n).Tot.(FN_in{k,2});
        Sum_norm.Hit.(FN_sum{k})(n)=Sum_ses.Hit.(FN_sum{k})(n)/Sum_tot.Hit.(FN_sum{k})(n);
        Sum_stat.Hit.(FN_sum{k})(n)=Sum_norm.Hit.(FN_sum{k})(n);
    end
    Sum_ses.Hit.p_Sample(n)=DAG_fisher_exakt_scalars(Sum_ses.Hit.SL(n),Sum_tot.Hit.SL(n)-Sum_ses.Hit.SL(n),Sum_ses.Hit.SR(n),Sum_tot.Hit.SR(n)-Sum_ses.Hit.SR(n));
    Sum_ses.Hit.p_Target(n)=DAG_fisher_exakt_scalars(Sum_ses.Hit.TL(n),Sum_tot.Hit.TL(n)-Sum_ses.Hit.TL(n),Sum_ses.Hit.TR(n),Sum_tot.Hit.TR(n)-Sum_ses.Hit.TR(n));
    
    Sum_ref.Hit(n)=Summary(n).Completed;
    Sum_tot.Hit.all(n)=Summary(n).Completed;
    Sum_norm.Hit.all(n)=(Sum_ses.Hit.SL(n)+Sum_ses.Hit.SR(n))/Summary(n).Completed;
    Sum_ses.Hit.all(n)=Sum_ses.Hit.SL(n)+Sum_ses.Hit.SR(n);
    
    % First response
    
    Sum_ref.First(n)=Summary(n).First.LL + Summary(n).First.LR + Summary(n).First.RL + Summary(n).First.RR;
    Sum_ses.First.all(n)=Summary(n).First.LL + Summary(n).First.LR + Summary(n).First.RL + Summary(n).First.RR;
    Sum_tot.First.all(n)=Summary(n).All;
    for k=1:numel(FN_sum)
        Sum_ses.First.(FN_sum{k})(n)=Summary(n).First.(FN_first{k,1}) + Summary(n).First.(FN_first{k,2}); %% note that SL means Same sidde, SR means opposite !!!
    end
    for k=1:numel(FN_sum)
        m=mod(k,2) + 2*ceil(k/2)-1;
        Sum_tot.First.(FN_sum{k})(n)=Sum_ses.First.(FN_sum{k})(n)+Sum_ses.First.(FN_sum{m})(n); %% note that SL means Same sidde, SR means opposite !!!
        Sum_norm.First.(FN_sum{k})(n)=Sum_ses.First.(FN_sum{k})(n)/Sum_tot.First.(FN_sum{k})(n);
        Sum_stat.First.(FN_sum{k})(n)=Sum_norm.First.(FN_sum{k})(n);
    end
    %     Sum_ses.First.p_Sample(n)=DAG_fisher_exakt_scalars(Sum_ses.First.SL(n),Sum_ses.First.SR(n),Sum_ses.First.SR(n),Sum_ses.First.SL(n));
    %     Sum_ses.First.p_Target(n)=DAG_fisher_exakt_scalars(Sum_ses.First.TL(n),Sum_ses.First.TR(n),Sum_ses.First.TR(n),Sum_ses.First.TL(n));
    Sum_ses.First.p_Sample(n)=DAG_fisher_exakt_scalars(Sum_ses.First.SL(n),Sum_ref.First(n)-Sum_ses.First.SL(n),Sum_ses.First.SR(n),Sum_ref.First(n)-Sum_ses.First.SR(n));
    Sum_ses.First.p_Target(n)=DAG_fisher_exakt_scalars(Sum_ses.First.TL(n),Sum_ref.First(n)-Sum_ses.First.TL(n),Sum_ses.First.TR(n),Sum_ref.First(n)-Sum_ses.First.TR(n));
    
    Sum_norm.First.all(n)=Sum_ref.First(n)/Summary(n).All;
    %Sum_ses.First.all(n)=Sum_ref.First(n)/Summary(n).All;
    
    
    % Returns and not inspected targets
    
    %     Sum_ref.Ret(n)=sum([Summary(n).Total_inspections.LL, Summary(n).Total_inspections.LR, Summary(n).Total_inspections.RL, Summary(n).Total_inspections.RR]);
    %     Sum_ref.Nex(n)=sum([Summary(n).Total_inspections.LL, Summary(n).Total_inspections.LR, Summary(n).Total_inspections.RL, Summary(n).Total_inspections.RR]);
    
    %     Sum_ref.Ret(n)=sum([Summary(n).Returns.LL, Summary(n).Returns.LR, Summary(n).Returns.RL, Summary(n).Returns.RR]);
    %     Sum_ref.Nex(n)=sum([Summary(n).Notexplored.LL, Summary(n).Notexplored.LR, Summary(n).Notexplored.RL, Summary(n).Notexplored.RR]);
    %
    
    Sum_ref.Ret(n)= Summary(n).All;
    Sum_ref.Nex(n)= Summary(n).All;
    
    Sum_ses.Ret.all(n)=sum([Summary(n).Returns.LL, Summary(n).Returns.LR, Summary(n).Returns.RL, Summary(n).Returns.RR]);
    Sum_ses.Nex.all(n)=sum([Summary(n).Notexplored.LL, Summary(n).Notexplored.LR, Summary(n).Notexplored.RL, Summary(n).Notexplored.RR]);
    
    for k=1:numel(FN_sum)
        Sum_ses.Ret.(FN_sum{k})(n)=sum([Summary(n).Returns.(FN_in{k,1}), Summary(n).Returns.(FN_in{k,2})]);
        Sum_ses.Nex.(FN_sum{k})(n)=sum([Summary(n).Notexplored.(FN_in{k,1}), Summary(n).Notexplored.(FN_in{k,2})]);
        
        Sum_tot.Nex.(FN_sum{k})(n)= sum([Summary(n).Notexplored.LL, Summary(n).Notexplored.LR, Summary(n).Notexplored.RL, Summary(n).Notexplored.RR]);
        Sum_tot.Ret.(FN_sum{k})(n)= sum([Summary(n).Returns.LL, Summary(n).Returns.LR, Summary(n).Returns.RL, Summary(n).Returns.RR]);
        
        Sum_ses.TotIns.(FN_sum{k})(n)=sum([Summary(n).Total_inspections.(FN_in{k,1}), Summary(n).Total_inspections.(FN_in{k,2})]);
        
        Sum_norm.Ret.(FN_sum{k})(n)         = Sum_ses.Ret.(FN_sum{k})(n)/Sum_ref.Ret(n);
        Sum_norm.Nex.(FN_sum{k})(n)         = Sum_ses.Nex.(FN_sum{k})(n)/Sum_ref.Nex(n);
        Sum_stat.Ret.(FN_sum{k})(n)         = Sum_norm.Ret.(FN_sum{k})(n);
        Sum_stat.Nex.(FN_sum{k})(n)         = Sum_norm.Nex.(FN_sum{k})(n);
        %         Sum_norm.Ret.(FN_sum{k})(n)         = Sum_ses.Ret.(FN_sum{k})(n)/Sum_ses.TotIns.(FN_sum{k})(n);
        %         Sum_norm.Nex.(FN_sum{k})(n)         = Sum_ses.Nex.(FN_sum{k})(n)/Sum_ses.TotIns.(FN_sum{k})(n);
        
    end
    Sum_tot.Ret.all(n)  =Summary(n).All;
    Sum_tot.Nex.all(n)  =Summary(n).All;
    Sum_norm.Ret.all(n) =Sum_ref.Ret(n)/Summary(n).All;
    Sum_norm.Nex.all(n) =Sum_ref.Nex(n)/Summary(n).All;
    
    %     Sum_ses.Ret.p_Sample(n) =DAG_fisher_exakt_scalars(Sum_ses.Ret.SL(n),Sum_ses.TotIns.SL(n)-Sum_ses.Ret.SL(n),Sum_ses.Ret.SR(n),Sum_ses.TotIns.SR(n)-Sum_ses.Ret.SR(n));
    %     Sum_ses.Ret.p_Target(n) =DAG_fisher_exakt_scalars(Sum_ses.Ret.TL(n),Sum_ses.TotIns.TL(n)-Sum_ses.Ret.TL(n),Sum_ses.Ret.TR(n),Sum_ses.TotIns.TR(n)-Sum_ses.Ret.TR(n));
    %     Sum_ses.Nex.p_Sample(n) =DAG_fisher_exakt_scalars(Sum_ses.Nex.SL(n),Sum_ses.TotIns.SL(n)-Sum_ses.Nex.SL(n),Sum_ses.Nex.SR(n),Sum_ses.TotIns.SR(n)-Sum_ses.Nex.SR(n));
    %     Sum_ses.Nex.p_Target(n) =DAG_fisher_exakt_scalars(Sum_ses.Nex.TL(n),Sum_ses.TotIns.TL(n)-Sum_ses.Nex.TL(n),Sum_ses.Nex.TR(n),Sum_ses.TotIns.TR(n)-Sum_ses.Nex.TR(n));
    %
    
    Sum_ses.Ret.p_Sample(n) =DAG_fisher_exakt_scalars(Sum_ses.Ret.SL(n),Sum_ses.TotIns.SL(n)-Sum_ses.Ret.SL(n),Sum_ses.Ret.SR(n),Sum_ses.TotIns.SR(n)-Sum_ses.Ret.SR(n));
    Sum_ses.Ret.p_Target(n) =DAG_fisher_exakt_scalars(Sum_ses.Ret.TL(n),Sum_ses.TotIns.TL(n)-Sum_ses.Ret.TL(n),Sum_ses.Ret.TR(n),Sum_ses.TotIns.TR(n)-Sum_ses.Ret.TR(n));
    Sum_ses.Nex.p_Sample(n) =DAG_fisher_exakt_scalars(Sum_ses.Nex.SL(n),Sum_ses.TotIns.SL(n)-Sum_ses.Nex.SL(n),Sum_ses.Nex.SR(n),Sum_ses.TotIns.SR(n)-Sum_ses.Nex.SR(n));
    Sum_ses.Nex.p_Target(n) =DAG_fisher_exakt_scalars(Sum_ses.Nex.TL(n),Sum_ses.TotIns.TL(n)-Sum_ses.Nex.TL(n),Sum_ses.Nex.TR(n),Sum_ses.TotIns.TR(n)-Sum_ses.Nex.TR(n));
    
    
    % Timeouts
    Sum_ref.Timeout(n)=Summary(n).TimeOut.All;
    Sum_ses.Timeout.all(n)=Summary(n).TimeOut.All;
    for k=1:numel(FN_sum)
        Sum_ses.Timeout.(FN_sum{k})(n)          = sum([Summary(n).TimeOut.(FN_in{k,1}), Summary(n).TimeOut.(FN_in{k,2})]);
        Sum_tot.Timeout.(FN_sum{k})(n)          = Summary(n).TimeOut.All;
        Sum_norm.Timeout.(FN_sum{k})(n)         = Sum_ses.Timeout.(FN_sum{k})(n)/Sum_tot.Timeout.(FN_sum{k})(n);
        Sum_stat.Timeout.(FN_sum{k})(n)         = Sum_norm.Timeout.(FN_sum{k})(n);
    end
    
    Sum_tot.Timeout.all(n)=Summary(n).All;
    Sum_norm.Timeout.all(n)=Summary(n).TimeOut.All/Summary(n).All;
    
    %     Sum_ses.Timeout.p_Sample(n)=DAG_fisher_exakt_scalars(Sum_ses.Timeout.SL(n),Summary(n).TimeOut.All-Sum_ses.Timeout.SL(n),Sum_ses.Timeout.SR(n),Summary(n).TimeOut.All-Sum_ses.Timeout.SR(n));
    %     Sum_ses.Timeout.p_Target(n)=DAG_fisher_exakt_scalars(Sum_ses.Timeout.TL(n),Summary(n).TimeOut.All-Sum_ses.Timeout.TL(n),Sum_ses.Timeout.TR(n),Summary(n).TimeOut.All-Sum_ses.Timeout.TR(n));
    
    Sum_ses.Timeout.p_Sample(n)=DAG_fisher_exakt_scalars(Sum_ses.Timeout.SL(n),Summary(n).TimeOut.All-Sum_ses.Timeout.SL(n),Sum_ses.Timeout.SR(n),Summary(n).TimeOut.All-Sum_ses.Timeout.SR(n));
    Sum_ses.Timeout.p_Target(n)=DAG_fisher_exakt_scalars(Sum_ses.Timeout.TL(n),Summary(n).TimeOut.All-Sum_ses.Timeout.TL(n),Sum_ses.Timeout.TR(n),Summary(n).TimeOut.All-Sum_ses.Timeout.TR(n));
    
    
    
    % raw dwell times
    Sum_ref.Raw(n)=nanmean(Summary(n).explore.total);
    Sum_tot.Raw.all{n}=(Summary(n).explore.total);
    for k=1:numel(FN_sum)
        if k<=2
            Sum_ses.Raw.(FN_sum{k}){n}          = [Summary(n).explore.(FN_in{k,1})] + [Summary(n).explore.(FN_in{k,2})];
        else
            Sum_ses.Raw.(FN_sum{k}){n}          = [Summary(n).explore.(FN_in{k,1}), Summary(n).explore.(FN_in{k,2})];
        end
        Sum_tot.Raw.(FN_sum{k}){n}          = Summary(n).explore.total;
        Sum_norm.Raw.(FN_sum{k})(n)         = nanmean(Sum_ses.Raw.(FN_sum{k}){n});
        Sum_stat.Raw.(FN_sum{k})(n)         = nanmean(Sum_ses.Raw.(FN_sum{k}){n})/nanmean(Sum_tot.Raw.(FN_sum{k}){n});
    end
    
    [~,Sum_ses.Raw.p_Sample(n)]=ttest2([Summary(n).explore.LL, Summary(n).explore.LR],[Summary(n).explore.RL, Summary(n).explore.RR]);
    [~,Sum_ses.Raw.p_Target(n)]=ttest2([Summary(n).explore.LL, Summary(n).explore.RL],[Summary(n).explore.LR, Summary(n).explore.RR]);
    Sum_norm.Raw.all(n)=nanmean(Summary(n).explore.total); %why would you take median instead??
    Sum_ses.Raw.all{n}=Summary(n).explore.total;
    
    
    % Reaction times
    Sum_ref.RTs(n)=nanmean([Summary(n).RTs.LL, Summary(n).RTs.LR, Summary(n).RTs.RL, Summary(n).RTs.RR]);
    Sum_tot.RTs.all{n}=nanmean([Summary(n).RTs.LL, Summary(n).RTs.LR, Summary(n).RTs.RL, Summary(n).RTs.RR]);
    for k=1:numel(FN_sum)
        %Sum_ses.RTs.(FN_sum{k})(n)=nanmean([Summary(n).RTs.(FN_in{k,1}), Summary(n).RTs.(FN_in{k,2})]);
        Sum_ses.RTs.(FN_sum{k}){n}=[Summary(n).RTs.(FN_in{k,1}), Summary(n).RTs.(FN_in{k,2})];
        
        
        Sum_tot.RTs.(FN_sum{k}){n}=[Summary(n).RTs.LL, Summary(n).RTs.LR, Summary(n).RTs.RL, Summary(n).RTs.RR];
        %Sum_norm.RTs.(FN_sum{k})(n)=Sum_ses.RTs.(FN_sum{k})(n) - Sum_ref.RTs(n);
        
        Sum_stat.RTs.(FN_sum{k})(n)=nanmean(Sum_ses.RTs.(FN_sum{k}){n}) - nanmean(Sum_ref.RTs(n));
        Sum_norm.RTs.(FN_sum{k})(n)=nanmean(Sum_ses.RTs.(FN_sum{k}){n});
    end
    [~,Sum_ses.RTs.p_Sample(n)]=ttest2([Summary(n).RTs.LL, Summary(n).RTs.LR],[Summary(n).RTs.RL, Summary(n).RTs.RR]);
    [~,Sum_ses.RTs.p_Target(n)]=ttest2([Summary(n).RTs.LL, Summary(n).RTs.RL],[Summary(n).RTs.LR, Summary(n).RTs.RR]);
    Sum_norm.RTs.all(n)=Sum_ref.RTs(n);
    Sum_ses.RTs.all{n}=Sum_ref.RTs(n);
    
    % Dwell times
    
    Sum_ref.DTs(n)=nanmean([Summary(n).inspection_time_per_trial.LL', Summary(n).inspection_time_per_trial.LR',...
        Summary(n).inspection_time_per_trial.RL', Summary(n).inspection_time_per_trial.RR']);
    Sum_tot.DTs.all{n}=[Summary(n).inspection_time_per_trial.LL', Summary(n).inspection_time_per_trial.LR',...
        Summary(n).inspection_time_per_trial.RL', Summary(n).inspection_time_per_trial.RR'];
    for k=1:numel(FN_sum)
        %Sum_ses.DTs.(FN_sum{k})(n) =nanmean([Summary(n).inspection_time_per_trial.(FN_in{k,1})', Summary(n).inspection_time_per_trial.(FN_in{k,2})']);
        Sum_ses.DTs.(FN_sum{k}){n} =[Summary(n).inspection_time_per_trial.(FN_in{k,1})', Summary(n).inspection_time_per_trial.(FN_in{k,2})'];
        
        Sum_stat.DTs.(FN_sum{k})(n)= nanmean(Sum_ses.DTs.(FN_sum{k}){n}) - nanmean(Sum_ref.DTs(n));
        
        Sum_norm.DTs.(FN_sum{k})(n)= nanmean(double(Sum_ses.DTs.(FN_sum{k}){n}));  %% WHY IS THIS NOT DOUBLE???
        Sum_tot.DTs.(FN_sum{k}){n}= [Summary(n).inspection_time_per_trial.LL', Summary(n).inspection_time_per_trial.LR', Summary(n).inspection_time_per_trial.RL', Summary(n).inspection_time_per_trial.RR'];
    end
    Sum_norm.DTs.all(n)=double(Sum_ref.DTs(n)); %% WHY IS THIS NOT DOUBLE???
    Sum_ses.DTs.all{n}=Sum_ref.DTs(n);
    [h,Sum_ses.DTs.p_Sample(n)]=ttest2([Summary(n).inspection_time_per_trial.LL', Summary(n).inspection_time_per_trial.LR'],[Summary(n).inspection_time_per_trial.RL', Summary(n).inspection_time_per_trial.RR']);
    [h,Sum_ses.DTs.p_Target(n)]=ttest2([Summary(n).inspection_time_per_trial.LL', Summary(n).inspection_time_per_trial.RL'],[Summary(n).inspection_time_per_trial.LR', Summary(n).inspection_time_per_trial.RR']);
    
end



%% significance across sessions
fieldnames_perception_intention={'Hit','RTs','Raw','DTs','Timeout','First','Ret','Nex'};
for FN=fieldnames_perception_intention
    [h,Sum_norm.(FN{:}).p_Sample_across_Sessions] = ttest2(Sum_stat.(FN{:}).SL,Sum_stat.(FN{:}).SR);
    [h,Sum_norm.(FN{:}).p_Target_across_Sessions] = ttest2(Sum_stat.(FN{:}).TL,Sum_stat.(FN{:}).TR);
    for n=1:numel(Summary)
        Sum_norm.(FN{:}).diff_Sample_LR(n)=Sum_stat.(FN{:}).SL(n)-Sum_stat.(FN{:}).SR(n);
        Sum_norm.(FN{:}).diff_Target_LR(n)=Sum_stat.(FN{:}).TL(n)-Sum_stat.(FN{:}).TR(n);
    end
end


%% figure 1: Mean and STD of Normalized Values per Session (or subject), significance across sessions
plot_specs.figure_handle=figure('units','normalized','outerposition',[0 0 1 1]);

subplots_perception_intention=[1:8];
titles_perception_intention={'Success rate','Reaction time','Exploration time per trial','Dwell time per target','Timeouts','First response','Returns','Not explored one side'};
ylables_perception_intention={'Success rate','Reaction time [ms]','Exploration time [ms]','Dwell time [ms]','Fraction of Timeouts','Fraction of First Responses','Fraction of Returns','Fraction'};
colors_Perception_intention=[0 0 1; 0 1 0; 1 0 0; 1 0 1; 0.5 0.5 0.5];
FN={'SL','SR','TL','TR','all'};

for k=1:numel(fieldnames_perception_intention)
    subplot(3,4,subplots_perception_intention(k))
    xlim([0 8]);
    if strcmp(fieldnames_perception_intention{k},'First')
        set(gca,'xTick',[0:8],'xTickLabel',{'' 'Toward' 'Avert' '' 'L' 'R' '' 'All' ''},'fontsize',10);
    else
        set(gca,'xTick',[0:8],'xTickLabel',{'' 'L' 'R' '' 'L' 'R' '' 'All' ''},'fontsize',10);
    end
    hold on
    for m=1:numel(FN)
        bar(m+ceil(m/2)-1,nanmean(Sum_norm.(fieldnames_perception_intention{k}).(FN{m})),'facecolor',colors_Perception_intention(m,:));
        errorbar(m+ceil(m/2)-1,nanmean(Sum_norm.(fieldnames_perception_intention{k}).(FN{m})),nanstd(Sum_norm.(fieldnames_perception_intention{k}).(FN{m})),'k');
        text(m+ceil(m/2)-1,0,[num2str(round(nanmean(Sum_norm.(fieldnames_perception_intention{k}).(FN{m}))*100)/100) '+/-' num2str(round(nanstd(Sum_norm.(fieldnames_perception_intention{k}).(FN{m}))*100)/100)],'Rotation',90);
    end
    % significance bars (!)
    groups={[1,2],[4,5]};
    stats=[Sum_norm.(fieldnames_perception_intention{k}).p_Sample_across_Sessions, Sum_norm.(fieldnames_perception_intention{k}).p_Target_across_Sessions];
    nsort=0;
    H = sigstar(groups,stats,nsort);
    
    old_lim=get(gca,'ylim');
    set(gca,'ylim',[0 old_lim(2)]);
    %     text(0.8,old_lim(2)*0.9,'Sample');
    %     text(3.8,old_lim(2)*0.9,'Target');
    %ylim([0 old_lim(2)]);
    %xlabels(FN);
    title([titles_perception_intention{k} ' (ref = ' num2str(round(nanmean(Sum_ref.(fieldnames_perception_intention{k})))) ')'],'fontsize',16)
    ylabel(ylables_perception_intention{k},'fontsize',12);
    xlabel('Sample                      Target','fontsize',12)
end

fieldnames_confusion={'Hit','DTs'};
titles_confusion={'Confusions','Dwell time for difference in shape'};

subplots_object_related=[9,10];
%OBJ_FN
for k=1:numel(fieldnames_confusion)
    subplot(3,4,subplots_object_related(k))
    xlim([0 4]);
    set(gca,'xTick',[0:6],'xTickLabel',{'' 'L' 'R' 'LR'},'fontsize',10);
    hold on
    for m=1:numel(OBJ_FN)
        bar(m,nanmean(Sum_norm.(fieldnames_confusion{k}).(OBJ_FN{m})),'facecolor',colors_Perception_intention(m,:));
        errorbar(m,nanmean(Sum_norm.(fieldnames_confusion{k}).(OBJ_FN{m})),nanstd(Sum_norm.(fieldnames_confusion{k}).(OBJ_FN{m})),'k');
    end
    % significance bars (!)
    groups={[1,2],[4,5]};
    stats=[Sum_norm.(fieldnames_confusion{k}).p_Sample_across_Sessions, Sum_norm.(fieldnames_confusion{k}).p_Target_across_Sessions];
    nsort=0;
    H = sigstar(groups,stats,nsort);
    
    old_lim=get(gca,'ylim');
    %maxlim=max([old_lim(1)*-1 old_lim(2)]);
    %ylim([0 old_lim(2)]);
    %xlabels(FN);
    title(titles_confusion{k},'fontsize',16)
end


plot_specs.title='Group Summary';
plot_specs.figure_name_suffix='Group Summary 1';
title_and_save(plot_specs)

%% figure 2: N significant sessions
plot_specs.figure_handle=figure('units','normalized','outerposition',[0 0 1 1]);
Difference_fieldnames={'diff_Sample_LR','diff_Target_LR'};
P_fieldnames={'p_Sample','p_Target'};
subplots_perception_intention_significant_sessions=[1:8];

for k=1:numel(fieldnames_perception_intention)
    
    subplot(3,4,subplots_perception_intention_significant_sessions(k))
    xlim([0 3]);
    set(gca,'xTick',[0:3],'xTickLabel',{'' 'Sample' 'Target' ''},'fontsize',10);
    hold on
    col_idx=1;
    for m=1:numel(P_fieldnames)
        for s=[1,-1]
            idx_tot=sign([Sum_norm.(fieldnames_perception_intention{k}).(Difference_fieldnames{m})])==s;
            idx_sig=[Sum_ses.(fieldnames_perception_intention{k}).(P_fieldnames{m})]<= 0.05 & idx_tot;
            n_tot=s*sum(idx_tot);
            n_sig=s*sum(idx_sig);
            bar(m,n_tot,'edgecolor',colors_Perception_intention(col_idx,:),'facecolor','none');
            bar(m,n_sig,'facecolor',colors_Perception_intention(col_idx,:));
            col_idx=col_idx+1;
            tot_subj_labels=plot_specs.individual_labels(idx_tot & ~idx_sig);
            sig_subj_labels=plot_specs.individual_labels(idx_sig);
            for t1=1:n_sig*s
                text(m,s*(t1 - 0.5),sig_subj_labels{t1});
            end
            for t2=1:(n_tot-n_sig)*s
                text(m,s*(n_sig*s+t2-0.5),tot_subj_labels{t2});
            end
        end
    end
    
    
    old_lim=get(gca,'ylim');
    if strcmp(fieldnames_perception_intention{k},'First')
        text(0.3,old_lim(2)*0.9,'More Toward                    More Left');
        text(0.3,old_lim(1)*0.9,'More Avert                     More Right');
    else
        %         text(0.3,old_lim(2)*0.9,'More Left');
        %         text(0.3,old_lim(1)*0.9,'More Right');
    end
    %ylim([0 old_lim(2)]);
    %xlabels(FN);
    title(titles_perception_intention{k},'fontsize',16)
    ylabel('More Right     More Left','fontsize',12);
end

plot_specs.title='Group Summary';
plot_specs.figure_name_suffix='Group Summary 2';
title_and_save(plot_specs)
end


function title_and_save(plot_specs)
mtit(plot_specs.figure_handle,plot_specs.title, 'fontsize', plot_specs.main_title_font, 'xoff', -0.05, 'yoff', 0.03, 'color',[0 0 0]);
stampit;
if plot_specs.export_pdf
    wanted_size=[50 30];
    set(plot_specs.figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
    export_fig([plot_specs.path plot_specs.figure_name '_' plot_specs.figure_name_suffix], '-pdf','-transparent')
    close gcf
end
end


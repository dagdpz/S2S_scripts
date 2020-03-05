function S2S_group_comparison(Sum_ses,Sum_norm,Sum_tot,stats_test,comparison,subjects)
FN_C={'SL','SR','TL','TR','all'};
%FN_C={'SL','SR','diff_Sample_LR','TL','TR','diff_Target_LR','all'};
FN_CD={'SL','SR','diff_Sample_LR','TL','TR','diff_Target_LR','all'};
FN_D={'diff_Sample_LR','diff_Target_LR'};
FN_P={'Hit','RTs','Raw','DTs','Timeout','First','Ret','Nex'};
titles={'Success rate','Reaction time','Exploration per trial','Dwell time per target','Timeout','First response','Returns','Not explored one side'};
ylables={'Success rate','Reaction time [ms]','Exploration time [ms]','Dwell time [ms]','Fraction of trials','Fraction of trials','Fraction of trials','Fraction of trials'};
Colors=[0 0 1; 0 1 0; 1 0 0; 1 0 1; 0.3 0.3 0.3; 0 0.5 1; 0 1 0.5; 1 0.5 0; 1 0.5 1; 0.6 0.6 0.6];
ColorsD=[0 0.4 0.4; 0.4 0 0.4; 0 0.8 0.8; 0.8 0 0.8];
%Colors=[0 0 1; 0 1 0; 0 0.3 0.3; 1 0 0; 1 0 1; 0.3 0 0.3; 0.3 0.3 0.3; 0 0.5 1; 0 1 0.5; 0 0.5 0.5; 1 0.5 0; 1 0.5 1; 0.5 0 0.5; 0.6 0.6 0.6];
Parameters_ttest={'RTs','Raw','DTs'};
Parameters_fexact={'Timeout','Hit','First','Ret','Nex'};
if size(Sum_ses,1) >1
    %% this is for TMS
    % p this is unfortunately g in that condition
    for p=1:size(Sum_ses,1) 
        for FN=FN_C
            for k=1:numel(FN_P)
                means.(FN_P{k})(p).(FN{:})     = nanmean(Sum_norm(p,2).(FN_P{k}).(FN{:})-Sum_norm(p,1).(FN_P{k}).(FN{:}));
                sems.(FN_P{k})(p).(FN{:})      = sterr(Sum_norm(p,2).(FN_P{k}).(FN{:})-Sum_norm(p,1).(FN_P{k}).(FN{:}));
            end
        end
    end
else
    for g=1:size(Sum_ses,2)
        for FN=FN_CD
            for k=1:numel(FN_P)
                means.(FN_P{k})(g).(FN{:})     = nanmean(Sum_norm(g).(FN_P{k}).(FN{:}));
                sems.(FN_P{k})(g).(FN{:})      = sterr(Sum_norm(g).(FN_P{k}).(FN{:}));
            end
        end
    end
    Dif=Sum_ses;
    Tot=Sum_tot;
end
for FN=FN_CD
    for k=1:numel(FN_P)
        switch stats_test
            case 'paired_ttest'
                [~,p_across_subjects.(FN_P{k}).(FN{:})]=ttest([Sum_norm(1).(FN_P{k}).(FN{:})],[Sum_norm(2).(FN_P{k}).(FN{:})]);
            case 'unpaired_ttest'
                [~,p_across_subjects.(FN_P{k}).(FN{:})]=ttest2([Sum_norm(1).(FN_P{k}).(FN{:})],[Sum_norm(2).(FN_P{k}).(FN{:})]);
            case 'paired_nonparametric'
                p_across_subjects.(FN_P{k}).(FN{:})=signrank([Sum_norm(1).(FN_P{k}).(FN{:})],[Sum_norm(2).(FN_P{k}).(FN{:})]);
            case 'unpaired_nonparametric'
                p_across_subjects.(FN_P{k}).(FN{:})=ranksum([Sum_norm(1).(FN_P{k}).(FN{:})],[Sum_norm(2).(FN_P{k}).(FN{:})]);
        end
        if numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
            for s=1:numel(Sum_ses(1).Hit.SL)
                
                %% Disregard, only relevant for TMS
%                 if size(Sum_ses,1) >1
%                     if ismember(FN_P{k},Parameters_ttest)
%                         Dif(1).(FN_P{k}).(FN{:}){s}=Sum_ses(2,1).(FN_P{k}).(FN{:}){s}-nanmean([Sum_ses(1,1).(FN_P{k}).(FN{:}){s}]);
%                         Dif(2).(FN_P{k}).(FN{:}){s}=Sum_ses(2,2).(FN_P{k}).(FN{:}){s}-nanmean([Sum_ses(1,2).(FN_P{k}).(FN{:}){s}]);
%                     else
%                         Dif(1).(FN_P{k}).(FN{:})(s)=Sum_ses(2,1).(FN_P{k}).(FN{:})(s)-Sum_ses(1,1).(FN_P{k}).(FN{:})(s);
%                         Dif(2).(FN_P{k}).(FN{:})(s)=Sum_ses(2,2).(FN_P{k}).(FN{:})(s)-Sum_ses(1,2).(FN_P{k}).(FN{:})(s);
%                         Tot(1).(FN_P{k}).(FN{:})(s)=Sum_tot(2,1).(FN_P{k}).(FN{:})(s)+Sum_tot(1,1).(FN_P{k}).(FN{:})(s);
%                         Tot(2).(FN_P{k}).(FN{:})(s)=Sum_tot(2,2).(FN_P{k}).(FN{:})(s)+Sum_tot(1,2).(FN_P{k}).(FN{:})(s);
%                     end
%                 end
                if ismember(FN,FN_C)
                %% Within subject/session comparison
                if ismember(FN_P{k},Parameters_ttest)
                    p_within_subjects.(FN_P{k}).(FN{:})(s)=ttest2([Dif(1).(FN_P{k}).(FN{:}){s}],[Dif(2).(FN_P{k}).(FN{:}){s}]);
                else
                    p_within_subjects.(FN_P{k}).(FN{:})(s)=DAG_fisher_exakt_scalars(Dif(1).(FN_P{k}).(FN{:})(s),Tot(1).(FN_P{k}).(FN{:})(s)-Dif(1).(FN_P{k}).(FN{:})(s),...
                        Dif(2).(FN_P{k}).(FN{:})(s),Tot(2).(FN_P{k}).(FN{:})(s)-Dif(2).(FN_P{k}).(FN{:})(s));
                end
                end
            end
        end
    end
    
end




%setting xticklabels
xitcklabels={'' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' ''};
grouplabels={'SHAM','TMS'};
switch comparison
    case 'Inactivation'
        xitcklabels={'' 'C' 'I' '' 'C' 'I' '' 'C' 'I' '' 'C' 'I' '' 'C' 'I' ''};
        grouplabels={'Control',' Inactivation'};
    case '3 targets Cornelius 1 cue vs 2 cues'
        xitcklabels={'' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' ''};
        grouplabels={'1 cue',' 2 cues'};
    case 'SHAM vs TMS young'
        xitcklabels={'' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' ''};
        grouplabels={'SHAM','TMS'};
    case 'preSHAM vs preTMS young'
        xitcklabels={'' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' ''};
        grouplabels={'preSHAM','preTMS'};
    case 'preTMS vs postTMS young'
        xitcklabels={'' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' ''};
        grouplabels={'preTMS','postTMS'};
    case 'preSHAM vs postSHAM young'
        xitcklabels={'' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' '' 'Pre' 'Post' ''};
        grouplabels={'preSHAM','postSHAM'};
    case 'postSHAM vs postTMS young'
        xitcklabels={'' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' '' 'S' 'T' ''};
        grouplabels={'postSHAM','postTMS'};
    case '7 targets'
        xitcklabels={''};
        grouplabels={'',''};
    case '3 targets 1 cue Y vs O'
        xitcklabels={'' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' ''};
    case '3 targets 1 cue Y vs Cornelius'
        xitcklabels={'' 'Y' 'M' '' 'Y' 'M' '' 'Y' 'M' '' 'Y' 'M' '' 'Y' 'M' ''};
        grouplabels={'Human','Monkey'};
    case '3 targets 1 cue ideal vs perceptual neglect'
        xitcklabels={'' 'I' 'P' '' 'I' 'P' '' 'I' 'P' '' 'I' 'P' '' 'I' 'P' ''};
        grouplabels={'Control','Neglect'};
    case '3 targets 1 cue ideal vs intentional neglect'
        xitcklabels={'' 'I' 'M' '' 'I' 'M' '' 'I' 'M' '' 'I' 'M' '' 'I' 'M' ''};
        grouplabels={'Control','Neglect'};
    case '3 targets 1 cue Y vs perceptual neglect'
        xitcklabels={'' 'Y' 'P' '' 'Y' 'P' '' 'Y' 'P' '' 'Y' 'P' '' 'Y' 'P' ''};
        grouplabels={'Control','Neglect'};
    case '3 targets 2 cues Y vs O'
        xitcklabels={'' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' ''};
        grouplabels={'Young','Old'};
    case '3 targets naive Y vs O'
        xitcklabels={'' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' '' 'Y' 'O' ''};
        grouplabels={'Young','Old'};
    case '3 targets Y 1 vs 2 cues'
        xitcklabels={'' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' ''};
        grouplabels={'1 cue','2 cues'};
    case '3 targets O 1 vs 2 cues'
        xitcklabels={'' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' '' '1' '2' ''};
        grouplabels={'1 cue','2 cues'};
    case 'GE pre vs post mem_per'
        xitcklabels={''};
        grouplabels={'',''};
    case 'Patient'
        xitcklabels={''};
        grouplabels={'',''};
end


% Mean and STD of Normalized Values per Session (or subject), significance
% across sessions (or subjects)
summary_figure1=figure('units','normalized','outerposition',[0 0 1 1]);
for k=1:numel(FN_P)
    subplot(3,4,k)
    %xlim([0 15]);
    xlim([0 15]);
    set(gca,'xTick',[0:15],'xTickLabel',xitcklabels,'fontsize',10);
    hold on
    for m=1:numel(FN_C)
        x=3*m-2;
        bar(x,means.(FN_P{k})(1).(FN_C{m}),'facecolor',Colors(m,:));
        errorbar(x,means.(FN_P{k})(1).(FN_C{m}),sems.(FN_P{k})(1).(FN_C{m}),'k');
        text(x,0,[num2str(round(means.(FN_P{k})(1).(FN_C{m})*100)/100) '+/-' num2str(round(sems.(FN_P{k})(1).(FN_C{m})*100)/100)],'Rotation',90);
        bar(x+1,means.(FN_P{k})(2).(FN_C{m}),'facecolor',Colors(m+5,:));
        errorbar(x+1,means.(FN_P{k})(2).(FN_C{m}),sems.(FN_P{k})(2).(FN_C{m}),'k');
        text(x+1,0,[num2str(round(means.(FN_P{k})(2).(FN_C{m})*100)/100) '+/-' num2str(round(sems.(FN_P{k})(2).(FN_C{m})*100)/100)],'Rotation',90);
    end
    % significance bars (!)
    groups={[1,2],[4,5],[7,8],[10,11],[13,14]};
    stats=[p_across_subjects.(FN_P{k}).(FN_C{1}),...
        p_across_subjects.(FN_P{k}).(FN_C{2}),...
        p_across_subjects.(FN_P{k}).(FN_C{3}),...
        p_across_subjects.(FN_P{k}).(FN_C{4}),...
        p_across_subjects.(FN_P{k}).(FN_C{5})];
    nsort=0;
    H = sigstar(groups,stats,nsort);
    title([titles{k}],'fontsize',16)
    ylabel(ylables{k},'fontsize',12);
    if strcmp(FN_P{k},'First')
        xlabel('To S    Avert S   L Side     R Side     All   ','fontsize',12)
    elseif strcmp(FN_P{k},'Hit')
        xlabel('SL       SR       ML         MR         All   ','fontsize',12)
    else
        xlabel('SL       SR       L Side     R Side     All   ','fontsize',12)
    end
end

% Mean and STD of Normalized Differences (Left - right) per Session (or subject), significance
% across sessions (or subjects)
summary_figure2=figure('units','normalized','outerposition',[0 0 1 1]);
for k=1:numel(FN_P)
    subplot(3,4,k)
    %xlim([0 15]);
    xlim([0 6]);
    set(gca,'xTick',[0:15],'xTickLabel',xitcklabels,'fontsize',10);
    hold on
    for m=1:numel(FN_D)
        x=3*m-2;
        bar(x,means.(FN_P{k})(1).(FN_D{m}),'facecolor',ColorsD(m,:));
        errorbar(x,means.(FN_P{k})(1).(FN_D{m}),sems.(FN_P{k})(1).(FN_D{m}),'k');
        ypos=double(min([0,means.(FN_P{k})(1).(FN_D{m})]));
        text(x,ypos,[num2str(round(means.(FN_P{k})(1).(FN_D{m})*100)/100) '+/-' num2str(round(sems.(FN_P{k})(1).(FN_D{m})*100)/100)],'Rotation',90);
        bar(x+1,means.(FN_P{k})(2).(FN_D{m}),'facecolor',ColorsD(m+2,:));
        errorbar(x+1,means.(FN_P{k})(2).(FN_D{m}),sems.(FN_P{k})(2).(FN_D{m}),'k');
        ypos=double(min([0,means.(FN_P{k})(2).(FN_D{m})]));
        text(x+1,ypos,[num2str(round(means.(FN_P{k})(2).(FN_D{m})*100)/100) '+/-' num2str(round(sems.(FN_P{k})(2).(FN_D{m})*100)/100)],'Rotation',90);
    end
    % significance bars (!)
    groups={[1,2],[4,5],};
    stats=[p_across_subjects.(FN_P{k}).(FN_D{1}),...
        p_across_subjects.(FN_P{k}).(FN_D{2})];
    nsort=0;
    H = sigstar(groups,stats,nsort);
    title([titles{k}],'fontsize',16)
    ylabel(ylables{k},'fontsize',12);
    if strcmp(FN_P{k},'First')
        xlabel('To-Avert S   Side(L-R)','fontsize',12)
    elseif strcmp(FN_P{k},'Hit')
        xlabel('Sample (L-R)   Match(L-R)   ','fontsize',12)
    else
        xlabel('Sample (L-R)   Side(L-R)   ','fontsize',12)
    end
end



% MNumber of significant Sessions (or subjects)
if numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
    summary_figure3=figure('units','normalized','outerposition',[0 0 1 1]);
    for k=1:numel(FN_P)
        subplot(3,4,k)
        xlim([0 6]);
        set(gca,'xTick',[0:6],'xTickLabel',{'' 'SL' 'SR' 'TL' 'TR' 'ALL' ''},'fontsize',10);
        hold on
        
        col_idx=1;
        for s=[-1,1]
            for m=1:numel(FN_C)
                FN=FN_C{m};
                idx_tot=sign([Sum_norm(2).(FN_P{k}).(FN)]-[Sum_norm(1).(FN_P{k}).(FN)])==s;
                idx_sig=p_within_subjects.(FN_P{k}).(FN) <= 0.05 & sign([Sum_norm(2).(FN_P{k}).(FN)]-[Sum_norm(1).(FN_P{k}).(FN)])==s;
                n_tot=s*sum(idx_tot);
                n_sig=s*sum(idx_sig);
                bar(m,n_tot,'edgecolor',Colors(col_idx,:),'facecolor','none');
                bar(m,n_sig,'facecolor',Colors(col_idx,:));
                col_idx=col_idx+1;
                
                tot_subj_labels=subjects(idx_tot & ~idx_sig);
                sig_subj_labels=subjects(idx_sig);
                for t1=1:n_sig*s
                    text(m,s*(t1 - 0.5),sig_subj_labels{t1});
                end
                for t2=1:(n_tot-n_sig)*s
                    text(m,s*(n_sig*s+t2-0.5),tot_subj_labels{t2});
                end
            end
        end
        old_lim=get(gca,'ylim');
        text(2,old_lim(2)*0.9,['More ' grouplabels{2}]);
        text(2,old_lim(1)*0.9,['More ' grouplabels{1}]);
        title(titles{k},'fontsize',16)
        ylabel('N','fontsize',12);
    end
end


a=1;
end







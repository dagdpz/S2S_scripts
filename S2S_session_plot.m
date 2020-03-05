function S2S_session_plot(M2S_output)
global GLO
GLO.linewidth=3;
GLO.main_title_font=15;

GLO.title_font=16;
GLO.ylabel_font=14;
GLO.xlabel_font=14;
GLO.tick_font=5;
GLO.text_font=8;


Bins    =M2S_output.Bins;
Sel     =M2S_output.Sel;
Tot     =M2S_output.Tot;
RTs     =M2S_output.RTs;

inspection_time         =M2S_output.inspection_time;
possible_convexities    =M2S_output.possible_convexities;
possible_positions      =M2S_output.possible_positions;
possible_sides          =M2S_output.possible_sides;
returns_to_revealed     =M2S_output.returns_to_revealed;
% revealed_tar_pos_LR     =M2S_output.revealed_tar_pos_LR;
revealed_tar_pos        =M2S_output.revealed_tar_pos;
Summary                 =M2S_output.Summary;
Disp                    =M2S_output.Disp;
Pattern                 =M2S_output.Pattern;
all_n_positions         =1:numel(possible_positions);

shape_lables            =[cellstr(num2str(possible_convexities)),repmat({' '},size(possible_sides)),possible_sides];


GLO.shapecolors                 =copper(numel(possible_convexities));
GLO.samplecolors                =[0 0 1; 0 1 0];
%GLO.space_colors_LRN           =[0,0.5,0.8;0,0.8,0.5;0.5,0.5,0.5];
GLO.space_colors_LRN            =[0,0,1;0,1,0;0.5,0.5,0.5];
GLO.space_colors_RT             =[0,0,1;0,1,0;1,0,1;1,1,0;];

colors_used                     =jet(numel(possible_positions));
GLO.space_colors_allpositions   =[colors_used;0.5,0.5,0.5];
GLO.search_intervals            =jet(numel(inspection_time.mean_search_to_reveal));
GLO.fix_y = Disp.fix_y;



print_out=([ Disp.sessions_runs ' ' Disp.S_or_R ': ' Disp.n_targets_string ' targets, ' Disp.n_cues ' at ' Disp.samplepos ' deg, ' Disp.masked_or_not]);

%% Figure 1

%alpha=0.7;
n_rows=min(numel(possible_convexities)+2,6);
n_columns=numel(possible_convexities)+2;

summary_figure1=figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'renderer','painters');

colormap('hot');
%colormap(gca,jet)

%% Confusion matrix
for k=1:numel(possible_convexities)
    for m=1:numel(possible_convexities)
        LL{k,m}=Sel.LL{k,m}/Tot.LL{k,m};
        LR{k,m}=Sel.LR{k,m}/Tot.LR{k,m};
        
        RL{k,m}=Sel.RL{k,m}/Tot.RL{k,m};
        RR{k,m}=Sel.RR{k,m}/Tot.RR{k,m};
        
        sum_cor{k,m}=Sel.LL{k,m}+Sel.LR{k,m}+Sel.RL{k,m}+Sel.RR{k,m};
        sum_tot{k,m}=Tot.all{k,m};
    end
end

for k=1:numel(possible_convexities)
    for m=1:numel(possible_convexities)
        subplot(n_rows,n_columns,(k)*n_columns+m+1)
        hold on
        
        bH(1)=barh(2,-(Sel.LL{k,m}/Tot.LL{k,m}),'Facecolor',GLO.samplecolors(1,:));
        bH(2)=barh(2, (Sel.LR{k,m}/Tot.LR{k,m}),'Facecolor',GLO.samplecolors(1,:));
        bH(4)=barh(1,-(Sel.RL{k,m}/Tot.RL{k,m}),'Facecolor',GLO.samplecolors(2,:));
        bH(2)=barh(1, (Sel.RR{k,m}/Tot.RR{k,m}),'Facecolor',GLO.samplecolors(2,:));
        
        text(-0.5,2,[num2str(Sel.LL{k,m}) '/' num2str(Tot.LL{k,m})]);
        text( 0.5,2,[num2str(Sel.LR{k,m}) '/' num2str(Tot.LR{k,m})]);
        text(-0.5,1,[num2str(Sel.RL{k,m}) '/' num2str(Tot.RL{k,m})]);
        text( 0.5,1,[num2str(Sel.RR{k,m}) '/' num2str(Tot.RR{k,m})]);
        
        if (Tot.LL{k,m})==0 || isnan(Sel.LL{k,m});barh(2,-1,'k'); end
        if (Tot.LR{k,m})==0 || isnan(Sel.LR{k,m});barh(2,1,'k');  end
        if (Tot.RL{k,m})==0 || isnan(Sel.RL{k,m});barh(1,-1,'k'); end
        if (Tot.RR{k,m})==0 || isnan(Sel.RR{k,m});barh(1,1,'k');  end
        title('Target left  Target right','fontsize',GLO.title_font)
        ylabel('Sample','fontsize',GLO.ylabel_font)
        Ylabels= {'R','L'};
        %set(gca,'ytick',[1:2],'ylim',[0 3],'xlim',[-1 1],'yticklabel',[Ylabels],'xtick',[-1:0.5:1],'xticklabel',[1 0.5 0 0.5 1]);
        set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[]);
    end
end

%% sample shapes
center=[0 0];
target_size=5;
for k=1:numel(possible_convexities)
    subplot(n_rows,n_columns,k+1)
    pointList=CalculateConvexPointList(center,target_size,possible_convexities(k),possible_sides{k});
    
    fill(pointList(:,1),pointList(:,2),[1,0,0]);
    text(-target_size*0.8, target_size*1.3 ,'Target','fontsize',GLO.title_font);
    text(-target_size*0.55, target_size*0,[shape_lables{k,:}],'fontsize',GLO.text_font+5);
    
    xlim([center(1)-2*target_size center(1)+2*target_size]);
    axis equal
    axis off
    
    subplot(n_rows,n_columns,k*n_columns+1,'Color',GLO.shapecolors(k,:))
    fill(pointList(:,1),pointList(:,2),[1 0 0]);
    text(-target_size*0.8,   target_size*1.3 ,'Sample','fontsize',GLO.title_font,'Color',GLO.shapecolors(k,:));
    text(-target_size*0.55, target_size*0,[shape_lables{k,:}],'fontsize',GLO.text_font+5);
    
    xlim([center(1)-2*target_size center(1)+2*target_size]);
    axis equal
    axis off
end

%% Confusion Matrix Summary
subplot(n_rows,n_columns,n_columns)
hold off
Confusions=cell2mat(Summary.confusion.Sel)./cell2mat(Summary.confusion.Tot);
imagesc(Confusions,[0 1]);
title(['Confusion Matrix'],'fontsize',GLO.title_font);
add_normalized_colorbar(100,'%');
axis equal

%% Summary Hitrates
subplot(n_rows,n_columns,1)
hold on

%bar(1,Summary.All/(Summary.All + Summary.TimeOut.All),'r');
bar(1,Summary.All_success/Summary.All_match_revealed_completed,'r');
bar(2,Summary.Succes_rate,'r');
bar(3,Summary.Succes_rate_cue_L,'Facecolor',GLO.samplecolors(1,:));
bar(4,Summary.Succes_rate_cue_R,'Facecolor',GLO.samplecolors(2,:));

To_display_succes={Summary.All_success,Summary.All_success,Summary.cue_L_success,Summary.cue_R_success};
To_display_all={Summary.All_match_revealed_completed,Summary.All,Summary.cue_L,Summary.cue_R};
for k=1:4
    text(k-0.3,0.32,num2str(To_display_succes{k}));
    text(k-0.3,0.2 ,'of');
    text(k-0.3,0.08,num2str(To_display_all{k}));
end

% text(m-0.5,0.9,[num2str(sum_cor{k,m}) '/' num2str(sum_tot{k,m})]);
set(gca,'xtick',[1:4],'xlim',[0 5],'xticklabel',{'M rev','total','sampL','sampR'})
title('Over all hitrate','fontsize',GLO.title_font)

%% TimeOuts

subplot(n_rows,n_columns,n_rows*(n_rows-1)+1)
hold on

FN={'RL','RR','LL','LR'};
CO={GLO.samplecolors(2,:),GLO.samplecolors(2,:),GLO.samplecolors(1,:),GLO.samplecolors(1,:)};
hold on

N=Summary.TimeOut.All;
for k=1:numel(FN)/2
    m=k*2-1;
    bbH(1)=barh(k,-Summary.TimeOut.(FN{m})/N,'Facecolor',CO{m});
    bbH(2)=barh(k, Summary.TimeOut.(FN{m+1})/N,'Facecolor',CO{m});
    
    text(-0.5,k,[num2str(Summary.TimeOut.(FN{m})) '/' num2str(N)]);
    text( 0.5,k,[num2str(Summary.TimeOut.(FN{m+1})) '/' num2str(N)]);
end
title('Match L |Timeouts| Match R','fontsize',GLO.title_font)
ylabel('Sample','fontsize',GLO.ylabel_font)
Ylabels= {'R','L'};%,'RC','LC','RU','LU'};
set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[]);

%% Fixation breaks
subplot(n_rows,n_columns,2*n_columns)
FN={'RL','RR','LL','LR'};
CO={GLO.samplecolors(2,:),GLO.samplecolors(2,:),GLO.samplecolors(1,:),GLO.samplecolors(1,:)};
hold on
for k=1:numel(FN)/2
    m=k*2-1;
    bbH(1)=barh(k,-Summary.Fix.(FN{m})  /Summary.Fixation_breaks,'Facecolor',CO{m});
    bbH(2)=barh(k, Summary.Fix.(FN{m+1})/Summary.Fixation_breaks,'Facecolor',CO{m});
    
    text(-0.5,k,num2str(Summary.Fix.(FN{m})));
    text( 0.5,k,num2str(Summary.Fix.(FN{m+1})));
    
    if Summary.Fixation_breaks==0; barh(k,-1,'k') ;barh(k,1, 'k'); end
end
title('Breaks left  Breaks right','fontsize',GLO.title_font)
ylabel('Sample','fontsize',GLO.ylabel_font)
Ylabels= {'R','L'};%,'RC','LC','RU','LU'};
set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[]);

%% Fixation break reaction times

subplot_assignments=[3*n_columns,4*n_columns];
ALL_RT_fieldnames={'LL','RL';'LR','RR'};
titles={'left','right'};
for k=1:2
    subplot(n_rows,n_columns,subplot_assignments(k))
    RT_fieldnames=ALL_RT_fieldnames(k,:);
    RT_lables={};
    hold on
    for m=1:numel(RT_fieldnames)
        current_structure=RTs.histograms_fb.(RT_fieldnames{m});
        current_median=RTs.median_fb.(RT_fieldnames{m});
        hfb(m)=plot(Bins.reaction_time,norm_hist(current_structure),'color',GLO.space_colors_RT(m,:),'linewidth',GLO.linewidth);
        line([current_median current_median],[0 1],'color',GLO.space_colors_RT(m,:),'linewidth',GLO.linewidth-1.5);
        RT_lables=[RT_lables; {[num2str(round(current_median)) ' ms']}];
    end
    legend(hfb,RT_lables,'Location','best');
    title(['Breaks ' titles{k}],'fontsize',GLO.title_font)
    ylim([0 1]);
    ylabel('Probability','fontsize',GLO.ylabel_font);
    xlabel('time after cue [ms]','fontsize',GLO.xlabel_font);
end

%% Errors
subplot(n_rows,n_columns,(n_rows-1)*n_columns)
FN={'RL','RR','LL','LR'};
CO={GLO.samplecolors(2,:),GLO.samplecolors(2,:),GLO.samplecolors(1,:),GLO.samplecolors(1,:)};
hold on
for k=1:numel(FN)/2
    m=k*2-1;
    %     bbH(1)=barh(k,-Summary.Err.(FN{m})  /Summary.TWr.(FN{m}) *numel(possible_convexities)/2,'Facecolor',CO{m});
    %     bbH(2)=barh(k, Summary.Err.(FN{m+1})/Summary.TWr.(FN{m+1})*numel(possible_convexities)/2,'Facecolor',CO{m});
    bbH(1)=barh(k,-Summary.Err.(FN{m})  /Summary.TWr.(FN{m}),'Facecolor',CO{m});
    bbH(2)=barh(k, Summary.Err.(FN{m+1})/Summary.TWr.(FN{m+1}),'Facecolor',CO{m});
    
    text(-0.5,k,num2str(Summary.Err.(FN{m})));
    text( 0.5,k,num2str(Summary.Err.(FN{m+1})));
end
title('Errors','fontsize',GLO.title_font)
ylabel('Sample','fontsize',GLO.ylabel_font)
Ylabels= {'R','L'};%,'RC','LC','RU','LU'};
set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[]);

%% Hitrates
subplot(n_rows,n_columns,n_rows*n_columns)
%FN={'LR','RR','LL','RL'};
FN={'RL','RR','LL','LR'};
CO={GLO.samplecolors(2,:),GLO.samplecolors(2,:),GLO.samplecolors(1,:),GLO.samplecolors(1,:)};
hold on
for k=1:numel(FN)/2
    m=k*2-1;
    bbH(1)=barh(k,-sum([Sel.(FN{m}){eye(numel(possible_convexities))==1}])  /sum([Tot.(FN{m}){eye(numel(possible_convexities))==1}]),'Facecolor',CO{m});
    bbH(2)=barh(k, sum([Sel.(FN{m+1}){eye(numel(possible_convexities))==1}])/sum([Tot.(FN{m+1}){eye(numel(possible_convexities))==1}]),'Facecolor',CO{m});
    
    text(-0.5,k,[num2str(sum([Sel.(FN{m}){eye(numel(possible_convexities))==1}])) '/' num2str(sum([Tot.(FN{m}){eye(numel(possible_convexities))==1}]))]);
    text( 0.5,k,[num2str(sum([Sel.(FN{m+1}){eye(numel(possible_convexities))==1}])) '/' num2str(sum([Tot.(FN{m+1}){eye(numel(possible_convexities))==1}]))]);
    
    if sum([Tot.(FN{m}){eye(numel(possible_convexities))==1}])==0;barh(k,-1,'k'); end
    if sum([Tot.(FN{m+1}){eye(numel(possible_convexities))==1}])==0;barh(k,1, 'k'); end
end
title('Hitrate left  Hitrate right','fontsize',GLO.title_font)
ylabel('Sample','fontsize',GLO.ylabel_font)
Ylabels= {'R','L'};%,'RC','LC','RU','LU'};
set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[]);
%set(gca,'ytick',[1:2],'ylim',[0 3],'yticklabel',[Ylabels],'xlim',[-1 1],'xtick',[-1:0.5:1],'xticklabel',[1 0.5 0 0.5 1]);


%% DT per target shape
for m=1:numel(possible_convexities)
    subplot(n_rows,n_columns,n_rows*(n_columns-1)+1+m)
    current_structure=inspection_time.hist_per_convexity;
    current_median=inspection_time.median_max_per_convexity;
    
    hold on
    Inspection_lables={};
    hi=NaN(1,numel(possible_convexities));
    for k=1:numel(possible_convexities)
        if m==k
            continue
        end
        hi(k)=plot(Bins.inspection_mean,norm_hist(current_structure{k,m}),'color',GLO.shapecolors(k,:),'linewidth',GLO.linewidth);
        line([current_median{k,m} current_median{k,m}],[0 1],'color',GLO.shapecolors(k,:),'linewidth',GLO.linewidth-1.5);
        Inspection_lables= [Inspection_lables; {[num2str(round(current_median{k,m})) ' ms']}];
    end
    hi(m)=[];
    
    title(['Dwell time ' [shape_lables{m,:}]],'fontsize',GLO.title_font)
    ylim([0 1]);
    ylabel('Probability','fontsize',GLO.ylabel_font);
    xlabel('time [s]','fontsize',GLO.xlabel_font);
    
    legend(hi,Inspection_lables,'Location','best');
end

%% Figure 1 end
% mtit(summary_figure1,print_out, 'fontsize', GLO.main_title_font, 'xoff', -0.05, 'yoff', 0.03, 'color',[0 0 0]);
% pdf_export(Disp,'confusion_matrix');

title_and_save(summary_figure1,print_out,'confusion_matrix',Disp);

%% Figure 2: Exploration

summary_figure2=figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'renderer','painters');
n_columns=5;
n_rows=5;

it_subplot_in_row=[2,3];
pie_chart_lables={'1st pos','2nd pos','3rd pos','4th pos','5th pos','6th pos','7th pos','8th pos','9th pos','10th pos',...
    '11th pos','12th pos','13th pos','14th pos','15th pos','16th pos','17th pos','18th pos','19th pos','20th pos',...
    '21th pos','22th pos','23th pos','24th pos','25th pos','26th pos','27th pos','28th pos','29th pos','30th pos',...
    '31th pos','32th pos','33th pos','34th pos','35th pos','36th pos','37th pos','38th pos','39th pos','40th pos',...
    '41th pos','42th pos','43th pos','44th pos','45th pos','46th pos','47th pos','48th pos','49th pos','50th pos','not continued'};
pie_chart_data=revealed_tar_pos.selected;

%% positions plot
subplot(n_rows,n_columns,[4,n_columns+5])
plot(0,Disp.fix_y,'+','MarkerSize',20,'Color','k');
for pos_idx=1:numel(possible_positions)
    hold on
    plot(real(possible_positions(pos_idx)),imag(possible_positions(pos_idx)),'o','MarkerSize',20,'MarkerFaceColor',GLO.space_colors_allpositions(pos_idx,:),'MarkerEdgeColor',GLO.space_colors_allpositions(pos_idx,:))
    %text(real(possible_positions(pos_idx))-max(real(possible_positions))/9,imag(possible_positions(pos_idx)),num2str(pos_idx),'fontsize',20)
    text(real(possible_positions(pos_idx))-max(real(possible_positions))/9,imag(possible_positions(pos_idx)),num2str(round(inspection_time.mean_spot_in_revealing_sequence(pos_idx)*10)/10),'fontsize',20)
end
title('Target Positions','fontsize',GLO.title_font)
axis equal

y_range_min=20;
y_range=max(max(imag(possible_positions))- min(imag(possible_positions)),y_range_min);
set(gca,'ylim',[min(imag(possible_positions))-y_range/5 max(imag(possible_positions))+y_range/5]);

x_range=max(real(possible_positions))- min(real(possible_positions));
set(gca,'xlim',[min(real(possible_positions))-x_range/5 max(real(possible_positions))+x_range/5]);

maximum_ratio=1;
first_target_nonnan{1}=revealed_tar_pos.nonnan{1}.*maximum_ratio;
first_target_data=revealed_tar_pos.selected;

%% pie chart(s)
%pie_chart_pos=[2:3,2+n_columns:3+n_columns,2+2*n_columns:3+2*n_columns];
pie_chart_pos=[3*n_columns-1:3*n_columns,4*n_columns-1:4*n_columns,5*n_columns-1:5*n_columns];
subplot(n_rows,n_columns,pie_chart_pos)
if Disp.plot_only_LR
    %   title('Exploration (more targets left)','fontsize',GLO.title_font);
else
    title('Exploration pattern','fontsize',GLO.title_font);
    pie_chart_plot(Disp.plot_only_LR,pie_chart_data,pie_chart_lables)
end


%% Heat map for % sucess rate per position
subplot(n_rows,n_columns,1)
Normalizing_time{1}=ones(size(possible_positions));
heat_map_plot(revealed_tar_pos.success_rate_per_target,Normalizing_time,possible_positions)
title('Success rate per position','fontsize',GLO.title_font);

%% Heat map for % first revealed
subplot(n_rows,n_columns,2)
heat_map_plot(first_target_data,first_target_nonnan,possible_positions)
title('1st revealed [% of present]','fontsize',GLO.title_font);

%% Complete heat map (!) for total dwell time
subplot(n_rows,n_columns,3)
Time_to_normalize=max(inspection_time.sum_per_position{1});
Normalizing_time{1}=repmat(Time_to_normalize,size(possible_positions));

%heat_map_plot(inspection_time.sum_per_position,Normalizing_time,possible_positions,round(inspection_time.mean_spot_in_revealing_sequence*10)/10);
heat_map_plot(inspection_time.sum_per_position,Normalizing_time,possible_positions);
title(['Total DT / position [% of ' num2str(round(Time_to_normalize*1000/Summary.All)) 'ms]'],'fontsize',GLO.title_font);

%% Probability N revealed plot
subplot(n_rows,n_columns,n_columns+1)
hold on
bar(Summary.total_revealed/sum(Summary.total_revealed), 'FaceColor','r');
bar(Summary.completed_revealed/sum(Summary.total_revealed), 'FaceColor',[1 0.5 0]);
bar(Summary.completed_revealed_including_match/sum(Summary.total_revealed), 'FaceColor',[1 1 0]);
bar(Summary.successful_revealed/sum(Summary.total_revealed), 'FaceColor',[0 1 0]);

%bar(Summary.successful_revealed/sum(Summary.total_revealed), 'FaceColor','g');
title('N revealed','fontsize',GLO.title_font)
set(gca,'xlim',[0.5 numel(Summary.total_revealed)+0.5],'ytick',[0:0.1:1],'xtick',[1:numel(possible_positions)]);%'ytick',[1:6],'ylim',[0 7]
legend({'Timeouts','M hidden', 'M neglected', 'M selected'},'location','Northeast');
xlabel('N revealed','fontsize',GLO.xlabel_font)
ylabel('Probability','fontsize',GLO.ylabel_font)

%% Hitrate per N revealed plot
subplot(n_rows,n_columns,n_columns+2)
hold on
%bar(Summary.total_revealed/sum(Summary.total_revealed), 'FaceColor','r');
bar(Summary.successful_revealed./Summary.total_revealed, 'FaceColor','g');
title('Success','fontsize',GLO.title_font)
set(gca,'xlim',[0.5 numel(Summary.successful_revealed)+0.5],'ytick',[0:0.2:1],'xtick',[1:numel(possible_positions)]);%'ytick',[1:6],'ylim',[0 7]
xlabel('N revealed','fontsize',GLO.xlabel_font)
ylabel('Hitrate','fontsize',GLO.ylabel_font)

%% Return histogram plot
subplot(n_rows,n_columns,n_columns+3)
title(['Return histograms'],'fontsize',GLO.title_font)
xlabel('Max times inspected [#]','fontsize',GLO.xlabel_font);
ylabel('N_t_r_i_a_l_s','fontsize',GLO.ylabel_font);
bar_min=min([-returns_to_revealed.hist.RL,returns_to_revealed.hist.RL,-returns_to_revealed.hist.RR,returns_to_revealed.hist.RR]);
bar_max=max([-returns_to_revealed.hist.LL,returns_to_revealed.hist.LL,-returns_to_revealed.hist.LR,returns_to_revealed.hist.LR]);
set(gca,'ylim',[bar_min*1.5, bar_max*1.5],'xlim',[Bins.n_return(1)-0.5 Bins.n_return(end)*2+2.5],'Xtick',[0:numel(Bins.n_return)*2],'Xticklabel',[ num2cell(Bins.n_return), {''} num2cell(Bins.n_return)]);
PT={'L','R'};
for k=1:2
    %subplot(n_rows,n_columns*2,k)
    current_bins=Bins.n_return + (k-1)*(Bins.n_return(end)+2);
    hold on
    B_H_L=bar(current_bins,+returns_to_revealed.hist.(['L' PT{k}]),1);
    set(B_H_L, 'FaceColor', [0 0 1])
    B_H_R=bar(current_bins,-returns_to_revealed.hist.(['R' PT{k}]),1);
    set(B_H_R, 'FaceColor', [0 1 0])
    text((current_bins(end)-current_bins(1))*0.25 + current_bins(1),bar_min*1.2, ['mean_R = ' num2str(round(returns_to_revealed.mean.(['R' PT{k}])*10)/10)]);
    text((current_bins(end)-current_bins(1))*0.25 + current_bins(1),bar_max*1.2, ['mean_L = ' num2str(round(returns_to_revealed.mean.(['L' PT{k}])*10)/10)]);
end


%% RT plots
subplot_assignments=[n_rows*2+2,n_rows*2+3];
ALL_RT_fieldnames={'LL','RL';'LR','RR'};
for k=1:2
    subplot(n_rows,n_columns,subplot_assignments(k))
    RT_fieldnames=ALL_RT_fieldnames(k,:);
    RT_lables={};
    hold on
    for m=1:numel(RT_fieldnames)
        current_structure=RTs.histograms.(RT_fieldnames{m});
        current_median=RTs.median.(RT_fieldnames{m});
        hRT(m)=plot(Bins.reaction_time,norm_hist(current_structure),'color',GLO.space_colors_RT(m,:),'linewidth',GLO.linewidth);
        line([current_median current_median],[0 1],'color',GLO.space_colors_RT(m,:),'linewidth',GLO.linewidth-1.5);
        RT_lables=[RT_lables; {[num2str(round(current_median)) ' ms']}];
    end
    legend(hRT,RT_lables,'Location','best');
    title('Reaction time','fontsize',GLO.title_font)
    ylim([0 1]);
    ylabel('Probability','fontsize',GLO.ylabel_font);
    xlabel('time [ms]','fontsize',GLO.xlabel_font);
end

%% First response plots
subplot(n_rows,n_columns,n_columns*2+1)
All_N_L=sum([revealed_tar_pos.first.LL,  revealed_tar_pos.first.LR]);
All_N_R=sum([revealed_tar_pos.first.RL, revealed_tar_pos.first.RR]);
hold on
bH(1)=barh(2,-(revealed_tar_pos.first.LL/All_N_L),'Facecolor',GLO.samplecolors(1,:));
bH(2)=barh(2, (revealed_tar_pos.first.LR/All_N_L),'Facecolor',GLO.samplecolors(1,:));
bH(3)=barh(1,-(revealed_tar_pos.first.RL/All_N_R),'Facecolor',GLO.samplecolors(2,:));
bH(4)=barh(1, (revealed_tar_pos.first.RR/All_N_R),'Facecolor',GLO.samplecolors(2,:));

text(-0.8,2,[num2str(revealed_tar_pos.first.LL)]);
text( 0.8,2,[num2str(revealed_tar_pos.first.LR)]);
text(-0.8,1,[num2str(revealed_tar_pos.first.RL)]);
text( 0.8,1,[num2str(revealed_tar_pos.first.RR)]);

title('1st response L   1st response R','fontsize',GLO.title_font)
ylabel('Sample','fontsize',GLO.ylabel_font)
Ylabels= {'R','L'};
set(gca,'ytick',[1:2],'ylim',[0 3],'xlim',[-1 1],'yticklabel',[Ylabels],'xtick',[-1:0.5:1],'xticklabel',[1 0.5 0 0.5 1]);


%% Distribution of search intervals
subplot(n_rows,n_columns,n_columns*3+1)
%ALL_RT_fieldnames={'LL','RL';'LR','RR'}
hold on
sid_lables={};
for m=1:numel(inspection_time.mean_search_to_reveal)
    %hRT(m)=plot(Bins.reaction_time,norm_hist(current_structure),'color',GLO.space_colors_RT(m,:),'linewidth',GLO.linewidth);
    sid(m)=plot(Bins.search_time,norm_hist(inspection_time.hist_search_to_reveal{m}),'color',GLO.search_intervals(m,:),'linewidth',GLO.linewidth);
    line([inspection_time.mean_search_to_reveal(m) inspection_time.mean_search_to_reveal(m)],[0 1],'color',GLO.search_intervals(m,:),'linewidth',GLO.linewidth-1.5);
    sid_lables=[sid_lables; {[num2str(round(inspection_time.mean_search_to_reveal(m))) ' ms']}];
end
legend(sid,sid_lables,'Location','best');
title('Time to reveal positions','fontsize',GLO.title_font);
ylim([0 1]);
ylabel('Probability','fontsize',GLO.ylabel_font);
xlabel('time [ms]','fontsize',GLO.xlabel_font);

%% Dwell time per position plots

ALL_RT_fieldnames={'LL','RL';'LR','RR'};
for k=1:numel(it_subplot_in_row)
    DT_fieldnames=ALL_RT_fieldnames(k,:);
    start_idx=n_columns*(3);
    subplot(n_rows,n_columns,start_idx+it_subplot_in_row(k))
    current_structure=inspection_time.hist_per_position;
    current_median=inspection_time.median_max_per_position;
    lables={'left','right'};
    current_title=['Max dwell time ' lables{k} ' targets'];
    hold on
    for m=1:numel(DT_fieldnames)
        hr(m)=plot(Bins.inspection_mean,norm_hist(current_structure.(DT_fieldnames{m})),'color',GLO.samplecolors(m,:),'linewidth',GLO.linewidth);
        line([current_median.(DT_fieldnames{m}) current_median.(DT_fieldnames{m})],[0 1],'color',GLO.samplecolors(m,:),'linewidth',GLO.linewidth-1.5);
    end
    Inspection_lables={[num2str(round(current_median.(DT_fieldnames{1}))) ' ms'],[num2str(round(current_median.(DT_fieldnames{2}))) ' ms']};
    legend(hr,Inspection_lables,'Location','best');
    
    title(current_title,'fontsize',GLO.title_font)
    ylim([0 1]);
    ylabel('Probability','fontsize',GLO.ylabel_font);
    xlabel('time [ms]','fontsize',GLO.xlabel_font);
end



%% Histogram success over time
subplot(n_rows,n_columns,n_columns*4+1)
plot(Bins.exploration_time,Summary.successrate_over_time,'color',[0 1 0],'linewidth',GLO.linewidth)
title('Successrate across time','fontsize',GLO.title_font);
ylim([0 1]);
ylabel('Hitrate','fontsize',GLO.ylabel_font);
xlabel('Exploration time [ms]','fontsize',GLO.xlabel_font);


%% Exploration time
if isfield(inspection_time,'raw')
    FN_exp={'explore','ITIexpl'};
    Exp_titles={'during trial','during ITI'};
    for k=1:2
        fn=FN_exp{k};
    %% total exploration
    subplot(n_rows,n_columns,n_columns*4+1+k)
    hold on
    bH(1)=barh(2,-(nanmean([Summary.(fn).LL])/nanmean([Summary.(fn).total])),'Facecolor',GLO.samplecolors(1,:));
    bH(2)=barh(2, (nanmean([Summary.(fn).LR])/nanmean([Summary.(fn).total])),'Facecolor',GLO.samplecolors(1,:));
    bH(3)=barh(1,-(nanmean([Summary.(fn).RL])/nanmean([Summary.(fn).total])),'Facecolor',GLO.samplecolors(2,:));
    bH(4)=barh(1, (nanmean([Summary.(fn).RR])/nanmean([Summary.(fn).total])),'Facecolor',GLO.samplecolors(2,:));
    text(-0.8,2,[num2str(round(nanmean([Summary.(fn).LL])))]);
    text( 0.8,2,[num2str(round(nanmean([Summary.(fn).LR])))]);
    text(-0.8,1,[num2str(round(nanmean([Summary.(fn).RL])))]);
    text( 0.8,1,[num2str(round(nanmean([Summary.(fn).RR])))]);
    
    title(['Exploration ' Exp_titles{k} ', total = ' num2str(round(nanmean([Summary.(fn).total]))) ' ms'],'fontsize',GLO.title_font)
    ylabel('Sample','fontsize',GLO.ylabel_font)
    Ylabels= {'R','L'};
    set(gca,'ytick',[1:2],'ylim',[0 3],'xlim',[-1 1],'yticklabel',[Ylabels],'xtick',[-1:0.5:1],'xticklabel',[1 0.5 0 0.5 1]);
    end
end

% mtit(summary_figure2,print_out, 'fontsize', GLO.main_title_font, 'xoff', -0.05, 'yoff', 0.03, 'color',[0 0 0]);
% pdf_export(Disp,'positions')

title_and_save(summary_figure2,print_out,'positions',Disp);

%% Figure 3: time weighted raw eye positions (Exploration pattern)
if isfield(inspection_time,'raw')
    summary_figure3=figure('units','normalized','outerposition',[0 0 1 1]);
    colormap(gca,flipud(jet));
    for k=1:8
        subplot(4,2, k)
        image(Pattern{k});
        steps_circle=0:0.01:2*pi;
        circle_x=Disp.tar_radius/Bins.stepsize*cos(steps_circle);
        circle_y=Disp.tar_radius/Bins.stepsize*sin(steps_circle);
        box off
        set(gca,...
            'xtick', ([Bins.raw_2d{1}(1):10:Bins.raw_2d{1}(end)]'-Bins.raw_2d{1}(1))/Bins.stepsize+1,...
            'ytick', ([Bins.raw_2d{2}(1):10:Bins.raw_2d{2}(end)]'-Bins.raw_2d{2}(1))/Bins.stepsize+1,...
            'xticklabel',[Bins.raw_2d{1}(1):10:Bins.raw_2d{1}(end) ],...
            'yticklabel',fliplr([Bins.raw_2d{2}(1):10:Bins.raw_2d{2}(end) ]))%, 'YAxisLocation', 'right', 'XAxisLocation', 'top')
        %axis equal
        for pos=1:size(possible_positions,1)
            if ismember(possible_positions(pos),Disp.pattern_target_positions{k})
                line(circle_x+(real(possible_positions(pos))-Bins.raw_2d{1}(1))/Bins.stepsize,circle_y+(-imag(possible_positions(pos))-Bins.raw_2d{2}(1))/Bins.stepsize,'Color','w','linewidth',1);
                %             else
                %             line(circle_x+(real(possible_positions(pos))-Bins.raw_2d{1}(1))/Bins.stepsize,circle_y+(imag(possible_positions(pos))-Bins.raw_2d{2}(1))/Bins.stepsize,'Color','w','linewidth',1);
            end
        end
        add_normalized_colorbar(Summary.max_pattern_samples,'ms');
        title(Disp.pattern_labels{k,2},'fontsize',GLO.title_font, 'interpreter', 'none');
        ylabel(Disp.pattern_labels{k,1},'fontsize',GLO.title_font, 'interpreter', 'none');
    end
    
    title_and_save(summary_figure3,print_out,'Pollock',Disp);
    %     mtit(summary_figure3,print_out, 'fontsize', GLO.main_title_font, 'xoff', -0.05, 'yoff', 0.03, 'color',[0 0 0]);
    %     pdf_export(Disp,'Pollock')
end

%% Figure 4: Raw heatmap
if isfield(inspection_time,'raw')
    summary_figure4=figure('units','normalized','outerposition',[0 0 1 1]);
    %image(rot90(inspection_time.raw)/max(max(inspection_time.raw))*64.5)
    steps_circle=0:0.01:2*pi;
    circle_x=Disp.tar_radius/Bins.stepsize*cos(steps_circle);
    circle_y=Disp.tar_radius/Bins.stepsize*sin(steps_circle);
    FN_EX={'heat_success','heat_error','ITI'};
    raw_titles={'Sucess','Errors','ITI'};
    for f=1:numel(FN_EX)
        FN=FN_EX{f};
        subplot(3,2,2*f-1);
        imagesc(single(rot90(inspection_time.(FN))/max(max(inspection_time.(FN)))),[0 1]);
        box off
        set(gca,...
            'xtick', ([Bins.raw_2d{1}(1):10:Bins.raw_2d{1}(end)]'-Bins.raw_2d{1}(1))/Bins.stepsize+1,...
            'ytick', ([Bins.raw_2d{2}(1):10:Bins.raw_2d{2}(end)]'-Bins.raw_2d{2}(1))/Bins.stepsize+1,...
            'xticklabel',[Bins.raw_2d{1}(1):10:Bins.raw_2d{1}(end) ],...
            'yticklabel',fliplr([Bins.raw_2d{2}(1):10:Bins.raw_2d{2}(end) ]))%, 'YAxisLocation', 'right', 'XAxisLocation', 'top')
        
        for pos=1:size(possible_positions,1)
            line(circle_x+(real(possible_positions(pos))-Bins.raw_2d{1}(1))/Bins.stepsize,circle_y+(-imag(possible_positions(pos))-Bins.raw_2d{2}(1))/Bins.stepsize,'Color','g','linewidth',1);
        end
        colormap(gca,'hot');
        
        max_IT=round(max(max(inspection_time.raw)));
        add_normalized_colorbar(max_IT,'ms');
        title(raw_titles{f})
        
        subplot(3,2,2*f);
        histogram_in_x=sum(rot90(inspection_time.(FN)),1)/sum(sum(inspection_time.(FN)));
        plot([Bins.raw_2d{1}],histogram_in_x,'k');
        ylabel('Fraction of total time');
        xlabel('Horizontal position [°]');
        
    end
    
    title_and_save(summary_figure4,print_out,'heatmap',Disp);
    
end
end

function normalized_hist=norm_hist(histogram_to_normalize)
normalized_hist=histogram_to_normalize/sum(histogram_to_normalize);
end

function pointList=CalculateConvexPointList(center,a_ellipse,convexity,convex_sides)
convexity_sign=sign(convexity);
b_ellipse=abs(convexity)*a_ellipse;
b_ellipse_reference=0.3*a_ellipse;

steps_for_interpolation=100;
euclidian_distance=(-a_ellipse):2*a_ellipse/steps_for_interpolation:(a_ellipse);
%corner_positions=[-1,-1;1,-1;1,1;-1,1].*a_ellipse;

Half_circle_Area=a_ellipse^2*pi/2;
rect_b=(Half_circle_Area-a_ellipse*b_ellipse*pi*convexity_sign/2)/(2*a_ellipse);
rect_b_reference=(Half_circle_Area-a_ellipse*b_ellipse_reference*pi*convexity_sign/2)/(2*a_ellipse);
%rect_ratio=Half_circle_Area/(2*a_ellipse^2);

tmp_bow_parameter=sqrt((b_ellipse)^2.*(1-euclidian_distance'.^2/a_ellipse^2));
tmp_bow_reference=sqrt((b_ellipse_reference)^2.*(1-euclidian_distance'.^2/a_ellipse^2));

bow_vector=[euclidian_distance',(tmp_bow_parameter.*convexity_sign.*-1-rect_b)];
bow_vector_reference=[euclidian_distance',(tmp_bow_reference.*convexity_sign.*-1-rect_b_reference)];

if strcmp(convex_sides,'LR')|| strcmp(convex_sides,'R') || strcmp(convex_sides,'L')
    bow_vector=[bow_vector(:,2)*-1,bow_vector(:,1)];
    bow_vector_reference=[bow_vector_reference(:,2)*-1,bow_vector_reference(:,1)];
end

switch convex_sides
    case 'T'
        pointList=[bow_vector;bow_vector_reference*-1];
    case 'B'
        pointList=[bow_vector_reference;bow_vector*-1];
    case 'TB'
        pointList=[bow_vector;bow_vector*-1];
    case 'R'
        pointList=[bow_vector;bow_vector_reference*-1];
    case 'L'
        pointList=[bow_vector_reference;bow_vector*-1];
    case 'LR'
        pointList=[bow_vector;bow_vector*-1];
end
pointList=pointList+repmat(center,size(pointList,1),1);
end

function pie_chart_plot(plot_only_LR,pie_chart_data,pie_chart_lables)
global GLO
hold on
r=1;
r_dec=1/numel(pie_chart_data);
p_vec=[numel(pie_chart_data):-1:1];
for circles=numel(pie_chart_data):-1:1
    if circles>1
        x=permute(pie_chart_data{circles},p_vec);
    else
        x=pie_chart_data{circles};
    end
    x=x(:);
    if plot_only_LR
        current_colors=repmat(GLO.space_colors_LRN,numel(x)/3,1);
    else
        current_colors=repmat(GLO.space_colors_allpositions,numel(x)/(size(GLO.space_colors_allpositions,1)),1);
    end
    piechart(circles).handle=DAG_pie_chart(x,r,current_colors,[0,0],(1-r_dec/2));
    r=r-r_dec;
    p_vec(1)=[];
end
n_pie_chart_positions=numel(pie_chart_data{1})-1;
leg_patterns=legend(piechart(1).handle(1:2:2*n_pie_chart_positions+1),pie_chart_lables([1:n_pie_chart_positions,end]),'Location','NorthEast');
end

function heat_map_plot(first_target_data,first_target_nonnan,possible_positions,varargin)
global GLO
N_steps=255;

colormap(gca,cool);
colors_heat_map=cool(N_steps);

plot(0,GLO.fix_y,'+','MarkerSize',20,'Color','k');
hold on
for pos_idx=1:numel(possible_positions)
    selected_pos_rate=round(first_target_data{1}(pos_idx)/first_target_nonnan{1}(pos_idx)*(N_steps-1))+1;
    if ~isnan(selected_pos_rate)
        plot(real(possible_positions(pos_idx)),imag(possible_positions(pos_idx)),'o','MarkerSize',16,'MarkerFaceColor',colors_heat_map(selected_pos_rate,:),'MarkerEdgeColor',colors_heat_map(selected_pos_rate,:))
        
        
        % plot(real(possible_positions(pos_idx)),imag(possible_positions(pos_idx)),'o','MarkerSize',20,'MarkerFaceColor',colors_used(pos_idx,:),'MarkerEdgeColor',colors_used(pos_idx,:))
        if numel(varargin)>0
            text(real(possible_positions(pos_idx))-max(real(possible_positions))/9,imag(possible_positions(pos_idx)),num2str(varargin{1}(pos_idx)),'fontsize',20)
        end
        
    end
end
%%,'ytick',[min(imag(possible_positions)):5:max(imag(possible_positions))],'xtick',[min(real(possible_positions)):5:max(real(possible_positions))]);%'ytick',[1:6],'ylim',[0 7]
%xlabel('X [deg]')
%ylabel('Y [deg]')


%automatic_ylim=get(gca,'ylim');
axis equal
y_range_min=20;
y_range=max(max(imag(possible_positions))- min(imag(possible_positions)),y_range_min);
set(gca,'ylim',[min(imag(possible_positions))-y_range/5 max(imag(possible_positions))+y_range/5]);

x_range=max(real(possible_positions))- min(real(possible_positions));
set(gca,'xlim',[min(real(possible_positions))-x_range/5 max(real(possible_positions))+x_range/5]);

%ylabels=0:round(max(first_target_nonnan{1})/5):round(max(first_target_nonnan{1}));

%colormap(gca,flipud(jet));
add_normalized_colorbar(100,'');
%
% ylabels=0:20:100;
% colormap('autumn');
% a=colorbar('eastoutside');
% ylims=get(a,'Ylim');
% %ylabelticks=ylabels/max(first_target_nonnan{1})*max(ylims);
% ylabelticks=ylabels/100*max(ylims);
%
% set(a,'Ytick',ylabelticks,'Yticklabel',ylabels);
% %set(get(a,'Children'),'Ydata',[0 100])

end

function add_normalized_colorbar(maximum,units)
ylabels=round([maximum/10:maximum/10:maximum])';
cbh=colorbar('eastoutside');
ylims=get(cbh,'Ylim');
ylabelticks=ylabels/maximum*max(ylims);
ylabels=cellstr([num2str(ylabels) repmat([' ' units],size(ylabels,1),1)]);

ylabels{end}=['>' ylabels{end}];
set(cbh,'Ytick',ylabelticks,'Yticklabel',ylabels)
end

function title_and_save(figure_handle,plot_title,title,Disp)
global GLO
mtit(figure_handle,plot_title, 'fontsize', GLO.main_title_font, 'xoff', -0.05, 'yoff', 0.03, 'color',[0 0 0]);
stampit;
if Disp.settings.export_pdf
    wanted_size=[50 30];
    set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
    export_fig([Disp.settings.export_pdf_filename '_' title], '-pdf','-transparent') % append to existing pdf
    close gcf
end
end

% function pdf_export(Disp,title)
% if Disp.settings.export_pdf && Disp.settings.append_pdf
%     export_fig([Disp.settings.export_pdf_filename '_' title], '-pdf','-transparent','-append') % append to existing pdf
%     close gcf
% elseif Disp.settings.export_pdf
%     export_fig([Disp.settings.export_pdf_filename '_' title], '-pdf','-transparent') % append to existing pdf
%     close gcf
% end
% end

function par=mtit(varargin)
%MTIT		creates a major title in a figure with many axes
%
%		MTIT
%		- creates a major title above all
%		  axes in a figure
%		- preserves the stack order of
%		  the axis handles
%
%SYNTAX
%-------------------------------------------------------------------------------
%		P = MTIT(TXT,[OPT1,...,OPTn])
%		P = MTIT(FH,TXT,[OPT1,...,OPTn])
%
%INPUT
%-------------------------------------------------------------------------------
%    FH	:	a valid figure handle		[def: gcf]
%   TXT	:	title string
%
% OPT	:	argument
% -------------------------------------------
%  xoff	:	+/- displacement along X axis
%  yoff	:	+/- displacement along Y axis
%  zoff	:	+/- displacement along Z axis
%
%		title modifier pair(s)
% -------------------------------------------
%   TPx	:	TVx
%		see: get(text) for possible
%		     parameters/values
%
%OUTPUT
%-------------------------------------------------------------------------------
% par	:	parameter structure
%  .pos :	position of surrounding axis
%   .oh	:	handle of last used axis
%   .ah :	handle of invisible surrounding axis
%   .th :	handle of main title
%
%EXAMPLE
%-------------------------------------------------------------------------------
%	subplot(2,3,[1 3]);		title('PLOT 1');
%	subplot(2,3,4); 		title('PLOT 2');
%	subplot(2,3,5); 		title('PLOT 3');
%	axes('units','inches',...
%	     'color',[0 1 .5],...
%	     'position',[.5 .5 2 2]);	title('PLOT 41');
%	axes('units','inches',...
%	     'color',[0 .5 1],...
%	     'position',[3.5 .5 2 2]);	title('PLOT 42');
%	shg;
%	p=mtit('the BIG title',...
%	     'fontsize',14,'color',[1 0 0],...
%	     'xoff',-.1,'yoff',.025);
% % refine title using its handle <p.th>
%	set(p.th,'edgecolor',.5*[1 1 1]);

% created:
%	us	24-Feb-2003		/ R13
% modified:
%	us	24-Feb-2003		/ CSSM
%	us	06-Apr-2003		/ TMW
%	us	13-Nov-2009 17:38:17

defunit='normalized';
if	nargout
    par=[];
end

% check input
if	nargin < 1
    help(mfilename);
    return;
end
if	isempty(get(0,'currentfigure'))
    disp('MTIT> no figure');
    return;
end

vl=true(size(varargin));
if	ischar(varargin{1})
    vl(1)=false;
    figh=gcf;
    txt=varargin{1};
elseif	any(ishandle(varargin{1}(:)))		&&...
        ischar(varargin{2})
    vl(1:2)=false;
    figh=varargin{1};
    txt=varargin{2};
else
    error('MTIT> invalid input');
end
vin=varargin(vl);
[off,vout]=get_off(vin{:});

% find surrounding box
ah=findall(figh,'type','axes');
if	isempty(ah)
    disp('MTIT> no axis');
    return;
end
oah=ah(1);

ou=get(ah,'units');
set(ah,'units',defunit);
ap=get(ah,'position');
if	iscell(ap)
    ap=cell2mat(get(ah,'position'));
end
ap=[	min(ap(:,1)),max(ap(:,1)+ap(:,3)),...
    min(ap(:,2)),max(ap(:,2)+ap(:,4))];
ap=[	ap(1),ap(3),...
    ap(2)-ap(1),ap(4)-ap(3)];

% create axis...
xh=axes('position',ap);
% ...and title
th=title(txt,'interpreter', 'none',vout{:});
tp=get(th,'position');
set(th,'position',tp+off);
set(xh,'visible','off','hittest','on');
set(th,'visible','on');

% reset original units
ix=find(~strcmpi(ou,defunit));
if	~isempty(ix)
    for	i=ix(:).'
        set(ah(i),'units',ou{i});
    end
end

% ...and axis' order
uistack(xh,'bottom');
axes(oah);				%#ok

if	nargout
    par.pos=ap;
    par.oh=oah;
    par.ah=xh;
    par.th=th;
end

    function	[off,vout]=get_off(varargin)
        
        % search for pairs <.off>/<value>
        
        off=zeros(1,3);
        io=0;
        for	mode={'xoff','yoff','zoff'};
            ix=strcmpi(varargin,mode);
            if	any(ix)
                io=io+1;
                yx=find(ix);
                ix(yx+1)=1;
                off(1,io)=varargin{yx(end)+1};
                varargin=varargin(xor(ix,1));
            end
        end
        vout=varargin;
    end
end




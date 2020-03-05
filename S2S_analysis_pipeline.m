clear
main_pdf_folder='Y:\Projects\STS_PUL_MEM_S2S\behavior\Human_Alex_2018\';

                wanted_size=[50 30];
% all_comparisons={'SHAM vs TMS young','preSHAM vs preTMS young','postSHAM vs postTMS young','preSHAM vs postSHAM young','preTMS vs postTMS young','7 targets','3 targets 1 cue Y vs Cornelius','3 targets 1 cue Y vs O',...
%     '3 targets 2 cues Y vs O','3 targets naive Y vs O','3 targets Y 1 vs 2 cues','3 targets O 1 vs 2 cues','GE pre vs post mem_per','Patient'};
%all_comparisons={'3 targets 1 cue Y vs O','3 targets Y 1 vs 2 cues','3 targets O 1 vs 2 cues','3 targets 2 cues Y vs O','3 targets naive Y vs O','3 targets 1 cue Y vs Cornelius'};

%all_comparisons={'3 targets 1 cue Y vs O'};
%all_comparisons={'3 targets 1 cue Y vs O','3 targets Y 1 vs 2 cues','3 targets O 1 vs 2 cues','3 targets 2 cues Y vs O','3 targets naive Y vs O'};
%all_comparisons={'3 targets 1 cue ideal vs perceptual neglect', '3 targets 1 cue ideal vs intentional neglect','3 targets 1 cue ideal vs object related neglect'};
all_comparisons={'Patient'};
%all_comparisons={'3 targets 1 cue Y vs O','3 targets Y 1 vs 2 cues','3 targets O 1 vs 2 cues','3 targets 2 cues Y vs O'};
for comparison_cell=all_comparisons
    
    clear Aiminput main_data_path Sum_ses Sum_norm Sum_tot subjects M2S_Analyze_settings_per_group Selection_cell; 
    M2S_Analyze_settings_per_group=struct;
comparison=comparison_cell{:}; 
 %   comparison='preTMS vs postTMS young';  % <<- CHANGE HERE WHAT YOU WANT TO ANALYZE

M2S_Analyze_settings.placeholders_visible=0;
M2S_Analyze_settings.export_pdf=1;
M2S_Analyze_settings.append_pdf=0;
M2S_Analyze_settings.effector=0; % always 0 (so far)
M2S_Analyze_settings.type=6; % always 6
M2S_Analyze_settings.plot_only_LR=0; % for 7 targets task
M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
%[M2S_Analyze_settings_per_group(1:2).limit_convexity]=deal(0.3);
Selection_cell{1}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',40,'sac_end_t',25,'downsampling',1,'eyetracker_sample_rate',60,'smoothing_samples',30,'i_sample_rate',1000,'correct_offset',1};
Selection_cell{2}=Selection_cell{1};
isrealdata={1,1};
[M2S_Analyze_settings_per_group(1:2).choice]=deal(0);

%Selection_cell={'display',0,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',40,'sac_end_t',25,'downsampling',1,'eyetracker_sample_rate',60,'smoothing_samples',6,'i_sample_rate',200,'correct_offset',0};

Inputsequal={};
Inputsrange={};
% Inputsequal={'Hand_dominance','r'};
% Inputsrange={'Varied_curvature_concave',-0.5,-0.4;'Session',20150000,20150831};
Inputslist={};



Subdirectories_to_evaluate = '1:numel(subdir)';
groups=1:2;
n_subjects_simulated=1;
switch comparison
    
    
    
    case 'SHAM vs TMS young'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1,1}=['preSHAM'];
        Aiminput{1,2}=['postSHAM'];
        Aiminput{2,1}=['preTMS'];
        Aiminput{2,2}=['postTMS'];
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings.placeholders_visible=1;
        
    case 'preSHAM vs preTMS young'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1}=['preSHAM'];
        Aiminput{2}=['preTMS'];
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings.placeholders_visible=1;
        
    case 'preTMS vs postTMS young'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1}=['preTMS'];
        Aiminput{2}=['postTMS'];
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings.placeholders_visible=1;
        
    case 'preSHAM vs postSHAM young'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1}=['preSHAM'];
        Aiminput{2}=['postSHAM'];
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings.placeholders_visible=1;
        
    case 'postSHAM vs postTMS young'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1}=['postSHAM'];
        Aiminput{2}=['postTMS'];
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings.placeholders_visible=1;
        
    case '7 targets'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\7 targets old subjects'];
        %main_data_path{2}=['Y:\Data\Human\Search2Sample\4 targets young TMS'];
        Aiminput{1}=['main task 1 cue'];
        %Aiminput{2}=['postTMS'];
        %stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        % ADDITIONAL CHANGES:
        % M2S_Analyze_settings.plot_only_LR=1;
        % g = 1
        % comment 'S2S_group_comparison' and the following 4 lines
        groups=1;
        M2S_Analyze_settings.plot_only_LR=1;
        M2S_Analyze_settings.placeholders_visible=1;
        
    case '3 targets 1 cue Y vs Cornelius'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 1 cue'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings_per_group(1).limit_convexity=0.3;
        M2S_Analyze_settings_per_group(2).limit_convexity=0.4;
        Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};
 
    case '3 targets Cornelius 1 cue vs 2 cues'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 2 cues'];
        M2S_Analyze_settings_per_group(1).choice=0;
        M2S_Analyze_settings_per_group(2).choice=1;
        stats_test='paired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings_per_group(1).limit_convexity=0.4;
        M2S_Analyze_settings_per_group(2).limit_convexity=0.4;
        Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};
        Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};

     case '3 targets 1 cue ideal vs perceptual neglect'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 1 cue'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings_per_group(1).limit_convexity=0.3;
        M2S_Analyze_settings_per_group(2).limit_convexity=0.3;
        %Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};
        impairment{1}='ideal';
        impairment{2}='sensory_neglect_L';
        isrealdata{1}=0;       
        isrealdata{2}=0;         
        
    case '3 targets 1 cue ideal vs intentional neglect'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 1 cue'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings_per_group(1).limit_convexity=0.3;
        M2S_Analyze_settings_per_group(2).limit_convexity=0.3;
        %Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};
        impairment{1}='ideal';
        impairment{2}='motivational_neglect_L';
        isrealdata{1}=0;       
        isrealdata{2}=0;       
        
     case '3 targets 1 cue ideal vs object related neglect'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets Cornelius'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 1 cue'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        M2S_Analyze_settings_per_group(1).limit_convexity=0.3;
        M2S_Analyze_settings_per_group(2).limit_convexity=0.3;
        %Selection_cell{2}={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',300,'sac_end_t',150,'downsampling',1,'eyetracker_sample_rate',220,'smoothing_samples',10,'i_sample_rate',1000,'correct_offset',1};
        impairment{1}='ideal';
        impairment{2}='object_related_neglect_L';
        isrealdata{1}=0;       
        isrealdata{2}=0;           
        
    case '3 targets 1 cue Y vs O'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets old subjects'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 1 cue'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        
    case '3 targets 2 cues Y vs O'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets old subjects'];
        Aiminput{1}=['main task 2 cues'];
        Aiminput{2}=['main task 2 cues'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        subjects_with_2cues{1} = [1,2,3,4,5,6,7,8];
        subjects_with_2cues{2} = [1,3,4,6];
        M2S_Analyze_settings_per_group(1).choice=1;
        M2S_Analyze_settings_per_group(2).choice=1;
        % ADDITIONAL CHANGES:
        % d = subjects_with_2cues{g}
        Subdirectories_to_evaluate = 'subjects_with_2cues{g}';    

    case '3 targets naive Y vs O'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets old subjects'];
        Aiminput{1}=['naïve run'];
        Aiminput{2}=['naïve run'];
        stats_test='unpaired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        subjects_with_naiverun{1} = [1,2,3,5,6,7,8];
        subjects_with_naiverun{2} = [1,2,3,4,6];
        % ADDITIONAL CHANGES:
        % d = subjects_with_naiverun{g}
        Subdirectories_to_evaluate = 'subjects_with_naiverun{g}';
        
    case '3 targets Y 1 vs 2 cues'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets young subjects'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 2 cues'];
        M2S_Analyze_settings_per_group(1).choice=0;
        M2S_Analyze_settings_per_group(2).choice=1;
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for same subjects (TMS, same group) / unpaired for young vs old
        
    case '3 targets O 1 vs 2 cues'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets old subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\3 black targets old subjects'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task 2 cues'];
        M2S_Analyze_settings_per_group(1).choice=0;
        M2S_Analyze_settings_per_group(2).choice=1;
        stats_test='paired_ttest'; % 'unpaired_ttest'   % paired for same subjects (TMS, same group) / unpaired for young vs old
        old_with_2cues = [1,3,4,6];
        % ADDITIONAL CHANGES:
        % d = old_with_2cues
        Subdirectories_to_evaluate = 'old_with_2cues';
        
    case 'GE pre vs post mem_per'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\7 targets old subjects'];
        main_data_path{2}=['Y:\Data\Human\Search2Sample\7 targets old subjects'];
        Aiminput{1}=['main task 1 cue'];
        Aiminput{2}=['main task mem 1 cue'];
        stats_test='paired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        folders_GE{1} = [1];
        folders_GE{2} = [1];
        % ADDITIONAL CHANGES:
        % d = folders_GE{g}
        Subdirectories_to_evaluate = 'folders_GE{g}';
        M2S_Analyze_settings.plot_only_LR=1;
        M2S_Analyze_settings.placeholders_visible=1;
        
        
    case 'Patient'
        main_data_path{1}=['Y:\Data\Human\Search2Sample\3 black targets patients'];
        %main_data_path{2}=['C:\Users\akratzenberg\Desktop\NSall\GE'];
        Aiminput{1}=['main task 1 cue'];
        %Aiminput{2}=['main task 1 cue'];
        %stats_test='paired_ttest'; % 'paired_ttest'   % paired for TMS (= same subjects) / unpaired for young vs old
        % ADDITIONAL CHANGES:
        % g = 1
        groups=1;
        
        
end

main_pdf_path{1}=[main_pdf_folder comparison];
main_pdf_path{2}=[main_pdf_folder comparison];

if ~exist(main_pdf_path{1},'file')
    mkdir(main_pdf_folder, comparison)
end
for p=1:size(Aiminput,1)
    for g=groups
        if isfield(M2S_Analyze_settings_per_group, 'limit_convexity')
        M2S_Analyze_settings.limit_convexity=M2S_Analyze_settings_per_group(g).limit_convexity;
        end
        
        M2S_Analyze_settings.choice=M2S_Analyze_settings_per_group(g).choice;
        clear M2S_output;
        if isrealdata{g}
            subdirectories=dir(main_data_path{g});
            subdirectories(1:2)=[];
            subdir={subdirectories([subdirectories.isdir]).name};
            subdir_indexes=eval(Subdirectories_to_evaluate);
            subjects{p}=subdir(subdir_indexes);
        else
            subdir_indexes=1:n_subjects_simulated;
            subdir=repmat({['simulated_' impairment{g} '_']},1,n_subjects_simulated)';            
            [subjects{p}]=deal(cellstr([cell2mat(subdir),num2str([1:numel(subdir)]')])');
        end
        for d=subdir_indexes %1:numel(subdir)  %subjects_with_naiverun{g} subjects_with_2cues{g}  ;insert here which subject(s) you want to analyze
           
            if isrealdata{g}
            datapath=[main_data_path{g} filesep subdir{d} filesep];
            xlspath=[main_data_path{g} filesep subdir{d} filesep subdir{d} '.xls'];
            %DAG_create_mastertable_sheet(xlspath)
            filelist_formatted = DAG_get_filelist_from_xls(Aiminput(p,g),Inputsequal,Inputsrange,Inputslist, datapath, xlspath);
                [out, ~, ~]= monkeypsych_analyze_working(filelist_formatted,Selection_cell{g});
            else
                out=S2S_simulate_deficits(impairment{g},1000);
            end
            M2S_Analyze_settings.export_pdf_filename=[main_pdf_path{g} filesep subdir{d} ' ' Aiminput{p,g}];
            M2S_output(d)=S2S_session_analysis(out{1},M2S_Analyze_settings);
            S2S_session_plot(M2S_output(d))
        end
        % M2S_output(4)=[];
        % M2S_output(1)=[];
        
        [Sum_ses(p,g),Sum_norm(p,g),Sum_tot(p,g)]=S2S_group_summary(M2S_output,subjects{p});
        if M2S_Analyze_settings.export_pdf
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([main_pdf_path{g} filesep comparison ' group ' num2str(g) ' ' Aiminput{p,g} ' ' 'GroupSummary2'], '-pdf', '-transparent') % pdf by run
            close(gcf);
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([main_pdf_path{g} filesep comparison ' group ' num2str(g) ' ' Aiminput{p,g} ' ' 'GroupSummary1'], '-pdf', '-transparent') % pdf by run
            close(gcf);
        end
    end
end


if numel(groups)==2
    S2S_group_comparison(Sum_ses,Sum_norm,Sum_tot,stats_test,comparison,subjects{p})
    if M2S_Analyze_settings.export_pdf && numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        export_fig([main_pdf_path{g} filesep comparison ' ' 'comparison3'], '-pdf', '-transparent') % pdf by run
        close(gcf);
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        export_fig([main_pdf_path{g} filesep comparison ' ' 'comparison2'], '-pdf', '-transparent') % pdf by run
        close(gcf);
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        export_fig([main_pdf_path{g} filesep comparison ' ' 'comparison1'], '-pdf', '-transparent') % pdf by run
        close(gcf);
    elseif M2S_Analyze_settings.export_pdf
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        export_fig([main_pdf_path{g} filesep comparison ' ' 'comparison2'], '-pdf', '-transparent') % pdf by run
        close(gcf);
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        export_fig([main_pdf_path{g} filesep comparison ' ' 'comparison1'], '-pdf', '-transparent') % pdf by run
        close(gcf);
    end
end

%
%'motivational_neglect_L','object_related_neglect_L','sensory_neglect_L'
% all_impairments={'sensory_neglect_L'};
% for IM=1:numel(all_impairments)
%     impairment=all_impairments{IM};
%     generated_out=S2S_simulate_deficits(impairment);
%     M2S_output_generated=S2S_session_analysis(generated_out{1},M2S_Analyze_settings);
%     M2S_output(d+IM)=M2S_output_generated;
% S2S_session_plot(M2S_output(d+IM))
% if M2S_Analyze_settings.export_pdf
%
%         export_fig([Drive '\Projects\STS\Human\fakedata' filesep impairment ' ' 'summary2'], '-pdf', '-transparent') % pdf by run
%         close(gcf);
%         export_fig([Drive '\Projects\STS\Human\fakedata' filesep impairment ' ' 'summary1'], '-pdf', '-transparent') % pdf by run
%         close(gcf);
% end
% end

end
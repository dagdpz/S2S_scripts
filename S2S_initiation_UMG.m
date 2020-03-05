clear

user=getUserName;
main_path='Y:\Data\HumanS2S';
%xlsfilename='Y:\Personal\Hannah\Probandenliste3.xlsx';
xlsfilename='Y:\Projects\STS_PUL_MEM_S2S\behavior\Hannah\Probandenliste4.xlsx';
pdfpath=['Y:\Projects\STS_PUL_MEM_S2S\behavior\Hannah\plots_2018\'];

%types={'unilateral_curvature','bilateral_curvature'};
types={'bilateral_curvature'};
groups={'Y','E'};
% types={'bilateral_curvature'};
% groups={'L','S'};


%% SETTINGS
M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.export_pdf=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.append_pdf=0; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.limit_convexity=0.1; % only relevant if concatinate_option==1
M2S_Analyze_settings.type=6; %6
M2S_Analyze_settings.effector=0; %
M2S_Analyze_settings.placeholders_visible=0; %
M2S_Analyze_settings.choice=0; %

plot_specs.main_title_font=12;
plot_specs.export_pdf=1;
plot_specs.path=pdfpath;

% MPA keys
% Keys={'display',0,'summary',0,'runs_as_batches',0,'sac_ini_t',80,'sac_end_t',30,'eyetracker_sample_rate',60,'inferential_on',0,...
%       'nsacc_max',3,'keep_raw_data',1,'saccade_definition',4,'downsampling',1};
Keys={'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'sac_ini_t',40,'sac_end_t',25,'downsampling',1,'eyetracker_sample_rate',60,'smoothing_samples',30,'i_sample_rate',1000,'correct_offset',1};

% group comparison settings
comparison='3 targets 1 cue Y vs O'; % this just defines labels etc
stats_test='unpaired_ttest';         % which test to use for group comparison

%% reading in table and relevant table indexes
[~,~,complete_table]=xlsread(xlsfilename,'mastertable');
idx.Subject=DAG_find_column_index(complete_table,'Subject');
idx.Group=DAG_find_column_index(complete_table,'Group');
idx.Date=DAG_find_column_index(complete_table,'Date');
idx.Run=DAG_find_column_index(complete_table,'Run');
idx.Type=DAG_find_column_index(complete_table,'Type');
idx.use=DAG_find_column_index(complete_table,'use');

for g=1:numel(groups)
    group=groups{g};
    %% reduce table rows --> contains only current group, only used runs (used==true), and used types
    xltable=complete_table([true;...
        ismember(vertcat(complete_table{2:end,idx.use}),1) &...
        ismember(complete_table(2:end,idx.Type),types) &...
        ismember(complete_table(2:end,idx.Group),group)],:);
    
    %% find unique subjects & runs per subject, separated already by type
    clear dates_per_subject runs_per_subject_per_date
    subjects=unique(xltable(2:end,idx.Subject));
    for s=1:numel(subjects)
        s_idx=ismember(xltable(:,idx.Subject),subjects(s));
        dates_per_subject{s}=unique([xltable{s_idx,idx.Date}]);
        for t=1:numel(types)
            t_idx=ismember(xltable(:,idx.Type),types(t));
            for d=1:numel(dates_per_subject{s})
                d_idx=[false; ismember(vertcat(xltable{2:end,idx.Date}),dates_per_subject{s}(d))];
                runs_per_subject_per_date{t}{s}{d}=strcat('_', xltable(t_idx & s_idx & d_idx,idx.Run), '.');
            end
        end
    end
    
    %monkeypsych_clean_data(main_path,[min(unique([dates_per_subject{:}])) max(unique([dates_per_subject{:}]))]);
    for t=1:numel(types)
        type=types{t};
        clear M2S_output
        for s=1:numel(subjects)
            
            %% create filelist for each subject (per type)
            subject=subjects{s};
            dates=dates_per_subject{s};
            allfiles={};
            for d=1:numel(dates)
                date=dates(d);
                runs=runs_per_subject_per_date{t}{s}{d};
                files_per_date=dir([main_path filesep num2str(date) filesep '*.mat']);
                files_per_date={files_per_date.name};
                files_per_date=files_per_date(cellfun(@(x) ~isempty(strfind(x,subject)),files_per_date));
                allfiles=[allfiles; files_per_date(cellfun(@(x) any(cellfun(@(y) ~isempty(strfind(x,y)),runs)),files_per_date))'];
            end
            if isempty(allfiles)
                continue;
            end
            %% run MPA, STS_session_analysis, and plot for each subject in the current group
            filelist={strcat([main_path filesep num2str(date) filesep], allfiles)};
            monkeypsych_analyze_input=[filelist; repmat({Keys},1,numel(filelist))];
            M2S_Analyze_settings.export_pdf_filename=[pdfpath group '_' subject '_' type ];
            [MPA_out, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
            M2S_output(s)=S2S_session_analysis(MPA_out{1},M2S_Analyze_settings);
            S2S_session_plot(M2S_output(s));
            
        end
        
        %% Group summary 
        plot_specs.figure_name=[group '_' type];
        plot_specs.individual_labels=subjects;
        [Sum_ses(t,g),Sum_norm(t,g),Sum_tot(t,g)]=S2S_group_summary_UMG(M2S_output,plot_specs);
    end
end


%% Group comparison
for t=1:numel(types)
    plot_specs.figure_name=types{t};
    S2S_group_comparison_UMG(Sum_ses(t,:),Sum_norm(t,:),Sum_tot(t,:),stats_test,comparison,plot_specs)
end

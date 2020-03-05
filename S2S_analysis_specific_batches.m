%% settings
for Ta={'dPul','STP'}
    for BL={'use_S2S_pre'} %,'use_S2S_bas'}
        clear M2S_output Ses_ses Ses_norm Ses_tot Sum_ses Sum_norm Sum_tot
        Target=Ta{:};
        Baseline=BL{:};
        % Target='STP';%STP%'dPul'
        % Baseline='use_S2S_pre';%'use_S2S_pre' %'use_S2S_bas'
        
        monkey='Cornelius';
        dag_drive_IP=DAG_get_server_IP;
        current_folder=[dag_drive_IP 'Projects\STS_PUL_MEM_S2S\behavior'];
        
        subfolder=[monkey '_S2S_saccades_' Target '_' Baseline];
        pdfpath=    [current_folder filesep subfolder];
        if ~exist(pdfpath,'dir')
            mkdir(current_folder,subfolder)
        end
        
        M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
        M2S_Analyze_settings.export_pdf=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
        M2S_Analyze_settings.append_pdf=0; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
        M2S_Analyze_settings.limit_convexity=0.2; % only relevant if concatinate_option==1
        M2S_Analyze_settings.type=6; %
        M2S_Analyze_settings.effector=0; %
        M2S_Analyze_settings.placeholders_visible=0; %
        M2S_Analyze_settings.choice=0; %
        
        Keys={'display',0,'runs_as_batches',0,'inferential_on',0,'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'downsampling',1};
        
        
        %% file selection
        monkey='Cornelius';
        xlspath=[dag_drive_IP 'Logs\Inactivation\Cornelius.xlsx'];
        Inputsequal={};
        Inputslist={};
        Inputsrange={};
        Inputsequal_for_batch={'Session'};
        combine_all_sheets_to_mastertable(xlspath,'Date');
        Inputsequal={'Target',Target};
        filelist_formatted = DAG_get_filelist_from_xls({'use_S2S_ina'},Inputsequal,Inputsrange,Inputslist, '', xlspath); %Inputslist
        ina_filelist=DAG_get_batch_input_from_filelist(filelist_formatted,Inputsequal_for_batch,[monkey '_ina']);
        
        
        Inputsequal={};
        filelist_formatted = DAG_get_filelist_from_xls({Baseline},Inputsequal,Inputsrange,Inputslist, '', xlspath); %Inputslist
        bas_filelist=DAG_get_batch_input_from_filelist(filelist_formatted,Inputsequal_for_batch,[monkey '_ina']);
        
        %[out_files tolerated_indexes]=match_closest_dates_batching(ina_filelist,bas_filelist,0);
        [out_files tolerated_indexes]=match_closest_dates_batching(ina_filelist,bas_filelist,1);
        ina_filelist=ina_filelist(tolerated_indexes);
        bas_filelist=out_files(tolerated_indexes);
        
        
        
        %% Summaries per Session and per group
        p=1;
        for g=1:2
            if g==1
                filelist=bas_filelist;
            elseif g==2
                filelist=ina_filelist;
            end
            for s=1:numel(filelist)
                sessions{s}=filelist{s}{1,1}(end-7:end);
            end
            
            %% meta-batch
            monkeypsych_analyze_input=[{vertcat(filelist{:})}; {Keys}];
            [MPA_out, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
            M2S_Analyze_settings.export_pdf_filename=[pdfpath filesep ' group ' num2str(g) ' meta_batch'];
            M2S_meta=S2S_session_analysis(MPA_out{1},M2S_Analyze_settings);
            S2S_session_plot(M2S_meta)
            
            [Met_ses(p,g),Met_norm(p,g),Met_tot(p,g)]=S2S_group_summary(M2S_meta,sessions);
            if M2S_Analyze_settings.export_pdf
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                
                export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'MetaGroupSummary2'], '-pdf', '-transparent') % pdf by run
                close(gcf);
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'MetaGroupSummary1'], '-pdf', '-transparent') % pdf by run
                close(gcf);
            end
            
            %% meta-batch
            
            monkeypsych_analyze_input=[filelist; repmat({Keys},1,numel(filelist))];
            [MPA_out, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
            
            for d=1:numel(MPA_out) %1:numel(subdir)  %subjects_with_naiverun{g} subjects_with_2cues{g}  ;insert here which subject(s) you want to analyze
                M2S_Analyze_settings.export_pdf_filename=[pdfpath filesep ' group ' num2str(g) ' ' sessions{d} ];
                M2S_output(g,d)=S2S_session_analysis(MPA_out{d},M2S_Analyze_settings);
                S2S_session_plot(M2S_output(g,d))
                
                [Ses_ses(d,g),Ses_norm(d,g),Ses_tot(d,g)]=S2S_group_summary(M2S_output(g,d),sessions);
                close(gcf);
                close(gcf);
            end
            
            [Sum_ses(p,g),Sum_norm(p,g),Sum_tot(p,g)]=S2S_group_summary(M2S_output(g,:),sessions);
            if M2S_Analyze_settings.export_pdf
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                
                export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'GroupSummary2'], '-pdf', '-transparent') % pdf by run
                close(gcf);
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'GroupSummary1'], '-pdf', '-transparent') % pdf by run
                close(gcf);
            end
        end
        
        
        
        %% per_session group comparisson
        comparison='Inactivation';
        stats_test='paired_ttest';
        for d=1:size(Ses_ses,1)
            date=ina_filelist{d}{1}(end-7:end);
            S2S_group_comparison(Ses_ses(d,:),Ses_norm(d,:),Ses_tot(d,:),stats_test,comparison,sessions)
            if M2S_Analyze_settings.export_pdf %&& numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                export_fig([pdfpath filesep date ' effects 3'], '-pdf', '-transparent') % pdf by run
                close(gcf);
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                export_fig([pdfpath filesep date ' effects 2'], '-pdf', '-transparent') % pdf by run
                close(gcf);
                wanted_size=[50 30];
                set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                export_fig([pdfpath filesep date ' effects 1'], '-pdf', '-transparent') % pdf by run
                close(gcf);
            end
        end
        
        %% meta batch group comparisson
        S2S_group_comparison(Met_ses(p,:),Met_norm(p,:),Met_tot(p,:),stats_test,comparison,sessions)
        if M2S_Analyze_settings.export_pdf %&& numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep 'meta effects 3'], '-pdf', '-transparent') % pdf by run
            close(gcf);
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep 'meta effects 2'], '-pdf', '-transparent') % pdf by run
            close(gcf);
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep 'meta effects 1'], '-pdf', '-transparent') % pdf by run
            close(gcf);
        end
        
        %% Group comparison
        comparison='Inactivation';
        stats_test='paired_ttest';
        S2S_group_comparison(Sum_ses,Sum_norm,Sum_tot,stats_test,comparison,sessions)
        if M2S_Analyze_settings.export_pdf %&& numel(Sum_ses(1).Hit.SL) == numel(Sum_ses(2).Hit.SL)
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep ' ' 'comparison3'], '-pdf', '-transparent') % pdf by run
            close(gcf);
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep ' ' 'comparison2'], '-pdf', '-transparent') % pdf by run
            close(gcf);
            wanted_size=[50 30];
            set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
            export_fig([pdfpath filesep ' ' 'comparison1'], '-pdf', '-transparent') % pdf by run
            close(gcf);
        end
    end
end


% %% MEMORY
%
% %% file selection
% current_folder=[dag_drive_IP 'Projects\STS_memory_saccades\behavior'];
% xlspath=[dag_drive_IP 'Logs\Inactivation\Cornelius.xlsx'];
% Inputsequal={};
% Inputslist={};
% Inputsrange={};
% Inputsequal_for_batch={'Session'};
% combine_all_sheets_to_mastertable(xlspath,'Date');
%
% Keys={'display',0,'runs_as_batches',0,'inferential_on',0,'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',0,'saccade_definition',4,'downsampling',1};
%
%
% filelist_formatted = DAG_get_filelist_from_xls({'use_mem_ina'},Inputsequal,Inputsrange,Inputslist, '', xlspath); %Inputslist
% ina_filelist=DAG_get_batch_input_from_filelist(filelist_formatted,Inputsequal_for_batch,[monkey '_ina']);
%
% filelist_formatted = DAG_get_filelist_from_xls({'use_mem_bas'},Inputsequal,Inputsrange,Inputslist, '', xlspath); %Inputslist
% bas_filelist=DAG_get_batch_input_from_filelist(filelist_formatted,Inputsequal_for_batch,[monkey '_ina']);
%
% ina_filelist=ina_filelist(end)
% bas_filelist=bas_filelist(end-2)
%
%
% %% Summaries per Session and per group
%
% pdffilename_memory_saccades=             [current_folder filesep monkey '_mem'];
% p=1;
% for g=1:2
%     if g==1
%         filelist=bas_filelist;
%     elseif g==2
%         filelist=ina_filelist;
%     end
%     monkeypsych_analyze_input=[filelist; repmat({Keys},1,numel(filelist))];
%     for s=1:numel(filelist)
%         sessions{s}=filelist{s}{1,1}(end-7:end);
%     end
%     [MPA_out, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
%
%     A=1;
%     clear Mem_output
%     for d=1:numel(MPA_out) %1:numel(subdir)  %subjects_with_naiverun{g} subjects_with_2cues{g}  ;insert here which subject(s) you want to analyze
%         [Mem_output(d) Info(d)]=Simple_memory_analysis(MPA_out{d});
%         wanted_size=[50 30];
%         set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%         export_fig([pdffilename_memory_saccades filesep ' group ' num2str(g) ' ' 'session '  num2str(d) ], '-pdf', '-transparent') % pdf by run
%         close(gcf);
%     end
%     Sum_ses(p,g).mem=Mem_output;
%     Sum_ses(p,g).meminf=Info;
%     %     [Sum_ses(p,g),Sum_norm(p,g),Sum_tot(p,g)]=S2S_group_summary(M2S_output,sessions);
%     %     if M2S_Analyze_settings.export_pdf
%     %         wanted_size=[50 30];
%     %         set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%     %
%     %         export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'GroupSummary2'], '-pdf', '-transparent') % pdf by run
%     %         close(gcf);
%     %         wanted_size=[50 30];
%     %         set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%     %         export_fig([pdfpath filesep ' group ' num2str(g) ' ' 'GroupSummary1'], '-pdf', '-transparent') % pdf by run
%     %         close(gcf);
%     %     end
% end
%
% Mem=vertcat(Sum_ses.mem);
% Info=vertcat(Sum_ses.meminf);
% for m=1:size(Mem,2)
%     Simple_memory_group_comparison(Mem(:,m),Info(:,m))
%     wanted_size=[50 30];
%     set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%     export_fig([pdffilename_memory_saccades filesep ' session ' num2str(m)], '-pdf', '-transparent') % pdf by run
%     close(gcf);
% end
%
%
% % global analysis_parameters
% % analysis_parameters.folders.extended_data = [dag_drive_IP 'Data\' monkey filesep];
% % analysis_parameters.files.updated_parameters = [dag_drive_IP 'protocols\' monkey filesep monkey '_protocol'];
% %
% % Pathpart_c=[dag_drive_IP 'Data\' monkey filesep];
% % Pathpart_i=[dag_drive_IP 'Data\' monkey '_ina' filesep];
% %
% %
% % batch_filelist_formatted{1}={[Pathpart_c '20170118'],4}; %[Pathpart_c '20170118'],3;                                                     %control previous day
% % batch_filelist_formatted{2}={[Pathpart_i '20170119'],4};                                                                                %run with largest potential effects
% % batch_filelist_formatted{3}={[Pathpart_i '20170119'],5;[Pathpart_i '20170119'],6;[Pathpart_i '20170119'],7;[Pathpart_i '20170119'],8};  %other inactivation runs
% % batch_filelist_formatted{4}={[Pathpart_i '20170119'],1;[Pathpart_i '20170119'],2;[Pathpart_i '20170119'],3};                            %potential same day baseline
% % batch_filelist_formatted{5}={[Pathpart_c '20170118'],2};                                                                                %control previous day - slightly different parameters
% %
% % batch_filelist_formatted{1}={[Pathpart_c '20170118'],4;[Pathpart_c '20170118'],2}; %[Pathpart_c '20170118'],3;                                                     %control previous day
% % batch_filelist_formatted{2}={[Pathpart_i '20170119'],5;[Pathpart_i '20170119'],6;[Pathpart_i '20170119'],7;[Pathpart_i '20170119'],8};  %other inactivation runs
% % batch_filelist_formatted{3}={[Pathpart_c '20170228'],7;[Pathpart_c '20170228'],4}; %[Pathpart_c '20170118'],3;                                                     %control previous day
% % batch_filelist_formatted{4}={[Pathpart_i '20170301'],5;[Pathpart_i '20170301'],6;[Pathpart_i '20170301'],8};  %other inactivation runs
% % batch_filelist_formatted{5}={[Pathpart_c '20170328'],2;[Pathpart_c '20170328'],3;[Pathpart_c '20170328'],4};                                                      %control previous day
% % batch_filelist_formatted{6}={[Pathpart_i '20170329'],6;[Pathpart_i '20170329'],8;};
% % batch_filelist_formatted{1}={[Pathpart_c '20170406'],2;[Pathpart_c '20170406'],3;[Pathpart_c '20170406'],5}; %[Pathpart_c '20170118'],3;                                                     %control previous day
% % batch_filelist_formatted{2}={[Pathpart_i '20170407'],5;[Pathpart_i '20170407'],6;[Pathpart_i '20170407'],7};  %other inactivation runs
% %
% %
% % for g=1:numel(out_comp)
% %     %colors={out_comp{k}.saccades.col_dim};
% %     %placeholders_visible=cellfun(@(x) any(x(:)>0), colors); % ?
% %
% %     M2S_Analyze_settings.placeholders_visible   =0;
% %     M2S_Analyze_settings.effector               =0;
% %     M2S_Analyze_settings.type                   =6;
% %     M2S_Analyze_settings.choice                 =0;
% %     M2S_Analyze_settings.export_pdf_filename    =[pdfpath '_batch_' num2str(g)];
% %
% %     M2S_output=S2S_session_analysis(out_comp{g},M2S_Analyze_settings);
% %     S2S_session_plot(M2S_output)
% %
% %
% %     if g<3
% %         subjects{p}={'Cor'};
% %     [Sum_ses(p,g),Sum_norm(p,g),Sum_tot(p,g)]=S2S_group_summary(M2S_output,subjects{p});
% %     if M2S_Analyze_settings.export_pdf
% %         export_fig([pdfpath '_batch_' num2str(g) ' ' 'GroupSummary2'], '-pdf', '-transparent') % pdf by run
% %         close(gcf);
% %         export_fig([pdfpath '_batch_' num2str(g) ' ' 'GroupSummary1'], '-pdf', '-transparent') % pdf by run
% %         close(gcf);
% %     end
% %     end
% % end
% %
% %
% %
% %

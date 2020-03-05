
monkey='Cornelius';
dag_drive_IP=DAG_get_server_IP;
current_folder=[dag_drive_IP 'Projects\STS_S2S_saccades\behavior'];

M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.export_pdf=0; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.append_pdf=0; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.limit_convexity=0.2; % only relevant if concatinate_option==1
M2S_Analyze_settings.type=6; %
M2S_Analyze_settings.effector=0; %
M2S_Analyze_settings.placeholders_visible=0; %
M2S_Analyze_settings.choice=0; %
% pdffilename_nonmasked_saccades=             [current_folder filesep monkey '_S2S_saccades_nonmasked'];
% pdffilename_nonmasked_reaches=              [current_folder filesep monkey '_S2S_reaches_nonmasked'];
% pdffilename_masked_saccades=                [current_folder filesep monkey '_S2S_saccades_masked'];
% pdffilename_masked_reaches=                 [current_folder filesep monkey '_S2S_reaches_masked'];
pdffilename_invisible_saccades=             [current_folder filesep monkey '_S2S_saccades_no_placeholders'];
% pdffilename_masked_saccades_double_cue=     [current_folder filesep monkey '_S2S_saccades_masked_two_cues'];
% pdffilename_nonmasked_saccades_double_cue=  [current_folder filesep monkey '_S2S_saccades_nonmasked_two_cues'];
% pdffilename_invisible_saccades_double_cue=  [current_folder filesep monkey '_S2S_saccades_no_placeholders_two_cues'];

Keys={'display',0,'runs_as_batches',0,'inferential_on',0,'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'downsampling',1};

datapath=[dag_drive_IP 'Data\Cornelius\'];


        M2S_Analyze_settings.export_pdf_filename=[pdffilename_invisible_saccades filesep ' gggg'];
        
        
        
        
filelist_formatted = {[datapath '20170608'],[4]};
    [MPA_out, ~, ~]= monkeypsych_analyze_working(filelist_formatted,Keys);

M2S_output=S2S_session_analysis(MPA_out{1},M2S_Analyze_settings);
        S2S_session_plot(M2S_output)
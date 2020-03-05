clear

filelist={{'C:\Users\lschneider\Desktop\tuebingen\20170510',[2]}};
M2S_Analyze_settings.export_pdf_filename=['C:\Users\lschneider\Desktop\tuebingen\plots\results_' ];
M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.export_pdf=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.append_pdf=0; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.limit_convexity=0.2; % only relevant if concatinate_option==1
M2S_Analyze_settings.type=6; %6
M2S_Analyze_settings.effector=0; %
M2S_Analyze_settings.placeholders_visible=0; %
M2S_Analyze_settings.choice=0; %

Keys={'display',0,'runs_as_batches',0,'inferential_on',0,'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'downsampling',1};

monkeypsych_analyze_input=[filelist; repmat({Keys},1,numel(filelist))];
[MPA_out, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
M2S_output=S2S_session_analysis(MPA_out{1},M2S_Analyze_settings);
S2S_session_plot(M2S_output);
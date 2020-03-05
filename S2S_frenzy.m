function S2S_frenzy(varargin)
%DAG_protocol_update('Cornelius',[20150810 20150814]);
%S2S_frenzy([20170118 20170118],'Cornelius')
if nargin==0
    %dates=[20170301 20170803];
    dates=[20150501 20151201];
    monkey='Cornelius';
elseif nargin==1
    dates=varargin{1};
else
    dates=varargin{1};
    monkey=varargin{2};
end
dag_drive_IP=DAG_get_server_IP;
start_date=dates(1);
%opengl('save', 'software') 

% DAG_plot_performance_for_dates('Cornelius',[20140701 20151130],[20140702 20140809 20140902 20141001 20141218 20150128 20150325 20150428 20150708 20151016 20151124], {'2 targets','4 targets', 'Reaches', 'Placeholders', 'Basic tasks rep', 'Shifting Offsets', 'More targets', '7 targets', 'Unilateral shapes', 'Exploration in dark', 'bilateral cues'})
% export_fig([monkey '_performance_across_days'], '-pdf','-transparent') % append to existing pdf
% close gcf

while true
    month=floor((start_date-floor(start_date/10000)*10000)/100);
    if month == 12
        end_date=(floor(start_date/10000)+1)*10000+100;
    else
        end_date=(floor(start_date/100)+1)*100;
    end
    if end_date>dates(2)
        S2S_frenzy_per_month(dag_drive_IP,monkey,[start_date dates(2)])
        break
    end
    S2S_frenzy_per_month(dag_drive_IP,monkey,[start_date end_date])
    start_date=end_date;
end

function S2S_frenzy_per_month(dag_drive_IP,monkey,dates)
%clear all
min_trials=50;
load_data=0;
save_data=0;
current_folder=[dag_drive_IP 'Projects\STS_PUL_MEM_S2S\baseline_behavior'];
folder_with_session_days=[dag_drive_IP 'Data\' monkey];
pdffilename_nonmasked_saccades=             [current_folder filesep monkey '_S2S_saccades_nonmasked'];
pdffilename_nonmasked_reaches=              [current_folder filesep monkey '_S2S_reaches_nonmasked'];
pdffilename_masked_saccades=                [current_folder filesep monkey '_S2S_saccades_masked'];
pdffilename_masked_reaches=                 [current_folder filesep monkey '_S2S_reaches_masked'];
pdffilename_invisible_saccades=             [current_folder filesep monkey '_S2S_saccades_no_placeholders'];
pdffilename_masked_saccades_double_cue=     [current_folder filesep monkey '_S2S_saccades_masked_two_cues'];
pdffilename_nonmasked_saccades_double_cue=  [current_folder filesep monkey '_S2S_saccades_nonmasked_two_cues'];
pdffilename_invisible_saccades_double_cue=  [current_folder filesep monkey '_S2S_saccades_no_placeholders_two_cues'];

%M2S_Analyze_settings.export_pdf=1;
%M2S_Analyze_settings.plot_only_LR=0; % for 7 targets task
M2S_Analyze_settings.concatinate_option=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.export_pdf=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.append_pdf=1; % 0=no concatination, 1=twosided and onesided differences -> put when analyzing more than one run
M2S_Analyze_settings.limit_convexity=0.2; % only relevant if concatinate_option==1


global analysis_parameters
analysis_parameters.folders.extended_data = [dag_drive_IP 'Data\' monkey filesep];
analysis_parameters.files.updated_parameters = [dag_drive_IP 'protocols\' monkey filesep monkey '_protocol'];

Keys={'display',0,'runs_as_batches',0,'inferential_on',0,'display',0,'nsacc_max',3,'summary',0,'keep_raw_data',1,'saccade_definition',4,'downsampling',1};
Inputsequal={'Session'};
[~, run_by_run_filelist_formatted] = DAG_get_filelist_from_folder(folder_with_session_days,dates);
batch_filelist_formatted = DAG_get_batch_input_from_xls(run_by_run_filelist_formatted,Inputsequal,Keys,dag_drive_IP,monkey);
monkeypsych_analyze_input=[batch_filelist_formatted; repmat({Keys},1,numel(batch_filelist_formatted))];

if load_data
    load([current_folder filesep 'out_M2S_Corny_' num2str(dates(1)) '-' num2str(dates(2)) '.mat']);
else
    [out_comp, ~, ~]= monkeypsych_analyze_working(monkeypsych_analyze_input{:});
    if save_data
        save([current_folder filesep 'out_M2S_Corny_' num2str(dates(1)) '-' num2str(dates(2)) '.mat'], 'out_comp')
    end
end


for k=1:numel(out_comp)    
    colors={out_comp{k}.saccades.col_dim};
    placeholders_visible=cellfun(@(x) any(x(:)>0), colors); % ?
    
    
    options(1).placeholders_visible=0;
    options(1).effector=0;
    options(1).type=6;
    options(1).choice=0;
    options(1).export_pdf_filename=pdffilename_invisible_saccades;    
    
    options(2).placeholders_visible=1;
    options(2).effector=0;
    options(2).type=6;
    options(2).choice=0;
    options(2).export_pdf_filename=pdffilename_masked_saccades;   
    
    options(3).placeholders_visible=1;
    options(3).effector=0;
    options(3).type=5;
    options(3).choice=0;
    options(3).export_pdf_filename=pdffilename_nonmasked_saccades;  
    
    options(4).placeholders_visible=1;
    options(4).effector=1;
    options(4).type=6;
    options(4).choice=0;
    options(4).export_pdf_filename=pdffilename_masked_reaches; 
    
    options(5).placeholders_visible=1;
    options(5).effector=1;
    options(5).type=5;
    options(5).choice=0;
    options(5).export_pdf_filename=pdffilename_nonmasked_reaches;
        
    options(6).placeholders_visible=0;
    options(6).effector=0;
    options(6).type=6;
    options(6).choice=1;
    options(6).export_pdf_filename=pdffilename_invisible_saccades_double_cue;    
    
    options(7).placeholders_visible=1;
    options(7).effector=0;
    options(7).type=6;
    options(7).choice=1;
    options(7).export_pdf_filename=pdffilename_masked_saccades_double_cue;    
    
    options(8).placeholders_visible=1;
    options(8).effector=0;
    options(8).type=5;
    options(8).choice=1;
    options(8).export_pdf_filename=pdffilename_nonmasked_saccades_double_cue;    
        
    for m=1:numel(options)
        if sum(placeholders_visible==options(m).placeholders_visible &  [out_comp{k}.binary.completed] & [out_comp{k}.task.type]==options(m).type  & [out_comp{k}.task.effector]==options(m).effector & [out_comp{k}.binary.choice]==options(m).choice) >min_trials
           
            M2S_Analyze_settings.placeholders_visible=options(m).placeholders_visible;
            M2S_Analyze_settings.effector=options(m).effector;
            M2S_Analyze_settings.type=options(m).type;
            M2S_Analyze_settings.choice=options(m).choice;
            M2S_Analyze_settings.export_pdf_filename=[options(m).export_pdf_filename num2str(out_comp{k}.selected(1).session)];
            
            M2S_output=S2S_session_analysis(out_comp{k},M2S_Analyze_settings);
            S2S_session_plot(M2S_output)
        end
    end
end
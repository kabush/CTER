%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.fmri_mod_exfb_beta]);
    eval(['! rm -rf ',proj.path.betas.fmri_mod_exfb_beta]);
    disp(['Creating ',proj.path.betas.fmri_mod_exfb_beta]);
    eval(['! mkdir ',proj.path.betas.fmri_mod_exfb_beta]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calc. fMRI betas (FB Modulate) of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{2};    % modulate
Nscans = proj.param.mri.Nscans(2); % # scans
Nvol = proj.param.mri.Nvol(2);
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;

%% ----------------------------------------
%% Fit beta series for each subject
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% ----------------------------------------
    %% Construct onset times for regression (combined)
    onset = [];
    trial_type = [];
    feedback = [];
    avg_feedback = [];
    valence = [];
    arousal = [];
    for j = 1:Nscans

        %% Load bids events files for the task
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);

        %% Search through events for event type
        for k = 1:numel(events.onset)

            type = events.trial_type(k,:);

            %% Match event type
            if(strcmp(type,'ex_stim') | strcmp(type,'fb_stim') | strcmp(type,'em_stim'))

                %% Add in event to observe affect immediately previous
                %% to extrinsic generated stimulus
                if(strcmp(type,'ex_stim'))

                    type(1:numel('ex_prev')) = 'ex_prev';
                    trial_type = [trial_type;type];
                    trial_onset = events.onset(k)-proj.param.mod.fb_prev_t;
                    adj_trial_onset = trial_onset+(j-1)*Nvol*TR;
                    onset = [onset;adj_trial_onset];
                    
                    valence = [valence;nan]; 
                    arousal = [arousal;nan]; 
                    feedback = [feedback;nan];
                    avg_feedback = [avg_feedback;nan]; 
                    
                end


                %% Add in event to observe affect immediately previous
                %% to feedback generated stimulus
                if(strcmp(type,'fb_stim') | strcmp(type,'em_stim'))

                    type = events.trial_type(k,:);
                    
                    if(strcmp(type,'fb_stim'))
                        type(1:numel('fb_prev')) = 'fb_prev';
                    else
                        type(1:numel('em_prev')) = 'em_prev';
                    end
                    
                    trial_type = [trial_type;type];
                    trial_onset = events.onset(k)-proj.param.mod.fb_prev_t;
                    adj_trial_onset = trial_onset+(j-1)*Nvol*TR;
                    onset = [onset;adj_trial_onset];
                    
                    valence = [valence;nan]; 
                    arousal = [arousal;nan]; 
                    feedback = [feedback;nan];
                    avg_feedback = [avg_feedback;nan]; 
                    
                end

                % Find onset time (of matched event type)
                trial_type = [trial_type;events.trial_type(k,:)];
                trial_onset = events.onset(k);

                % Correct the onset to adjust for concatenated scans
                adj_trial_onset = trial_onset+(j-1)*Nvol*TR;
                onset = [onset;adj_trial_onset];

                % Store associated labels
                valence = [valence;str2double(events.valence(k,:))];
                arousal = [arousal;str2double(events.arousal(k,:))];
                feedback = [feedback;str2double(events.feedback(k,:))];
                avg_feedback = [avg_feedback;str2double(events.avg_feedback(k,:))];
            end

        end

    end
        
    %% Write onset time file to tmp for 3dlss
    dlmwrite([proj.path.tmp,'stim_times.1D'],onset);

    %% ----------------------------------------
    %% Write data table for future MVPA and modeling
    mod_exfb_table = table(onset,...
                           trial_type,...
                           feedback,...
                           avg_feedback,...
                           valence,...
                           arousal);
    filename = [proj.path.betas.fmri_mod_exfb_beta,'sub-',name,'_task-',task,'_exfb_trials.tsv'];
    writetable(mod_exfb_table,filename,'FileType','text','Delimiter','\t');
        
    %% ----------------------------------------
    %% Construct combined censor file
    cmb_censor = [];
    for j = 1:Nscans
        filename = [proj.path.mri.mri_clean,'sub-',name,'/sub-',name,'_task-',task,num2str(j),'_censor.1D'];
        censor = load(filename);
        cmb_censor = [cmb_censor;censor];
    end

    %% Write onset time file to tmp
    dlmwrite([proj.path.tmp,'censor.1D'],cmb_censor);

    %% ----------------------------------------
    %% Construct combined motion file
    cmb_motion = [];
    for j = 1:Nscans
        filename = [proj.path.mri.mri_clean,'sub-',name,'/sub-',name,'_task-',task,num2str(j),'_motion.1D'];
        motion = load(filename);
        cmb_motion = [cmb_motion;motion];
    end

    %% Write onset time file to tmp
    dlmwrite([proj.path.tmp,'motion.1D'],cmb_motion);

    %% ----------------------------------------
    %% Write path to BOLD, censor, and motion files
    dlmwrite([proj.path.tmp,'in_path.txt'],proj.path.mri.mri_clean,'');
    dlmwrite([proj.path.tmp,'out_path.txt'],proj.path.betas.fmri_mod_exfb_beta,'');
    dlmwrite([proj.path.tmp,'subject.txt'],name,'');
    dlmwrite([proj.path.tmp,'task.txt'],task,'');
    dlmwrite([proj.path.tmp,'Nvol.txt'],Nvol,'');
    dlmwrite([proj.path.tmp,'TR.txt'],TR,'');
    dlmwrite([proj.path.tmp,'stim_t.txt'],stim_t,'');        

    %% ----------------------------------------
    %% Do the preprocessing
    eval(['! ',proj.path.code,'source/betas/fmri_3dlss ',...
          proj.path.tmp]);
    
    %% Clean-up
    eval(['! rm ',proj.path.tmp,'*']);
        
end

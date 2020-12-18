%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.fmri_id_in_beta]);
    eval(['! rm -rf ',proj.path.betas.fmri_id_in_beta]);
    disp(['Creating ',proj.path.betas.fmri_id_in_beta]);
    eval(['! mkdir ',proj.path.betas.fmri_id_in_beta]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Calc. fMRI beta-series (IN Identify) of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{1};    % identify
Nscans = proj.param.mri.Nscans(1); % # scans
Nvol = proj.param.mri.Nvol(1);
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
            if(strcmp(type,'in_stim'))

                % find onset time
                trial_type = [trial_type;events.trial_type(k,:)];
                trial_onset = events.onset(k);

                % correct the onset to adjust for concatenated scans
                adj_trial_onset = trial_onset+(j-1)*Nvol*TR;
                onset = [onset;adj_trial_onset];

                % store associated labels
                valence = [valence;str2double(events.valence(k,:))];
                arousal = [arousal;str2double(events.arousal(k,:))];

            end

            %% Match event type
            if(strcmp(type,'in_feel'))

                for m = 1:proj.param.in.Nfeel

                    % find onset time
                    trial_type = [trial_type;events.trial_type(k,:)];
                    trial_onset = events.onset(k);

                    % correct the onset to adjust for concatenated scans
                    adj_trial_onset = trial_onset+(m-1)*TR+((j-1)*Nvol*TR);
                    onset = [onset;adj_trial_onset];
                    
                    % store associated labels
                    valence = [valence;str2double(events.valence(k,:))];
                    arousal = [arousal;str2double(events.arousal(k,:))];

                end

            end

        end

    end
        
    %% Write onset time file to tmp for 3dlss
    dlmwrite([proj.path.tmp,'stim_times.1D'],onset);

    %% ----------------------------------------
    %% Write scores to in_betas for future MVPA
    id_in_table = table(onset,...
                        trial_type,...
                        valence,...
                        arousal);
    filename = [proj.path.betas.fmri_id_in_beta,'sub-',name,'_task-',task,'_in_trials.tsv'];
    writetable(id_in_table,filename,'FileType','text','Delimiter','\t');
    
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
    dlmwrite([proj.path.tmp,'out_path.txt'],proj.path.betas.fmri_id_in_beta,'');
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

%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for EMG
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.emg_id_ex_beta]);
    eval(['! rm -rf ',proj.path.betas.emg_id_ex_beta]);
    disp(['Creating ',proj.path.betas.emg_id_ex_beta]);
    eval(['! mkdir ',proj.path.betas.emg_id_ex_beta]);
end

%% Load subjects
subjs = load_subjs(proj);

%% Define the task and type details
task = proj.param.mri.tasks{1};    % identify
Nscans = proj.param.mri.Nscans(1); % # scans
n_vols = proj.param.mri.Nvol(1);
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;

%% Fit beta series for each subject
for i=1:numel(subjs)
    
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Initialize beta structure
    ex_betas = struct();

    %% Construct onset times for regression (combined)
    for j = 1:Nscans

        %% Load bids events files for the task
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);

        %% Search through events for event type
        onset = [];
        for k = 1:numel(events.onset)

            type = events.trial_type(k,:);

            %% Match event type
            if(strcmp(type,'ex_stim'))
                onset = [onset;events.onset(k)];
            end

        end
        
        ex_betas.(['id',num2str(j)]) = [];
        path = [proj.path.physio.emg_clean,'sub-',name,'_task-',task,'_run',num2str(j),'.mat'];
        load(path);

        for k = 1:numel(onset)
            start_time = onset(k);
            start_id = start_time*proj.param.physio.hz_emg;
            stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
            end_id = start_id+stim_samples-1;
            ex_beta = sum(emg(round(start_id):round(end_id)));
            ex_betas.(['id',num2str(j)]) = [ex_betas.(['id',num2str(j)]);ex_beta];
        end
        
        %% Normalize fits by run
        ex_betas.(['id',num2str(j)]) = zscore(ex_betas.(['id',num2str(j)]));
        
    end

    %% Save individual betas
    save([proj.path.betas.emg_id_ex_beta,'sub-',name,'_id_ex_betas.mat'],'ex_betas');
        
end
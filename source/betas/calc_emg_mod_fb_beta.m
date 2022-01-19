%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for EMG
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.emg_mod_fb_beta]);
    eval(['! rm -rf ',proj.path.betas.emg_mod_fb_beta]);
    disp(['Creating ',proj.path.betas.emg_mod_fb_beta]);
    eval(['! mkdir ',proj.path.betas.emg_mod_fb_beta]);
end

%% Create the subjects to be analyzed
subjs = load_subjs(proj);

%% Define the task and type details
task = proj.param.mri.tasks{2};    % identify
Nscans = proj.param.mri.Nscans(2); % # scans
n_vols = proj.param.mri.Nvol(2);
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;

%% Fit beta series for each subject
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Initialize emg beta structure
    fb_betas = struct();
    fb_onsets = struct();

    for j = 1:Nscans

        %% Load bids events files for the task
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);

        %% Construct onset times for regression (combined)
        onset = [];
        for k = 1:numel(events.onset)

            type = events.trial_type(k,:);

            %% Match event type
            if(strcmp(type,'fb     '))
                onset = [onset;events.onset(k)];
            end

        end

        fb_betas.(['mod',num2str(j)]) = [];
        path = [proj.path.physio.emg_clean,'sub-',name,'_task-',task,'_run',num2str(j),'.mat'];
        load(path);
        
        for k = 1:numel(onset)
            start_time = onset(k);
            start_mod = start_time*proj.param.physio.hz_emg;
            stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
            end_mod = start_mod+stim_samples-1;
            fb_beta = sum(emg(round(start_mod):round(end_mod)));
            fb_betas.(['mod',num2str(j)]) = [fb_betas.(['mod',num2str(j)]);fb_beta];

        end
        
        %% Normalize fits by run
        fb_betas.(['mod',num2str(j)]) = zscore(fb_betas.(['mod',num2str(j)]));

        %% Store onset times
        fb_onsets.(['mod',num2str(j)]) = onset;

    end

    %% Save individual betas
    save([proj.path.betas.emg_mod_fb_beta,'sub-',name,'_mod_fb_betas.mat'],'fb_betas');
    save([proj.path.betas.emg_mod_fb_beta,'sub-',name,'_mod_fb_onsets.mat'],'fb_onsets');        

end

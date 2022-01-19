%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for EMG
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.emg_mod_exfb_beta]);
    eval(['! rm -rf ',proj.path.betas.emg_mod_exfb_beta]);
    disp(['Creating ',proj.path.betas.emg_mod_exfb_beta]);
    eval(['! mkdir ',proj.path.betas.emg_mod_exfb_beta]);
end

%% Create the subjects to be analyzed
subjs = load_subjs(proj);

%% Define the task and type details
task = proj.param.mri.tasks{2};  
Nscans = proj.param.mri.Nscans(2); 
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
    betas = struct();
    trial_types = struct();
    valences = struct();
    arousals = struct();

    for j = 1:Nscans

        %% Load bids events files for the task
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);

        %% Construct onset times for regression (combined)
        onset = [];
        trial_type = [];
        valence = [];
        arousal = []; 

        for k = 1:numel(events.onset)

            type = events.trial_type(k,:);
            val = events.valence(k,:);
            aro = events.arousal(k,:);

            %% Match event type
            if(strcmp(type,'ex_stim'))
                onset = [onset;events.onset(k)];
                trial_type = [trial_type;0];
                valence = [valence;str2double(val)];
                arousal = [arousal;str2double(aro)];
            end

            %% Match event type
            if(strcmp(type,'fb_stim') | strcmp(type,'em_stim'))
                onset = [onset;events.onset(k)];
                trial_type = [trial_type;1];
                valence = [valence;str2double(val)];
                arousal = [arousal;str2double(aro)];
            end

        end

        betas.(['mod',num2str(j)]) = [];
        path = [proj.path.physio.emg_clean,'sub-',name,'_task-',task,'_run',num2str(j),'.mat'];
        load(path);
        
        for k = 1:numel(onset)
            start_time = onset(k);
            start_mod = start_time*proj.param.physio.hz_emg;
            stim_samples = proj.param.mri.TR*proj.param.physio.hz_emg;
            end_mod = start_mod+stim_samples-1;
            ex_beta = sum(emg(round(start_mod):round(end_mod)));
            betas.(['mod',num2str(j)]) = [betas.(['mod',num2str(j)]);ex_beta];

        end
        
        %% Normalize fits by run
        betas.(['mod',num2str(j)]) = zscore(betas.(['mod',num2str(j)]));

        %% Save trial types
        trial_types.(['mod',num2str(j)]) = [];
        trial_types.(['mod',num2str(j)]) = trial_type;

        %% Save affect information
        valences.(['mod',num2str(j)]) = [];
        valences.(['mod',num2str(j)]) = valence;

        arousals.(['mod',num2str(j)]) = [];
        arousals.(['mod',num2str(j)]) = arousal;

    end

    %% Save individual betas
    save([proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_betas.mat'],'betas');
    save([proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_trial_types.mat'],'trial_types');        
    save([proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_valences.mat'],'valences');        
    save([proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_arousals.mat'],'arousals');        

end

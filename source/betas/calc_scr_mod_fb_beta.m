%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for SCR
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.scr_mod_fb_beta]);
    eval(['! rm -rf ',proj.path.betas.scr_mod_fb_beta]);
    disp(['Creating ',proj.path.betas.scr_mod_fb_beta]);
    eval(['! mkdir ',proj.path.betas.scr_mod_fb_beta]);
end

%% Create the subjects to be analyzed
subjs = load_subjs(proj);

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{2};    % identify
Nscans = proj.param.mri.Nscans(2); % # scans
n_vols = proj.param.mri.Nvol(2);
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

    %% Initialize scr beta structure
    fb_betas = struct();
    fb_onsets = [];

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

        %% build design(s)
        [prime_fb other_fb] = scr_dsgn_preproc(proj,n_vols,onset);

        %% LSS of scr signal (Mumford, 2012) - Modulate
        fb_betas.(['mod',num2str(j)]) = [];
        path = [proj.path.physio.scr_clean,'sub-',name,'_task-',task,'_run',num2str(j),'.mat'];
        load(path);
    
        for k=1:size(prime_fb,1)
            mdl_fb = regstats(scr,[prime_fb(k,:)',other_fb(k,:)']);
            fb_betas.(['mod',num2str(j)]) = [fb_betas.(['mod',num2str(j)]),mdl_fb.beta(2)];
        end
        
        %% Normalize fits by run
        fb_betas.(['mod',num2str(j)]) = zscore(fb_betas.(['mod', ...
                            num2str(j)]));

        %% Store onset times
        fb_onsets.(['mod',num2str(j)]) = onset;
            
    end

    %% ----------------------------------------
    %% Save individual betas
    save([proj.path.betas.scr_mod_fb_beta,'sub-',name,'_mod_fb_betas.mat'],'fb_betas');
    save([proj.path.betas.scr_mod_fb_beta,'sub-',name,'_mod_fb_onsets.mat'],'fb_onsets');

end

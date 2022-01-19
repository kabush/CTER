%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for SCR
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.physio.emg_clean]);
    eval(['! rm -rf ',proj.path.physio.emg_clean]);
    disp(['Creating ',proj.path.physio.emg_clean]);
    eval(['! mkdir ',proj.path.physio.emg_clean]);
end

%% Create the subjects to be analyzed
subjs = load_subjs(proj);

logger(['************************************'],proj.path.logfile);
logger(['Processing EMG (corrugator) of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************'],proj.path.logfile);

for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Define input/outputs paths
    in_path = [proj.path.bids,'sub-',name,'/func/sub-',name];
    out_path = [proj.path.physio.emg_clean];

    %% ----------------------------------------
    %% Build Identify 1 EMG (corrugator)
    try
        n_vols = proj.param.mri.Nvol(1);  %Identify length
        path = [in_path,'_task-identify1_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        physio_emg_corr = physio_raw(:,proj.param.physio.chan_emg_corr);
        emg = abs(emg_preproc(proj,n_vols,physio_emg_corr));
        save([out_path,'sub-',name,'_task-identify_run1.mat'],'emg');        
    catch
        logger(['  -Processing Error: EMG of Identify run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Identify 2 EMG (corrugator)
    try
        n_vols = proj.param.mri.Nvol(1);  %Identify length
        path = [in_path,'_task-identify2_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_emg_corr = load(zipfile{1});
        emg = abs(emg_preproc(proj,n_vols,physio_emg_corr));
        save([out_path,'sub-',name,'_task-identify_run2.mat'],'emg');        
    catch
        logger(['  -Processing Error: EMG of Identify run 2: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Rest EMG (corrugator)
    try
        n_vols = proj.param.mri.Nvol(3);  %Rest length
        path = [in_path,'_task-rest_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_emg_corr = load(zipfile{1});
        emg = abs(emg_preproc(proj,n_vols,physio_emg_corr));
        save([out_path,'sub-',name,'_task-rest.mat'],'emg');        
    catch
        logger(['  -Processing Error: EMG of Rest: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Modulate 1 EMG (corrugator)
    try
        n_vols = proj.param.mri.Nvol(2);  %Modulate length
        path = [in_path,'_task-modulate1_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_emg_corr = load(zipfile{1});
        emg = abs(emg_preproc(proj,n_vols,physio_emg_corr));
        save([out_path,'sub-',name,'_task-modulate_run1.mat'],'emg');        
    catch
        logger(['  -Processing Error: EMG of Modulate run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Modulate 2 EMG (corrugator)
    try
        n_vols = proj.param.mri.Nvol(2);  %Modulate length
        path = [in_path,'_task-modulate2_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_emg_corr = load(zipfile{1});
        emg = abs(emg_preproc(proj,n_vols,physio_emg_corr));
        save([out_path,'sub-',name,'_task-modulate_run2.mat'],'emg');        
    catch
        logger(['  -Processing Error: EMG of Modulate run 2: ',path],proj.path.logfile);
    end

end

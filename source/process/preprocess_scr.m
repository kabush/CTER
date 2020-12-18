%% Load in path data
load('proj.mat');

%% Set-up Directory Structure for SCR
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.physio.scr_clean]);
    eval(['! rm -rf ',proj.path.physio.scr_clean]);
    disp(['Creating ',proj.path.physio.scr_clean]);
    eval(['! mkdir ',proj.path.physio.scr_clean]);
end

%% Create the subjects to be analyzed
subjs = load_subjs(proj);

logger(['************************************'],proj.path.logfile);
logger(['Processing SCR of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************'],proj.path.logfile);

for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Define input/outputs paths
    in_path = [proj.path.bids,'sub-',name,'/func/sub-',name];
    out_path = [proj.path.physio.scr_clean];

    %% ----------------------------------------
    %% Build Identify 1 SCR
    try
        n_vols = proj.param.mri.Nvol(1);  %Identify length
        path = [in_path,'_task-identify1_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        scr = scr_preproc(proj,n_vols,physio_raw);
        save([out_path,'sub-',name,'_task-identify_run1.mat'],'scr');        
    catch
        logger(['  -Processing Error: SCR of Identify run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Identify 2 SCR
    try
        n_vols = proj.param.mri.Nvol(1);  %Identify length
        path = [in_path,'_task-identify2_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        scr = scr_preproc(proj,n_vols,physio_raw);
        save([out_path,'sub-',name,'_task-identify_run2.mat'],'scr');        
    catch
        logger(['  -Processing Error: SCR of Identify run 2: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Rest SCR
    try
        n_vols = proj.param.mri.Nvol(3);  %Rest length
        path = [in_path,'_task-rest_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        scr = scr_preproc(proj,n_vols,physio_raw);
        save([out_path,'sub-',name,'_task-rest.mat'],'scr');        
    catch
        logger(['  -Processing Error: SCR of Rest: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Modulate 1 SCR
    try
        n_vols = proj.param.mri.Nvol(2);  %Modulate length
        path = [in_path,'_task-modulate1_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        scr = scr_preproc(proj,n_vols,physio_raw);
        save([out_path,'sub-',name,'_task-modulate_run1.mat'],'scr');        
    catch
        logger(['  -Processing Error: SCR of Modulate run 1: ',path],proj.path.logfile);
    end

    %% ----------------------------------------
    %% Build Modulate 2 SCR
    try
        n_vols = proj.param.mri.Nvol(2);  %Modulate length
        path = [in_path,'_task-modulate2_physio.tsv.gz']
        zipfile = gunzip(path);
        physio_raw = load(zipfile{1});
        scr = scr_preproc(proj,n_vols,physio_raw);
        save([out_path,'sub-',name,'_task-modulate_run2.mat'],'scr');        
    catch
        logger(['  -Processing Error: SCR of Modulate run 2: ',path],proj.path.logfile);
    end

end

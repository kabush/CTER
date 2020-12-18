%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Construct Surrogate IN trials from Rest Affect '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.rest_in_affect]);
    eval(['! rm -rf ',proj.path.analysis.rest_in_affect]);
    disp(['Creating ',proj.path.analysis.rest_in_affect]);
    eval(['! mkdir ',proj.path.analysis.rest_in_affect]);
end

%% ----------------------------------------
%% Load Subjects
subjs = load_subjs(proj);

%% ----------------------------------------
%% Extract each subject's surrogate IN 
for i = 1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % debug
    logger([subj_study,':',name],proj.path.logfile);

    %% VALENCE
    try
        in_surr = calc_surrogate_in(proj,name,'v');
        save([proj.path.analysis.rest_in_affect,'sub-',name,'_v_in_surr.mat'],'in_surr');
    catch
        logger(['  -Error in extracting surrogate Val IN trials from rest'],proj.path.logfile);
    end

    %% AROUSAL
    try
        in_surr = calc_surrogate_in(proj,name,'a');
        save([proj.path.analysis.rest_in_affect,'sub-',name,'_a_in_surr.mat'],'in_surr');
    catch
        logger(['  -Error in extracting surrogate Aro IN trials from rest'],proj.path.logfile);
    end

end

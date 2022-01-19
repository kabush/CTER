%% Load in path data
load('proj.mat');

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mri.gm_mask]);
    eval(['! rm -rf ',proj.path.mri.gm_mask]);
    disp(['Creating ',proj.path.mri.gm_mask]);
    eval(['! mkdir ',proj.path.mri.gm_mask]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);

logger(['************************************'],proj.path.logfile);
logger(['build gm masks of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************'],proj.path.logfile);

%% Postprocess fmriprep
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    disp(['sub-',name,]);

    %% write processing information to tmp for afni code
    dlmwrite([proj.path.tmp,'in_path.txt'],proj.path.fmriprep,'');
    dlmwrite([proj.path.tmp,'template_path.txt'],proj.path.mri.mri_clean,'');
    dlmwrite([proj.path.tmp,'out_path.txt'],proj.path.mri.gm_mask,'');
    dlmwrite([proj.path.tmp,'subject.txt'],name,'');
    dlmwrite([proj.path.tmp,'space.txt'],num2str(proj.param.mri.anat_space),'');
    dlmwrite([proj.path.tmp,'gm_prob.txt'],num2str(proj.param.mri.gm_prob),'');

    %% Do the preprocessing
    eval(['! ',proj.path.code,'source/process/gm_mask_process_afni ',...
          proj.path.tmp]);
    
    %% Clean-up
    eval(['! rm ',proj.path.tmp,'*']);

end

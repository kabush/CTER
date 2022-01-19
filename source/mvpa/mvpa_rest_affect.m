%% Load in path data
load('proj.mat');


%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.rest_affect]);
    eval(['! rm -rf ',proj.path.mvpa.rest_affect]);
    disp(['Creating ',proj.path.mvpa.rest_affect]);
    eval(['! mkdir ',proj.path.mvpa.rest_affect]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Model REST surrogate IN affect (Platt space) '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{3};    % identify
Nvol = proj.param.mri.Nvol(3);
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;
Npseudo = proj.param.rest.Npseudo;      % Num. pseudo stims per sample
Nresample = proj.param.rest.Nresample;  % Num. times to resample
Ntrans =  proj.param.rest.Ntrs_trans;   
Ntail =  proj.param.rest.Ntrs_tail;

%% ----------------------------------------
%% Quality check prior to modeling
cnt = 0;
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    try

        % Check if there is enough usable data (censor file)
        censor_path = [proj.path.mri.mri_clean,'sub-',name,'/sub-',name,'_task-',task,'_censor.1D'];
        censor = load(censor_path);
        sum_censor = sum(censor);
        Npseudo = proj.param.rest.Npseudo;
        Nmotion = proj.param.rest.Nmotion;

        if(sum_censor>(Npseudo+Nmotion))

            % Test if the models exist
            cv_svm_path = [proj.path.mvpa.id_ex_gm_mdl];
            cv_v_svm_name = ['sub-',name,'_v_model.mat'];
            cv_a_svm_name = ['sub-',name,'_a_model.mat'];
            v_model = load([cv_svm_path,cv_v_svm_name]);
            a_model = load([cv_svm_path,cv_a_svm_name]);
            cnt = cnt + 1;
            good_subjs{cnt} = subjs{i};
        
        else
            logger([' ***Not enough usable data: sub-',name],proj.path.logfile);
        end

    catch
        logger([' *Model error: sub-',name,],proj.path.logfile);
    end

end
subjs = good_subjs; %%Subjects are now correctly sized

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    
    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% ----------------------------------------
    %% Load gray matter mask 
    gm_nii = load_untouch_nii([proj.path.mri.gm_mask,'sub-',name,'_gm_mask.nii']);
    mask = double(gm_nii.img);
    brain_size=size(mask);
    mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
    in_brain=find(mask==1);  
    
    %% ----------------------------------------
    %% Set-up paths/names
    tmp_path = proj.path.tmp;
    rest_path = [proj.path.mri.mri_clean,'sub-',name,'/'];
    censor_path = [rest_path,'sub-',name,'_task-',task,'_censor.1D'];
    motion_path = [rest_path,'sub-',name,'_task-',task,'_motion.1D'];
    censor = load(censor_path);
    motion = load(motion_path);
    
    %% ----------------------------------------
    %% Write path to BOLD, censor, and motion files
    dlmwrite([proj.path.tmp,'in_path.txt'],proj.path.mri.mri_clean,'');
    dlmwrite([proj.path.tmp,'subject.txt'],name,'');
    dlmwrite([proj.path.tmp,'task.txt'],task,'');
    dlmwrite([proj.path.tmp,'Nvol.txt'],Nvol,'');
    dlmwrite([proj.path.tmp,'TR.txt'],TR,'');
    dlmwrite([proj.path.tmp,'stim_t.txt'],stim_t,'');        
    dlmwrite([proj.path.tmp,'motion.1D'],motion);
    dlmwrite([proj.path.tmp,'censor.1D'],censor);

    %% ----------------------------------------
    %% Iteration based data from regression
    id_rest = zeros(Nvol,Nresample);
    
    %% Will become avg values (used to keep running sum also)
    mu_v_rest = zeros(Nvol,1);
    mu_a_rest = zeros(Nvol,1);

    %% Estimate beta and affect multiple times
    for j=1:Nresample
        
        %% Sample a stimulus timing set (50% of total TRs)
        stim_ids = randsample((Ntrans+1):(Nvol-Ntail),Npseudo);
        stim_times = stim_ids*TR;
        dlmwrite([proj.path.tmp,'stim_times.1D'],stim_times','');

        %% Execute beta-series regression
        eval(['! ',proj.path.code,'source/mvpa/fmri_3dlss ',...
              proj.path.tmp]);

        %% Load beta-series
        base_nii = load_untouch_nii([tmp_path,'sub-',name,'_lss.nii']);
        brain_size = size(base_nii.img);

        %% *** Post clean-up tmp ***
        eval(['! rm ',tmp_path,'sub-',name,'_lss.nii']);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));

        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';

        %% ----------------------------------------
        %% VALENCE Machine Learning Prediction
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_model.mat']);
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_rdx_ids.mat']);
        [label,hd]= predict(model,zscore(subj_img(:,rdx_ids)));
        v_label = label;
        v_hd = hd(:,2);
        v_platt = 1./(1+exp(-v_hd));

        %% ----------------------------------------
        %% AROUSAL Machine Learning Prediction
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_model.mat']);
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_rdx_ids.mat']);
        [label,hd]= predict(model,zscore(subj_img(:,rdx_ids)));
        a_label = label;
        a_hd = hd(:,2);
        a_platt = 1./(1+exp(-a_hd));

        %% ----------------------------------------
        %% STORE Predictions

        %% store the ids involved in these pseudo stims
        id_rest(stim_ids,j) = ones(numel(stim_ids),1);
        
        %% build total affect scores (via running sums)
        mu_v_rest(stim_ids) = mu_v_rest(stim_ids)+v_platt;
        mu_a_rest(stim_ids) = mu_a_rest(stim_ids)+a_platt;

    end

    %% *** Post clean-up tmp ***
    eval(['! rm ',tmp_path,'*']);

    %% Use id_count to compute averages scores|states
    id_count = sum(id_rest,2);
    good_ids = find(id_count>0);
    for j=1:numel(good_ids)
        id = good_ids(j);
        cnt = id_count(id);
        mu_v_rest(id) = mu_v_rest(id)/cnt;
        mu_a_rest(id) = mu_a_rest(id)/cnt;
    end

    %% ----------------------------------------
    %% SAVE OUT
    dlmwrite([proj.path.mvpa.rest_affect,'sub-',name,'_v_rest.1D'],mu_v_rest,'');    
    dlmwrite([proj.path.mvpa.rest_affect,'sub-',name,'_a_rest.1D'],mu_a_rest,'');    

end

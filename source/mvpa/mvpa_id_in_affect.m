%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Prediction Affect of Identify (IN) Trials  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.id_in_affect]);
    eval(['! rm -rf ',proj.path.mvpa.id_in_affect]);
    disp(['Creating ',proj.path.mvpa.id_in_affect]);
    eval(['! mkdir ',proj.path.mvpa.id_in_affect]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load gray matter mask 
        gm_nii = load_untouch_nii([proj.path.mri.gm_mask,'sub-',name,'_gm_mask.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        base_nii = load_untouch_nii([proj.path.betas.fmri_id_in_beta,'sub-',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));

        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';

        %% Load labels
        filename = [proj.path.betas.fmri_id_in_beta,'sub-',name,'_task-identify_in_trials.tsv'];
        events = tdfread(filename);

        %% ----------------------------------------
        %% VALENCE Machine Learning Prediction
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_model.mat']);
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_rdx_ids.mat']);
        [label,hd]= predict(model,zscore(subj_img(:,rdx_ids)));
        v_label = label;
        v_hd = hd(:,2);
        pred_valence = 1./(1+exp(-v_hd)); % Platt scale

        %% ----------------------------------------
        %% AROUSAL Machine Learning Prediction
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_model.mat']);
        load([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_rdx_ids.mat']);
        [label,hd]= predict(model,zscore(subj_img(:,rdx_ids)));
        a_label = label;
        a_hd = hd(:,2);
        pred_arousal = 1./(1+exp(-a_hd)); % Platt scale

        %% ----------------------------------------
        %% Create data table and save
        onset = events.onset;
        trial_type = events.trial_type;
        valence = 1./(1+exp(-(events.valence-proj.param.mvpa.likert))); % Ctr. & Platt scale
        arousal = 1./(1+exp(-(events.arousal-proj.param.mvpa.likert))); % Ctr. & Platt scale

        id_in_table = table(onset,...
                            trial_type,...
                            valence,...
                            arousal,...
                            pred_valence,...
                            pred_arousal);
        filename = [proj.path.mvpa.id_in_affect,'sub-',name,'_task-identify_in_predictions.tsv'];
        writetable(id_in_table,filename,'FileType','text','Delimiter','\t');

    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end

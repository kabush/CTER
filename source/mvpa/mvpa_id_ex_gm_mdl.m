%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Intra-subject LOOCV MVPA of Gray Matter Models  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.id_ex_gm_mdl]);
    eval(['! rm -rf ',proj.path.mvpa.id_ex_gm_mdl]);
    disp(['Creating ',proj.path.mvpa.id_ex_gm_mdl]);
    eval(['! mkdir ',proj.path.mvpa.id_ex_gm_mdl]);
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
        base_nii = load_untouch_nii([proj.path.betas.fmri_id_ex_beta,'sub-',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));

        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';

        %% Quality control
        % Find beta maps with nans (likely motion problem)
        mu_col = mean(abs(subj_img),2);
        rmv_beta_ids = find(isnan(mu_col));

        % Report problem
        disp(['   # NaN betas=',num2str(numel(rmv_beta_ids))]);
        
        % Find good beta maps
        keep_beta_ids = find(~isnan(mu_col));
        subj_img = subj_img(keep_beta_ids,:);

        disp(['   # kept betas=',num2str(numel(keep_beta_ids))]);

        %% Load labels
        filename = [proj.path.betas.fmri_id_ex_beta,'sub-',name,'_task-identify_ex_trials.tsv'];
        events = tdfread(filename);

        %% ----------------------------------------
        %% VALENCE Machine Learning

        %% Load/format labels
        scores = events.valence(keep_beta_ids);
        labels = 0*scores;
        labels(find(scores>=proj.param.mvpa.likert)) = proj.param.mvpa.pos_cls;
        labels(find(scores<proj.param.mvpa.likert)) = proj.param.mvpa.neg_cls;
        
        %% Balance pos/neg examples for one-time fit
        pos_ids = find(labels==proj.param.mvpa.pos_cls);
        neg_ids = find(labels==proj.param.mvpa.neg_cls);
        Npos = numel(pos_ids);
        Nneg = numel(neg_ids);
        if(Npos >= Nneg)
            Nsample = Nneg;
        else
            Nsample = Npos;
        end
            
        %% Randomly re-order and combined samples
        rnd_pos_ids = pos_ids(randsample(1:numel(pos_ids),Nsample));
        rnd_neg_ids = neg_ids(randsample(1:numel(neg_ids),Nsample));
        trn_ids = [rnd_pos_ids;rnd_neg_ids];
        
        %% Find voxels with nonzero beta values
        mu = mean(abs(subj_img),1);
        rdx_ids = find(mu>0); %%GM ids to be used in MVPA fit

        %% Z-score reducted image
        trn_img = zscore(subj_img(trn_ids,rdx_ids));
        trn_label = labels(trn_ids);

        disp(['   Val: # trn labels=',num2str(numel(trn_label))]);

        %% Fit classifier
        model = fitcsvm(trn_img,trn_label,'KernelFunction',proj.param.mvpa.kernel);

        %% Save out model and GM ids
        save([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_model.mat'],'model');
        save([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_v_rdx_ids.mat'],'rdx_ids');

        %% ----------------------------------------
        %% AROUSAL Machine Learning

        %% Load/format labels
        scores = events.arousal(keep_beta_ids);
        labels = 0*scores;
        labels(find(scores>=proj.param.mvpa.likert)) = proj.param.mvpa.pos_cls;
        labels(find(scores<proj.param.mvpa.likert)) = proj.param.mvpa.neg_cls;
        
        %% Balance pos/neg examples for one-time fit
        pos_ids = find(labels==proj.param.mvpa.pos_cls);
        neg_ids = find(labels==proj.param.mvpa.neg_cls);
        Npos = numel(pos_ids);
        Nneg = numel(neg_ids);
        if(Npos >= Nneg)
            Nsample = Nneg;
        else
            Nsample = Npos;
        end
            
        %% Randomly re-order and combined samples
        rnd_pos_ids = pos_ids(randsample(1:numel(pos_ids),Nsample));
        rnd_neg_ids = neg_ids(randsample(1:numel(neg_ids),Nsample));
        trn_ids = [rnd_pos_ids;rnd_neg_ids];
        
        %% Find voxels with nonzero beta values
        mu = mean(abs(subj_img),1);
        rdx_ids = find(mu>0); %%GM ids to be used in MVPA fit

        %% Z-score reducted image
        trn_img = zscore(subj_img(trn_ids,rdx_ids));
        trn_label = labels(trn_ids);

        disp(['   Aro: # trn labels=',num2str(numel(trn_label))]);

        %% Fit classifier
        model = fitcsvm(trn_img,trn_label,'KernelFunction',proj.param.mvpa.kernel);

        %% Save out model and GM ids
        save([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_model.mat'],'model');
        save([proj.path.mvpa.id_ex_gm_mdl,'sub-',name,'_a_rdx_ids.mat'],'rdx_ids');
        
    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end
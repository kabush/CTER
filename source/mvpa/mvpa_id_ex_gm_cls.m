%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Intra-subject LOOCV MVPA of Gray Matter Features'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.id_ex_gm_cls]);
    eval(['! rm -rf ',proj.path.mvpa.id_ex_gm_cls]);
    disp(['Creating ',proj.path.mvpa.id_ex_gm_cls]);
    eval(['! mkdir ',proj.path.mvpa.id_ex_gm_cls]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% allocate storage (logging)
all_v_cls_acc = [];
all_a_cls_acc = [];

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

        % Quality control
        
        % Find betamaps with nans (likely motion problem)
        mu_col = mean(abs(subj_img),2);
        rmv_beta_ids = find(isnan(mu_col));

        % Report problem
        disp(['  # NaN betas=',num2str(numel(rmv_beta_ids))]);

        % Find good beta maps
        keep_beta_ids = find(~isnan(mu_col));
        subj_img = subj_img(keep_beta_ids,:);

        disp(['  # kept betas=',num2str(numel(keep_beta_ids))]);

        %% Load labels
        filename = [proj.path.betas.fmri_id_ex_beta,'sub-',name,'_task-identify_ex_trials.tsv'];
        events = tdfread(filename);

        %% Prepare output structure
        prds = struct();
  
        %% ----------------------------------------
        %% VALENCE Machine Learning
        scores = events.valence(keep_beta_ids);
        labels = 0*scores;
        labels(find(scores>=proj.param.mvpa.likert)) = proj.param.mvpa.pos_cls;
        labels(find(scores<proj.param.mvpa.likert)) = proj.param.mvpa.neg_cls;
        
        %% Run machine learning        
        prds.v_cls_acc = [];
        prds.v_cls_hd = [];
        for j=1:proj.param.mvpa.Nresample
            
            %% Fit the data space
            [v_tst_hd,v_cls_stats] = classify_loocv2(proj,subj_img,labels);

            %% Store results
            prds.v_cls_acc = [prds.v_cls_acc;cell2mat(v_cls_stats.tst_acc)];
            prds.v_cls_hd = [prds.v_cls_hd;v_tst_hd'];
            prds.ids = keep_beta_ids;

        end

        %% debug
        all_v_cls_acc = [all_v_cls_acc;mean(mean(prds.v_cls_acc,1))];
        logger(['  v acc: ',num2str(mean(mean(prds.v_cls_acc,2)))],proj.path.logfile);

        %% ----------------------------------------
        %% AROUSAL Machine Learning

        %% Load/format labels
        scores = events.arousal(keep_beta_ids);
        labels = 0*scores;
        labels(find(scores>=proj.param.mvpa.likert)) = proj.param.mvpa.pos_cls;
        labels(find(scores<proj.param.mvpa.likert)) = proj.param.mvpa.neg_cls;
        
        %% Run machine learning        
        prds.a_cls_acc = [];
        prds.a_cls_hd = [];
        for j=1:proj.param.mvpa.Nresample
            
            %% Fit the data space
            [a_tst_hd,a_cls_stats] = classify_loocv2(proj,subj_img,labels);

            %% Store results
            prds.a_cls_acc = [prds.a_cls_acc;cell2mat(a_cls_stats.tst_acc)];
            prds.a_cls_hd = [prds.a_cls_hd;a_tst_hd'];
            prds.ids = keep_beta_ids;

        end

        %% debug
        all_a_cls_acc = [all_a_cls_acc;mean(mean(prds.a_cls_acc,1))];
        logger(['  a acc: ',num2str(mean(mean(prds.a_cls_acc,2)))],proj.path.logfile);

        %% ----------------------------------------
        %% SAVE Modeling Results
        save([proj.path.mvpa.id_ex_gm_cls,'sub-',name,'_prds.mat'],'prds');

    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end

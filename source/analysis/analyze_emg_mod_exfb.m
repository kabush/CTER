%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyze MOD_EXFB Valence Induction via EMG '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.emg_mod_exfb]);
    eval(['! rm -rf ',proj.path.analysis.emg_mod_exfb]);
    disp(['Creating ',proj.path.analysis.emg_mod_exfb]);
    eval(['! mkdir ',proj.path.analysis.emg_mod_exfb]);
end

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{2};    
Nscans = proj.param.mri.Nscans(2); 
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;

%% ----------------------------------------
%% Load Subjects
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load Subject Age/Sex 
filename = [proj.path.bids,'participants.tsv'];
participants = tdfread(filename);

%% ----------------------------------------
%% Extract Targets and Predictors
all_betas = [];
all_trial_types = [];
all_valences = [];
all_arousals = [];

scans = [];
subj_ids = [];
sexes = [];
ages = [];
v_acc = [];

for i = 1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % debug
    logger([subj_study,':',name],proj.path.logfile);

    for j = 1:Nscans
        
        %% Load betas for the task
        filename = [proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_betas.mat'];
        load(filename);
        betas = betas.(['mod',num2str(j)])';
        all_betas = [all_betas;betas'];

        %% Load trial types
        filename = [proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_trial_types.mat'];
        load(filename);
        trial_types = trial_types.(['mod',num2str(j)]);
        all_trial_types = [all_trial_types;trial_types];
        
        %% Load affect information
        filename = [proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_valences.mat'];
        load(filename);
        vals = valences.(['mod',num2str(j)]);
        all_valences = [all_valences;vals];

        filename = [proj.path.betas.emg_mod_exfb_beta,'sub-',name,'_mod_arousals.mat'];
        load(filename);
        aros = arousals.(['mod',num2str(j)]);
        all_arousals = [all_arousals;aros];

        %% Load subject's model performance estimate
        filename = [proj.path.mvpa.id_ex_gm_cls,'sub-',name,'_prds.mat'];
        load(filename);
        sbj_v_acc = mean(mean(prds.v_cls_acc));

        %% Gather all data: add in age, sex, and random effects
        search_id = find(participants.participant_id==str2double(name));
        age = participants.age(search_id);
        if(strcmp(participants.sex(search_id),'F'))
            sex = 0;
        else
            sex = 1;
        end
        ages = [ages;repmat(age,numel(betas),1)];
        sexes = [sexes;repmat(sex,numel(betas),1)];
        subj_ids = [subj_ids;repmat(i,numel(betas),1)];
        scans = [scans;repmat(j,numel(betas),1)];
        v_acc = [v_acc;repmat(sbj_v_acc,numel(betas),1)];

    end
    
end

%%  Covert all values to double precision (lme requirement)
beta = double(zscore(all_betas));
trial_type = double(zscore(all_trial_types));
scan = double(zscore(scans));
age = double(zscore(ages));
sex = double(zscore(sexes));
subj_id = double(zscore(subj_ids));
val = double(zscore(all_valences));
aro = double(zscore(all_arousals));
v_acc = double(zscore(v_acc));

%%  Construct data table
tbl = table(trial_type,beta,scan,age,sex,subj_id,val,aro,v_acc,...
            'VariableNames',{'trial_type','beta','scan','age','sex','subj','val','aro','v_acc'});

%%  Fit models (both GLM and GLMM)
mdl_fe = fitlme(tbl,['beta ~ 1 + trial_type + trial_type*val +' ...
                    'val + aro + age + sex + v_acc']);

mdl_re = fitlme(tbl,['beta ~ 1 + trial_type + trial_type*val +' ...
                    'val + aro + age + sex + v_acc + (1|subj)']);


%%  Compare GLM vs GLMM for model fit (control for complexity)
fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
re_flag = 0;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    re_flag = 1;
    logger('  *Random effects matter',proj.path.logfile);
end

%% Save Model
save([proj.path.analysis.emg_mod_exfb,'mdl.mat'],'mdl');

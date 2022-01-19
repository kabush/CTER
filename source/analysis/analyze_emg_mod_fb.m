%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyze MOD_FB Valence Induction via EMG '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.emg_mod_fb]);
    eval(['! rm -rf ',proj.path.analysis.emg_mod_fb]);
    disp(['Creating ',proj.path.analysis.emg_mod_fb]);
    eval(['! mkdir ',proj.path.analysis.emg_mod_fb]);
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
betas = [];
onsets = [];
subj_ids = [];
sexes = [];
ages = [];
scans = [];

for i = 1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % debug
    logger([subj_study,':',name],proj.path.logfile);

    for j = 1:Nscans
        
        %% Load bids events files for the task
        filename = [proj.path.betas.emg_mod_fb_beta,'sub-',name,'_mod_fb_betas.mat'];
        load(filename);
        fb_betas = fb_betas.(['mod',num2str(j)]);
        betas = [betas;zscore(fb_betas)];

        %% Load bids events files for the task
        filename = [proj.path.betas.emg_mod_fb_beta,'sub-',name,'_mod_fb_onsets.mat'];
        load(filename);
        fb_onsets = fb_onsets.(['mod',num2str(j)]);
        
        %% Process onsets to separate tasks
        diffs = diff(fb_onsets);
        start_idx = [1;find(diffs>2.5)+1];
        
        %% Extract FB tasks and shift start to 0
        for j = 1:numel(start_idx)
            if(j<numel(start_idx))
                seq = start_idx(j):(start_idx(j+1)-1);
            else
                seq = start_idx(j):numel(fb_onsets);
            end
            fb_onsets(seq) = fb_onsets(seq)-fb_onsets(start_idx(j));
        end
        
        
        % Add in age, sex, and random effects
        search_id = find(participants.participant_id==str2double(name));
        age = participants.age(search_id);
        if(strcmp(participants.sex(search_id),'F'))
            sex = 0;
        else
            sex = 1;
        end
        
        onsets = [onsets;fb_onsets];
        ages = [ages;repmat(age,numel(fb_onsets),1)];
        sexes = [sexes;repmat(sex,numel(fb_onsets),1)];
        subj_ids = [subj_ids;repmat(i,numel(fb_onsets),1)];
        scans = [scans;repmat(j,numel(fb_onsets),1)];

    end
    
end

logger(['Task time Predicting EMG (All Feedback trials)'],proj.path.logfile);

beta = double(betas);
onset = double(onsets);
scan = double(scans);
age = double(ages);
sex = double(sexes);
subj_id = double(subj_ids);

tbl = table(onset,beta,scan,age,sex,subj_id,...
            'VariableNames',{'onset','beta','scan','age','sex','subj'});

mdl_fe = fitlme(tbl,['beta ~ 1 + onset + scan + age + sex']);
mdl_re = fitlme(tbl,['beta ~ 1 + onset + scan + age + sex + (1|subj)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
re_flag = 0;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    re_flag = 1;
    logger('  *Random effects matter',proj.path.logfile);
end

%% Save Model
save([proj.path.analysis.emg_mod_fb,'mdl.mat'],'mdl');

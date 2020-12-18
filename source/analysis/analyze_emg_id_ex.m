%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyze EX Valence Induction via EMG (corrugator) '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.emg_id_ex]);
    eval(['! rm -rf ',proj.path.analysis.emg_id_ex]);
    disp(['Creating ',proj.path.analysis.emg_id_ex]);
    eval(['! mkdir ',proj.path.analysis.emg_id_ex]);
end

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{1};    
Nscans = proj.param.mri.Nscans(1); 
Nvol = proj.param.mri.Nvol(1);
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
valence = [];
arousal = [];
betas = [];
subj_ids = [];
sex = [];
age = [];

for i = 1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % debug
    logger([subj_study,':',name],proj.path.logfile);
    

    for j = 1:Nscans
        
        %% Load bids events files for the task (provides the labels)
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);
        
        %% Search through events for event type
        for k = 1:numel(events.onset)
            
            type = events.trial_type(k,:);
            
            %% Match event type
            if(strcmp(type,'ex_stim'))

                %Extract trial valence
                valence = [valence;str2double(events.valence(k,:))];
                arousal = [arousal;str2double(events.arousal(k,:))];

                %Extract subject id, age, and sex for control
                subj_ids = [subj_ids;i];
                search_id = find(participants.participant_id==str2double(name));
                age = [age;participants.age(search_id)];
                if(strcmp(participants.sex(search_id),'F'))
                    sex = [sex;0];
                else
                    sex = [sex;1];
                end

            end
            
        end
        
        %% Load bids events files for the task
        filename = [proj.path.betas.emg_id_ex_beta,'sub-',name,'_id_ex_betas.mat'];
        load(filename);

        %% Concatenate beta values
        betas = [betas;ex_betas.(['id',num2str(j)])];

    end
    
end

%% ----------------------------------------
%% Detect measure specific physiological signal

logger(['EMG Predicting Valence (All Identify trials)'],proj.path.logfile);

valence = double(valence); 
betas = double(betas); 
age = double(age);
sex = double(sex);
subj_ids = double(subj_ids);

tbl = table(valence,betas,age,sex,subj_ids,...
            'VariableNames',{'true_val','betas','age','sex','subjs'});

mdl_fe = fitlme(tbl,['true_val ~ 1 + betas + age + sex']);
mdl_re = fitlme(tbl,['true_val ~ 1 + betas + age + sex + (betas|subjs)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
re_flag = 0;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    re_flag = 1;
    logger('  *Random effects matter',proj.path.logfile);
end

%% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);
[~,~,RE] = randomEffects(mdl);

logger(['  EMG Valence beta=',num2str(FE.Estimate(2)),', p=',num2str(FE.pValue(2))]);

%% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

%% Log Model
logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

save([proj.path.analysis.emg_id_ex,'val_mdl.mat'],'mdl');

%% ----------------------------------------
%% Detect orthogonal physiological signal

logger(['EMG Predicting Arousal (All Identify trials)'],proj.path.logfile);

arousal = double(arousal); 
betas = double(betas); 
age = double(age);
sex = double(sex);
subj_ids = double(subj_ids);

tbl = table(arousal,betas,age,sex,subj_ids,...
            'VariableNames',{'true_aro','betas','age','sex','subjs'});

mdl_fe = fitlme(tbl,['true_aro ~ 1 + betas + age + sex']);
mdl_re = fitlme(tbl,['true_aro ~ 1 + betas + age + sex + (betas|subjs)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
re_flag = 0;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    re_flag = 1;
    logger('  *Random effects matter',proj.path.logfile);
end

%% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);
[~,~,RE] = randomEffects(mdl);

%% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

%% Log Model
logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

logger(['  EMG Arousal beta=',num2str(FE.Estimate(2)),', p=',num2str(FE.pValue(2))]);

save([proj.path.analysis.emg_id_ex,'aro_mdl.mat'],'mdl');


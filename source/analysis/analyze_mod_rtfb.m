%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyze Real-time experimental outcomes  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.mod_rtfb]);
    eval(['! rm -rf ',proj.path.analysis.mod_rtfb]);
    disp(['Creating ',proj.path.analysis.mod_rtfb]);
    eval(['! mkdir ',proj.path.analysis.mod_rtfb]);
end

%% ----------------------------------------
%% Load Subjects
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load Participant Information
filename = [proj.path.bids,'participants.tsv'];
participants = tdfread(filename);

%% ----------------------------------------
%% Iterate over study subjects
onset = [];
trial = [];
avg_fb = [];
thresh = [];
subject = [];
age = [];
sex = [];
scan = [];

in_val_b1 = [];
in_val_b0 = [];
in_aro_b1 = [];
in_aro_b0 = [];


for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load labels
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-modulate1_events.tsv'];
        events1 = tdfread(filename);
        events1.scan = repmat(1,numel(events1.onset),1);
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-modulate2_events.tsv'];
        events2 = tdfread(filename);
        events2.scan = repmat(2,numel(events2.onset),1);

        %% Merge strucutres
        events = events1;
        f = fieldnames(events);
        for j=1:length(f)
            events.(f{j}) = [events.(f{j});events2.(f{j})];
        end
        
        %%Find subject IN-task performance (measured as beta0 and beta1)
        load([proj.path.analysis.id_in_affect,'sub-',name,'_v_in.mat']);
        in_val_b1 = [in_val_b1;repmat(in_true.b1,numel(events.onset),1)];
        in_val_b0 = [in_val_b0;repmat(in_true.b0,numel(events.onset),1)];        

        load([proj.path.analysis.id_in_affect,'sub-',name,'_a_in.mat']);
        in_aro_b1 = [in_aro_b1;repmat(in_true.b1,numel(events.onset),1)];
        in_aro_b0 = [in_aro_b0;repmat(in_true.b0,numel(events.onset),1)];

        %% Convert trial types to numeric
        trial_tmp = [];
        avg_fb_tmp = [];
        thresh_tmp = [];
        for j =1:numel(events.onset)
  
            trial_type = events.trial_type(j,:);
            
            %% Match event type to index (makes indexing easier
            %% than text from BIDS events file)
            if(strcmp(trial_type,'fb_init'))
                trial_tmp  = [trial_tmp;1];
            else
                if(strcmp(trial_type,'fb     '))
                    trial_tmp  = [trial_tmp;2];
                else
                    if(strcmp(trial_type,'fb_stim') | strcmp(trial_type,'em_stim'))
                        trial_tmp = [trial_tmp;3];
                    else
                        trial_tmp = [trial_tmp;0];
                    end
                end
            end

            avg_fb_tmp = [avg_fb_tmp;str2double(events.avg_feedback(j,:))];
            thresh_tmp = [thresh_tmp;str2double(events.threshold(j,:))];

        end
        
        % Collect Experiment Data
        onset = [onset;events.onset];
        trial = [trial;trial_tmp];
        avg_fb = [avg_fb;avg_fb_tmp];
        thresh = [thresh;thresh_tmp];
        scan = [scan;events.scan];

        % Collect Subject Data
        subject = [subject;repmat(i,numel(events.onset),1)];
        
        %Find subject in participants.tsv
        subj_id = find(participants.participant_id==str2double(name));
        age = [age;repmat(participants.age(subj_id),numel(events.onset),1)];
        sex_tmp = participants.sex(subj_id);
        if(strcmp(sex_tmp,'F'))
            sex = [sex;repmat(0,numel(events.onset),1)];
        else
            sex = [sex;repmat(1,numel(events.onset),1)];
        end
        
    catch
        logger(['  -MVPA Error: possible data'],proj.path.logfile);
    end

end

%% ----------------------------------------
%% Construct GLMM to FB stim trigger speed

fb_diff = [];
fb_trial = [];
fb_scan = [];
fb_subj = [];
fb_age = [];
fb_sex = [];
fb_in_val_b1 = [];
fb_in_val_b0 = [];
fb_in_aro_b1 = [];
fb_in_aro_b0 = [];

sbj_ids = unique(subject)
for i = 1:numel(sbj_ids)

    %% Select data for this subject only
    sbj_id = sbj_ids(i);
    idx = find(subject==sbj_id);

    trialx = trial(idx);
    onsetx = onset(idx);
    scanx = scan(idx);
    sbjx = subject(idx);
    agex = age(idx);
    sexx = sex(idx);

    in_val_b1x = in_val_b1(idx);
    in_val_b0x = in_val_b0(idx);
    in_aro_b1x = in_aro_b1(idx);
    in_aro_b0x = in_aro_b0(idx);

    %% Extract out FB start and end times
    id_stm = find(trialx==3);
    id_init = find(trialx==1);

    % handle 3 codign errors (see CTER2bids code for details)
    if(numel(id_stm)~=numel(id_init))
        disp('*****different lengths of fb start/end*****')
        miss_stm_id = -1;
        for j=1:(numel(id_stm)-1)
            if(id_stm(j)>id_init(j+1))
        
                miss_stm_id = j;
                break;
            end
        end
        disp(['error at id_stm(',num2str(miss_stm_id),')']);
        id_init = [id_init(1:miss_stm_id-1);id_init(miss_stm_id+1:end)];
    end

    %% Compute time for FB trigger
    fb_tx = [onsetx(id_stm)-onsetx(id_init)];
    fb_diff = [fb_diff;fb_tx]; %% Target

    scan_tmp = scanx(id_stm);
    trial1 = 1:numel(find(scan_tmp==1));
    trial2 = 1:numel(find(scan_tmp==2));
    fb_trial = [fb_trial;trial1';trial2'];  %% FB trial number (per scan)

    fb_scan = [fb_scan;scanx(id_stm)]; %% FB scan (1 or 2)
    fb_subj = [fb_subj;sbjx(id_stm)]; %% Subject ID (for random effects)
    fb_age =  [fb_age;agex(id_stm)]; 
    fb_sex =  [fb_sex;sexx(id_stm)];

    %% Gather measures of IN skill
    fb_in_val_b1 = [fb_in_val_b1;in_val_b1x(id_stm)];
    fb_in_val_b0 = [fb_in_val_b0;in_val_b0x(id_stm)];
    fb_in_aro_b1 = [fb_in_aro_b1;in_aro_b1x(id_stm)];
    fb_in_aro_b0 = [fb_in_aro_b0;in_aro_b0x(id_stm)];

end

%% Format double for GLMM
fb_diff = double(fb_diff);
fb_trial = double(fb_trial);
fb_scan = double(fb_scan);
fb_subj = double(fb_subj);
fb_age = double(fb_age);
fb_sex = double(fb_sex);
fb_in_val_b1 = double(fb_in_val_b1);
fb_in_val_b0 = double(fb_in_val_b0);
fb_in_aro_b1 = double(fb_in_aro_b1);
fb_in_aro_b0 = double(fb_in_aro_b0);

%% Build and solve GLMM (Check for sign. Random effects - ignore if not)
tbl = table(fb_diff,fb_trial,fb_scan,fb_in_val_b1,fb_in_val_b0,...
            fb_in_aro_b1,fb_in_aro_b0,fb_subj,fb_age,fb_sex,...
            'VariableNames',{'trg','trial','scan','age','sex','invb1','invb0',...
                    'inab1','inab0','subj'});
mdl_fe = fitlme(tbl,['trg ~ 1 + trial + scan + age + sex + invb1 '...
                    '+ invb0 + inab1 + inab0']);
mdl_re = fitlme(tbl,['trg ~ 1 + trial + scan + age + sex + invb1 '...
                    '+ invb0 + inab1 + inab0 + (trial|subj)']);
fe_vs_re = compare(mdl_fe,mdl_re);
mdl = mdl_fe;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    logger('  *Random effects matter',proj.path.logfile);
end

%% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

%% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

%% Log Model
logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

logger(['  Intercept=',num2str(FE.Estimate(1))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(1))],proj.path.logfile);
logger(['  Trial Beta=',num2str(FE.Estimate(2))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(2))],proj.path.logfile);
logger(['  Scan Beta=',num2str(FE.Estimate(3))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(3))],proj.path.logfile);
logger(['  Age Beta=',num2str(FE.Estimate(4))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(4))],proj.path.logfile);
logger(['  Sex Beta=',num2str(FE.Estimate(5))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(5))],proj.path.logfile);
logger(['  INvb1 Beta=',num2str(FE.Estimate(6))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(6))],proj.path.logfile);
logger(['  INvb0 Beta=',num2str(FE.Estimate(7))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(7))],proj.path.logfile);
logger(['  INab1 Beta=',num2str(FE.Estimate(8))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(8))],proj.path.logfile);
logger(['  INab0 Beta=',num2str(FE.Estimate(9))],proj.path.logfile);
logger(['           p=',num2str(FE.pValue(9))],proj.path.logfile);

%% Store model
save([proj.path.analysis.mod_rtfb,'mdl.mat'],'mdl');

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyzing Modulate FB Perf. via IN Perf.       '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.mod_fb_mdl_via_in_perf]);
    eval(['! rm -rf ',proj.path.analysis.mod_fb_mdl_via_in_perf]);
    disp(['Creating ',proj.path.analysis.mod_fb_mdl_via_in_perf]);
    eval(['! mkdir ',proj.path.analysis.mod_fb_mdl_via_in_perf]);
end

%% ----------------------------------------
%% Load Subjects
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load Participant Information
filename = [proj.path.bids,'participants.tsv']
participants = tdfread(filename);

%% ----------------------------------------
%% Iterate over study subjects
onset = [];
trial = [];
valence = [];
arousal = [];
pred_valence = [];
pred_arousal = [];
in_valence_b1 = [];
in_arousal_b1 = [];
in_valence_b0 = [];
in_arousal_b0 = [];
subject = [];
age = [];
sex = [];
v_acc = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load subject's model performance estimate
        filename = [proj.path.mvpa.id_ex_gm_cls,'sub-',name,'_prds.mat'];
        load(filename);
        sbj_v_acc = mean(mean(prds.v_cls_acc));

        %% Load labels
        filename = [proj.path.mvpa.mod_exfb_affect,'sub-',name,'_task-modulate_exfb_predictions.tsv'];
        predictions = tdfread(filename);

        %% Convert trial types to numeric
        trial_tmp = [];
        for j =1:numel(predictions.onset)

            trial_type = predictions.trial_type(j,:);

            %% Match event type to index (makes indexing easier
            %% than text from BIDS events file)
            if(strcmp(trial_type,'ex_stim'))
                trial_tmp  = [trial_tmp;0];
            end

            if(strcmp(trial_type,'fb_stim') | strcmp(trial_type,'em_stim'))
                trial_tmp  = [trial_tmp;1];
            end

            if(strcmp(trial_type,'fb_prev'))
                trial_tmp  = [trial_tmp;2];
            end

            if(strcmp(trial_type,'ex_prev'))
                trial_tmp  = [trial_tmp;3];
            end

        end

        %Find subject IN-task performance (measured as beta0 and beta1)
        load([proj.path.analysis.id_in_affect,'sub-',name,'_v_in.mat']);
        in_valence_b1 = [in_valence_b1;repmat(in_true.b1,numel(predictions.onset),1)];
        in_valence_b0 = [in_valence_b0;repmat(in_true.b0,numel(predictions.onset),1)];        

        load([proj.path.analysis.id_in_affect,'sub-',name,'_a_in.mat']);
        in_arousal_b1 = [in_arousal_b1;repmat(in_true.b1,numel(predictions.onset),1)];
        in_arousal_b0 = [in_arousal_b0;repmat(in_true.b0,numel(predictions.onset),1)];

        % Collect Experiment Data
        onset = [onset;predictions.onset];
        trial = [trial;trial_tmp];
        valence = [valence;predictions.valence];
        arousal = [arousal;predictions.arousal];
        pred_valence = [pred_valence;predictions.pred_valence];
        pred_arousal = [pred_arousal;predictions.pred_arousal];
        size(pred_valence)

        % Collect Subject Data
        subject = [subject;repmat(i,numel(predictions.onset),1)];
        v_acc = [v_acc;repmat(sbj_v_acc,numel(predictions.onset),1)];
        
        %Find subject in participants.tsv
        subj_id = find(participants.participant_id==str2double(name));
        age = [age;repmat(participants.age(subj_id),numel(predictions.onset),1)];
        sex_tmp = participants.sex(subj_id);
        if(strcmp(sex_tmp,'F'))
            sex = [sex;repmat(0,numel(predictions.onset),1)];
        else
            sex = [sex;repmat(1,numel(predictions.onset),1)];
        end

    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end

%% ----------------------------------------
%% Extract Data to Conduct Estimates
ids = find(trial==2); % find predictions of fb_prev
trg = double(pred_valence(ids)-nanmean(pred_valence(ids)));
age = double(age(ids));
sex = double(sex(ids));
subjs = double(subject(ids));
in_val_b1 = double(in_valence_b1(ids));
in_aro_b1 = double(in_arousal_b1(ids));
in_val_b0 = double(in_valence_b0(ids));
in_aro_b0 = double(in_arousal_b0(ids));
v_acc = double(v_acc(ids));

%% ----------------------------------------
%% Combined Feedback Trials

logger(['-ALL TRIALS (FB previous as a function of IN skill params-'],proj.path.logfile);

tbl = table(trg,in_val_b1,in_aro_b1,in_val_b0,in_aro_b0,age,sex,v_acc,...
            'VariableNames',{'trg','invb1','inab1','invb0', ...
                    'inab0','age','sex','v_acc'});

mdl = fitlme(tbl,['trg ~ 1 + invb1 + inab1 + invb0 + inab0 + v_acc' ...
                 '+ invb1*inab1 ' ...
                 '+ invb1*inab0 ' ...
                 '+ invb1*invb0 ' ...
                 '+ inab1*inab0' ...
                 '+ invb1*age ' ...
                 '+ invb0*age ' ...
                 '+ invb1*sex ' ...
                 '+ invb0*sex ' ...
                 '+ invb1*v_acc ' ...
                 '+ invb0*v_acc ' ...
                 '+ age*sex']);

save([proj.path.analysis.mod_fb_mdl_via_in_perf,'mdl.mat'],'mdl');

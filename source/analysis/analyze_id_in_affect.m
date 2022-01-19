%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyzing Identify (IN) Task Affect Processing  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.id_in_affect]);
    eval(['! rm -rf ',proj.path.analysis.id_in_affect]);
    disp(['Creating ',proj.path.analysis.id_in_affect]);
    eval(['! mkdir ',proj.path.analysis.id_in_affect]);
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
subject = [];
age = [];
sex = [];

for i = 1:numel(subjs)

    % extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    % debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        % Load labels
        filename = [proj.path.mvpa.id_in_affect,'sub-',name,'_task-identify_in_predictions.tsv'];
        predictions = tdfread(filename);

        % Convert trial types to numeric
        trial_tmp = [];
        for j =1:numel(predictions.onset)
            trial_type = predictions.trial_type(j,:);
            % Match event type
            if(strcmp(trial_type,'in_stim'))
                trial_tmp  = [trial_tmp;0];
            else
                trial_tmp  = [trial_tmp;1];
            end
        end
        
        % Gather Experiment Data
        onset = [onset;predictions.onset];
        trial = [trial;trial_tmp];
        valence = [valence;predictions.valence];
        arousal = [arousal;predictions.arousal];
        pred_valence = [pred_valence;predictions.pred_valence];
        pred_arousal = [pred_arousal;predictions.pred_arousal];

        % Gather Subject Data
        subject = [subject;repmat(i,numel(predictions.onset),1)];
        
        % Find subject in participants.tsv
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

cue_ids = find(trial==0);
recall_ids = find(trial==1);

%% ----------------------------------------
%% Set-up IN-Cue out-of-sample prediction

true_val = double(valence(cue_ids));
pred_val = double(pred_valence(cue_ids));

true_aro = double(arousal(cue_ids));
pred_aro = double(pred_arousal(cue_ids));

cue_sex = double(sex(cue_ids));
cue_age = double(age(cue_ids));
cue_subj = double(subject(cue_ids));

%% ----------------------------------------
%% Set-up IN-Recall out-of-sample prediction

% Organize cues for repeated measures GLMM
base_cue_val = pred_valence(cue_ids);
cue_val = repmat(base_cue_val,1,proj.param.in.Nfeel);
cue_val = double(reshape(cue_val',1,prod(size(cue_val)))');

base_cue_aro = pred_arousal(cue_ids);
cue_aro = repmat(base_cue_aro,1,proj.param.in.Nfeel);
cue_aro = double(reshape(cue_aro',1,prod(size(cue_aro)))');

% Organize task lag for repeated measures GLMM
traj = 1:proj.param.in.Nfeel;
traj = repmat(traj,numel(cue_ids),1);
traj = double(reshape(traj',1,prod(size(traj)))');

% Organize recall for repeated measures GLMM
recall_val = double(pred_valence(recall_ids));
recall_aro = double(pred_arousal(recall_ids));

% Organize subject ID, sex, age for repeated measures GLMM
recall_sex = double(sex(recall_ids));
recall_age = double(age(recall_ids));
recall_subj = double(subject(recall_ids));

%% ----------------------------------------
%% VALENCE IN-Cue prediction

logger(['-VAL IN-Cue Prediction------------'],proj.path.logfile);

tbl = table(true_val,pred_val,cue_subj,cue_sex,... 
            cue_age,'VariableNames',{'trg','pred','subj','sex','age'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred + sex + age']);
mdl_re = fitlme(tbl,['trg ~ 1 + pred + sex + age + (pred|subj)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    logger('  *Random effects matter',proj.path.logfile);
end

% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(['  Cue Beta=',num2str(FE.Estimate(2))],proj.path.logfile);
logger(['  Cue Beta p=',num2str(FE.pValue(2))],proj.path.logfile);
logger(['  Sex Beta=',num2str(FE.Estimate(3))],proj.path.logfile);
logger(['  Sex Beta p=',num2str(FE.pValue(3))],proj.path.logfile);
logger(['  Age Beta=',num2str(FE.Estimate(4))],proj.path.logfile);
logger(['  Age Beta p=',num2str(FE.pValue(4))],proj.path.logfile);

% Save group effect model
save([proj.path.analysis.id_in_affect,'group_v_in_cue_pred.mat'],'mdl');

%% ----------------------------------------
%% AROUSAL IN-Cue prediction

logger(['-ARO IN-Cue Prediction------------'],proj.path.logfile);

tbl = table(true_aro,pred_aro,cue_subj,cue_sex,... 
            cue_age,'VariableNames',{'trg','pred','subj','sex','age'});

mdl_fe = fitlme(tbl,['trg ~ 1 + pred + sex + age']);
mdl_re = fitlme(tbl,['trg ~ 1 + pred + sex + age + (pred|subj)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    logger('  *Random effects matter',proj.path.logfile);
end

% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(['  Cue Beta=',num2str(FE.Estimate(2))],proj.path.logfile);
logger(['  Cue Beta p=',num2str(FE.pValue(2))],proj.path.logfile);
logger(['  Sex Beta=',num2str(FE.Estimate(3))],proj.path.logfile);
logger(['  Sex Beta p=',num2str(FE.pValue(3))],proj.path.logfile);
logger(['  Age Beta=',num2str(FE.Estimate(4))],proj.path.logfile);
logger(['  Age Beta p=',num2str(FE.pValue(4))],proj.path.logfile);

% Save group effect model
save([proj.path.analysis.id_in_affect,'group_a_in_cue_pred.mat'],'mdl');

%% ----------------------------------------
%% VALENCE IN-recall effect size

logger(['-VAL IN-Recall Prediction--------'],proj.path.logfile);

tbl = table(recall_val,cue_val,traj,recall_subj,recall_sex,... 
            recall_age,'VariableNames',{'recall','cue','traj',...
                    'subj','sex','age'});

mdl_fe = fitlme(tbl,['recall ~ 1 + traj + sex + age + cue']);
mdl_re = fitlme(tbl,['recall ~ 1 + traj + sex + age + cue + (cue|subj)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    logger('  *Random effects matter',proj.path.logfile);
end

% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(['  Cue Beta=',num2str(FE.Estimate(2))],proj.path.logfile);
logger(['  Cue Beta p=',num2str(FE.pValue(2))],proj.path.logfile);
logger(['  Traj Beta=',num2str(FE.Estimate(3))],proj.path.logfile);
logger(['  Traj Beta p=',num2str(FE.pValue(3))],proj.path.logfile);

% Save group effect model
save([proj.path.analysis.id_in_affect,'group_v_in_recall.mat'],'mdl');

%% ----------------------------------------
%% AROUSAL IN-recall effect size

logger(['-ARO IN-Recall Prediction--------'],proj.path.logfile);

tbl = table(recall_aro,cue_aro,traj,recall_subj,recall_sex,... 
            recall_age,'VariableNames',{'recall','cue','traj',...
                    'subj','sex','age'});

mdl_fe = fitlme(tbl,['recall ~ 1 + traj + sex + age + cue']);
mdl_re = fitlme(tbl,['recall ~ 1 + traj + sex + age + cue + (cue|subj)']);

fe_vs_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_vs_re.pValue<.05)
    mdl = mdl_re;
    logger('  *Random effects matter',proj.path.logfile);
end

% Extract Fixed effects
[~,~,FE] = fixedEffects(mdl);

% Compute effect size
Rsqr = mdl.Rsquared.Adjusted;
Fsqr = Rsqr/(1-Rsqr);

logger(['  Rsqr_adj=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(['  Cue Beta=',num2str(FE.Estimate(2))],proj.path.logfile);
logger(['  Cue Beta p=',num2str(FE.pValue(2))],proj.path.logfile);
logger(['  Traj Beta=',num2str(FE.Estimate(3))],proj.path.logfile);
logger(['  Traj Beta p=',num2str(FE.pValue(3))],proj.path.logfile);

logger([' '],proj.path.logfile);

% Save group effect model
save([proj.path.analysis.id_in_affect,'group_a_in_recall.mat'],'mdl');

%% ----------------------------------------
%% Compare VALENCE IN-recall vs surrogate
in_surr_cmp_b_v = [];
in_surr_cmp_b_a = [];
for i = 1:numel(subjs)

    name = subjs{i}.name;
    subj_ids = find(subject==i);
    subj_trial = trial(subj_ids);
    
    try
        % ----------------------------------------
        % VALENCE IN Task fit
        subj_affect = pred_valence(subj_ids);
        in_true = calc_in(proj,subj_trial,subj_affect);
        save([proj.path.analysis.id_in_affect,'sub-',name,'_v_in.mat'],'in_true');

        % Compare surrogate
        load([proj.path.analysis.rest_in_affect,'sub-',name,'_v_in_surr.mat']);
        in_surr_cmp_b_v = [in_surr_cmp_b_v,in_true.b1-in_surr.b1];

    catch
        logger(['  -GLM Error: possible missing data (VAL)'],proj.path.logfile);
    end


    try
        % ----------------------------------------
        % AROUSAL IN Task fit
        subj_affect = pred_arousal(subj_ids);
        in_true = calc_in(proj,subj_trial,subj_affect);
        save([proj.path.analysis.id_in_affect,'sub-',name,'_a_in.mat'],'in_true');

        % Compare surrogate
        load([proj.path.analysis.rest_in_affect,'sub-',name,'_a_in_surr.mat']);
        in_surr_cmp_b_a = [in_surr_cmp_b_a,in_true.b1-in_surr.b1];
    catch
        logger(['  -GLM Error: possible missing data (ARO)'],proj.path.logfile);
    end

end

logger(['-VAL IN-RECALL Entrainment'],proj.path.logfile);
p = signrank(in_surr_cmp_b_v);
logger(['  Entrain p=',num2str(p)],proj.path.logfile);

logger(['-ARO IN-RECALL Entrainment'],proj.path.logfile);
p = signrank(in_surr_cmp_b_a);
logger(['  Entrain p=',num2str(p)],proj.path.logfile);
logger(['  ']);


%% ----------------------------------------
%% Plotting VALENCE IN-recall 

figure(1)
set(gcf,'color','w');
hold on;

%% First scatter all the raw data
Smark = 10;
scatter(cue_val,recall_val,Smark,'MarkerFaceColor',...
        proj.param.plot.white,'MarkerEdgeColor',...
        proj.param.plot.very_light_grey);

plt_tbl = table(cue_val,recall_val,'VariableNames',{'cue','recall'});
writetable(plt_tbl,'val_scatter.txt','Delimiter','\t');

%% Plot non-significant single subjs
for i = 1:numel(subjs)
    name = subjs{i}.name;
    try
        load([proj.path.analysis.id_in_affect,'sub-',name,'_v_in.mat']);
        if(in_true.p1>=0.05)
            plot(sort(in_true.cue),sort(in_true.cue)*in_true.b1+in_true.b0,'Color',...
                 proj.param.plot.light_grey,'LineWidth',1);
        end
    catch
        logger(['  -Plot: missing mdl'],proj.path.logfile);
    end
end

%% Plot significant signle subjs
for i = 1:numel(subjs)
    name = subjs{i}.name;
    try
        load([proj.path.analysis.id_in_affect,'sub-',name,'_v_in.mat']);
        if(in_true.p1<0.05)
            plot(sort(in_true.cue),sort(in_true.cue)*in_true.b1+in_true.b0,'Color',...
                 proj.param.plot.dark_grey,'LineWidth',2);
        end
    catch
        logger(['  -Plot: missing mdl'],proj.path.logfile);
    end
end

%% Plot group effect
load([proj.path.analysis.id_in_affect,'group_v_in_recall.mat'],'mdl');
[~,~,FE] = fixedEffects(mdl);
b1 = FE.Estimate(2);
b0 = FE.Estimate(1);

plot(sort(cue_val),sort(cue_val)*b1+b0,'r-','LineWidth',3);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
set(gca,'Layer','top');


%% Export hi-resolution figure
export_fig cue_recall_perf_val.png -r300
eval(['! mv ',proj.path.code,'cue_recall_perf_val.png ',proj.path.fig,'cue_recall_perf_val.png']);

%% ----------------------------------------
%% Plotting VALENCE IN-entrain

figure(2)
set(gcf,'color','w');
hold on;

colors = [1.0,0.0,0.0;
          0.8,0.8,0.8];
x = [in_surr_cmp_b_v'];
g = repmat(1,numel(in_surr_cmp_b_v),1);
g_labels = {' '};
kabBoxPlot(x,g,g_labels,colors,0.8);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = 22; proj.param.plot.axisLabelFontSize;

set(gcf,'position',[0,0,300,600]);

%% Export hi-resolution figure
export_fig cue_recall_entrain_hist_val.png -r300
eval(['! mv ',proj.path.code,'cue_recall_entrain_hist_val.png ',...
      proj.path.fig,'cue_recall_entrain_hist_val.png']);

writetable(table(x),'val_entrain.txt');

%% ----------------------------------------
%% Plotting AROUSAL IN-recall 

figure(3)
set(gcf,'color','w');
hold on;

%% First scatter all the raw data
Smark = 10;
scatter(cue_aro,recall_aro,Smark,'MarkerFaceColor',...
        proj.param.plot.white,'MarkerEdgeColor',...
        proj.param.plot.very_light_grey);

plt_tbl = table(cue_aro,recall_aro,'VariableNames',{'cue','recall'});
writetable(plt_tbl,'aro_scatter.txt','Delimiter','\t');

%% Plot non-significant single subjs
for i = 1:numel(subjs)
    name = subjs{i}.name;
    try
        load([proj.path.analysis.id_in_affect,'sub-',name,'_a_in.mat']);
        if(in_true.p1>=0.05)
            plot(sort(in_true.cue),sort(in_true.cue)*in_true.b1+in_true.b0,'Color',...
                 proj.param.plot.light_grey,'LineWidth',1);
        end
    catch
        logger(['  -Plot: missing mdl'],proj.path.logfile);
    end
end

%% Plot significant single subjs
for i = 1:numel(subjs)
    name = subjs{i}.name;
    try
        load([proj.path.analysis.id_in_affect,'sub-',name,'_a_in.mat']);
        if(in_true.p1<0.05)
            plot(sort(in_true.cue),sort(in_true.cue)*in_true.b1+in_true.b0,'Color',...
                 proj.param.plot.dark_grey,'LineWidth',2);
        end
    catch
        logger(['  -Plot: missing mdl'],proj.path.logfile);
    end
end

%% Plot group effect
load([proj.path.analysis.id_in_affect,'group_a_in_recall.mat'],'mdl');
[~,~,FE] = fixedEffects(mdl);
b1 = FE.Estimate(2);
b0 = FE.Estimate(1);

plot(sort(cue_val),sort(cue_val)*b1+b0,'r-','LineWidth',3);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
set(gca,'Layer','top');

%% Export hi-resolution figure
export_fig cue_recall_perf_aro.png -r300
eval(['! mv ',proj.path.code,'cue_recall_perf_aro.png ',proj.path.fig,'cue_recall_perf_aro.png']);


%% ----------------------------------------
%% Plotting AROUSAL IN-entrain

figure(4)
set(gcf,'color','w');
hold on;

colors = [1.0,0.0,0.0;
          0.8,0.8,0.8];
x = [in_surr_cmp_b_a'];
g = repmat(1,numel(in_surr_cmp_b_a),1);
g_labels = {' '};
kabBoxPlot(x,g,g_labels,colors,0.8);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = 22; 

set(gcf,'position',[0,0,300,600]);

%% Export hi-resolution figure
export_fig cue_recall_entrain_hist_aro.png -r300
eval(['! mv ',proj.path.code,'cue_recall_entrain_hist_aro.png ',...
      proj.path.fig,'cue_recall_entrain_hist_aro.png']);


writetable(table(x),'aro_entrain.txt');


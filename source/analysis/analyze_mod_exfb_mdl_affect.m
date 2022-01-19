%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyzing Modulate Task Affect Processing  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.mod_exfb_mdl_affect]);
    eval(['! rm -rf ',proj.path.analysis.mod_exfb_mdl_affect]);
    disp(['Creating ',proj.path.analysis.mod_exfb_mdl_affect]);
    eval(['! mkdir ',proj.path.analysis.mod_exfb_mdl_affect]);
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

            if(strcmp(trial_type,'fb_stim')) 
                trial_tmp  = [trial_tmp;1];
            end

            if(strcmp(trial_type,'fb_prev')) 
                trial_tmp  = [trial_tmp;2];
            end

            if(strcmp(trial_type,'em_stim'))
                trial_tmp  = [trial_tmp;3];
            end

            if(strcmp(trial_type,'em_prev'))
                trial_tmp  = [trial_tmp;4];
            end

            if(strcmp(trial_type,'ex_prev'))
                trial_tmp  = [trial_tmp;5];
            end


        end
        
        % Collect Experiment Data
        onset = [onset;predictions.onset];
        trial = [trial;trial_tmp];
        valence = [valence;predictions.valence];
        arousal = [arousal;predictions.arousal];
        pred_valence = [pred_valence;predictions.pred_valence];
        pred_arousal = [pred_arousal;predictions.pred_arousal];

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
%% Effect of model accuracy on fb_prev performance
fb_raw_ids = find(trial==1);
fb_prev_raw_ids = find(trial==2);
fb_prev_ids = fb_raw_ids(find(~isnan(pred_valence(fb_prev_raw_ids))));
fb_prev_val = pred_valence(fb_prev_ids);
fb_v_acc = v_acc(fb_prev_ids);


%% ----------------------------------------
%% Extract Data to Conduct Estimates

%% Extract out-of-sample ex_stim predictions
ex_raw_ids = find(trial==0);
ex_ids = ex_raw_ids(find(~isnan(pred_valence(ex_raw_ids))));
ex_true_val = double(valence(ex_ids));
ex_pred_val = double(pred_valence(ex_ids));
ex_true_aro = double(arousal(ex_ids));
ex_pred_aro = double(pred_arousal(ex_ids));
ex_age = double(age(ex_ids));
ex_sex = double(sex(ex_ids));
ex_subjs = double(subject(ex_ids));
ex_v_acc = double(v_acc(ex_ids));
ex_prev_raw_ids = find(trial==5);
ex_prev_ids = ex_prev_raw_ids(find(~isnan(pred_valence(ex_prev_raw_ids))));
ex_prev_pred_val = pred_valence(ex_prev_ids);

%% Extract out-of-sample fb_stim predictions
fb_raw_ids = find(trial==1);
fb_ids = fb_raw_ids(find(~isnan(pred_valence(fb_raw_ids))));
fb_true_val = double(valence(fb_ids));
fb_pred_val = double(pred_valence(fb_ids));
fb_true_aro = double(arousal(fb_ids));
fb_pred_aro = double(pred_arousal(fb_ids));
fb_age = double(age(fb_ids));
fb_sex = double(sex(fb_ids));
fb_subjs = double(subject(fb_ids));
fb_v_acc = double(v_acc(fb_ids));

%% Scaled Combined modulation data
exfb_trial = double(zscore([trial(ex_ids);trial(fb_ids)]));
true_val = double(zscore([ex_true_val;fb_true_val]));
pred_val = double(zscore([ex_pred_val;fb_pred_val]));
true_aro = double(zscore([ex_true_aro;fb_true_aro]));
pred_aro = double(zscore([ex_pred_aro;fb_pred_aro]));
age = double(zscore([ex_age;fb_age]));
sex = double(zscore([ex_sex;fb_sex]));
subjs = double(zscore([ex_subjs;fb_subjs]));
v_acc = double(zscore([ex_v_acc;fb_v_acc]));
prev_val = double(zscore([ex_prev_pred_val;fb_prev_val]));

%% ----------------------------------------
%% PREDICTION VS TRUE (Valence)

logger(['PREDICTION VS TRUE (Valence: EX trials only)'],proj.path.logfile);

tbl = table(ex_true_val, ex_pred_val, ex_pred_aro, ex_age, ex_sex, ex_subjs,...
            'VariableNames',{'true_val','pred_val','pred_aro','age','sex','subjs'});

mdl_fe = fitlme(tbl,['true_val ~ 1 + pred_val + pred_aro ' ...
                    '+ age + sex']);
mdl_re = fitlme(tbl,['true_val ~ 1 + pred_val + pred_aro ' ...
                    '+ age + sex + (pred_val|subjs) + (pred_aro|subjs)']);

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

%% Save model
save([proj.path.analysis.mod_exfb_mdl_affect,'group_v_mod_ex_pred.mat'],'mdl');

%% ----------------------------------------
%% PREDICTION VS TRUE (Arousal)

logger(['PREDICTION VS TRUE (Arousal: EX trials only)'],proj.path.logfile);

tbl = table(ex_true_aro, ex_pred_aro, ex_pred_val, ex_age, ex_sex, ex_subjs,...
            'VariableNames',{'true_aro','pred_aro','pred_val','age','sex','subjs'});

mdl_fe = fitlme(tbl,['true_aro ~ 1 + pred_aro + pred_val ' ...
                    '+ age + sex']);
mdl_re = fitlme(tbl,['true_aro ~ 1 + pred_aro + pred_val ' ...
                    '+ age + sex + (pred_aro|subjs) + (pred_val|subjs)']);

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

%% Save model
save([proj.path.analysis.mod_exfb_mdl_affect,'group_a_mod_ex_pred.mat'],'mdl');

%% ----------------------------------------
%% FEEDBACK VS EX (All Modulation trials)

logger(['FEEDBACK VS EX (All Modulation trials)'],proj.path.logfile);

tbl = table(exfb_trial, true_val, pred_val, prev_val, true_aro, pred_aro, v_acc, age, sex, subjs,...
            'VariableNames',{'trial','true_val','pred_val','prev_val' ...
                    'true_aro','pred_aro','v_acc', 'age','sex','subjs'});

mdl_fe = fitlme(tbl,['pred_val ~ 1 + true_val + prev_val + true_aro + v_acc + trial + ' ...
                    'trial*prev_val + trial*true_val + age + sex']);
mdl_re = fitlme(tbl,['pred_val ~ 1 + true_val + prev_val + true_aro + v_acc + trial + ' ...
                    'trial*prev_val + trial*true_val + age + sex + ' ...
                    '(1|subjs)']);


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

%% Save model
save([proj.path.analysis.mod_exfb_mdl_affect,'group_v_mod_fb.mat'],'mdl');

%% ----------------------------------------
%% 

figure(1)
set(gcf,'color','w');
hold on;

plot([0,0],[0;12],'k--','LineWidth',1);

Smark = 60;
y_val = fliplr(1:9);

cnt = 1;
for i=1:numel(y_val)

    if(FE.pValue(i+1)<0.05)
        
        if(FE.Estimate(i+1)>0)
            scatter(FE.Estimate(i+1),y_val(i),Smark,'MarkerFaceColor',...
                    'b','MarkerEdgeColor','b');
            lo = FE.Lower(i+1);
            hi = FE.Upper(i+1);
            plot([lo;hi],[y_val(i);y_val(i)],'b-','LineWidth',3);
            
        else
            
            i
            
            scatter(FE.Estimate(i+1),y_val(i),Smark,'MarkerFaceColor',...
                    'r','MarkerEdgeColor','r');
            lo = FE.Lower(i+1);
            hi = FE.Upper(i+1);
            plot([lo;hi],[y_val(i);y_val(i)],'r-','LineWidth',3);
        end
        
    else
        
        scatter(FE.Estimate(i+1),y_val(i),Smark,'MarkerFaceColor',...
                proj.param.plot.light_grey,'MarkerEdgeColor',...
                proj.param.plot.light_grey);
        lo = FE.Lower(i+1);
        hi = FE.Upper(i+1);
        plot([lo;hi],[y_val(i);y_val(i)],'Color',proj.param.plot.light_grey,'LineWidth',3);
        
    end
    
end

xlabel('Effect Estimate');

xlim([-0.1,0.15]);
ylim([0,10]);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
set(gca,'Layer','top');

names = {
'self-induct.:dcd. val. prev.',
'self-induct.:norm. val. stim.',
'sex',
'age',
'dcd. valence acc.',
'norm. arousal stim.',
'*dcd. valence prev.',
'norm. valence stim.',
'self-induction'}

set(gca,'ytick',[1:9],'yticklabel',names)

%% Export hi-resolution figure
export_fig mod_exfb_mdl_affect -r300
eval(['! mv ',proj.path.code,'mod_exfb_mdl_affect.png ',...
      proj.path.fig,'mod_exfb_mdl_affect.png']);

%% ----------------------------------------
%% FEEDBACK INDUCTION VS PASSIVE VALENCE (Induction effects)

logger(['FEEDBACK VS PASSIVE (Val. Induction)'],proj.path.logfile);

%% Control for significant decodings
sig_ids = find(v_acc>=.607);

%% Combined modulation data
exfb_trial = double([trial(ex_ids);trial(fb_ids)]);
subjs = double(([ex_subjs;fb_subjs]));
v_acc = double(([ex_v_acc;fb_v_acc]));
prev_val = double([ex_prev_pred_val;fb_prev_val]);


%% Extract only data for significant decodings
sig_ids = find(v_acc>=.607);

exfb_trial = exfb_trial(sig_ids);
subjs = subjs(sig_ids);
v_acc = v_acc(sig_ids);
prev_val = prev_val(sig_ids);

fb_prev_val = prev_val(find(exfb_trial==1));
fb_prev_subjs = subjs(find(exfb_trial==1));
unique_subjs = unique(fb_prev_subjs);

p = signrank(fb_prev_val,.5,'tail','right');
disp(p)

%% BOOTSTRAP (Measuring random effects
all_strap_val = [];

%% Sampling parameters
fnull = proj.param.bootstrap.fnull;
Nboot = proj.param.bootstrap.Nboot;

for i = 1:Nboot

    strap_val = [];

    % Nsubjs = numel(unique_subjs);
    boot_sbjs = unique_subjs(randsample(1:numel(unique_subjs), ...
                                        numel(unique_subjs)),true);

    for j=1:numel(unique_subjs);

        % Extract subject values
        sbj = boot_sbjs(j); %unique_subjs(j);
        sbj_idx = find(fb_prev_subjs==sbj);

        % Gather valence
        sbj_prev_val = fb_prev_val(sbj_idx);

        % Sample with replacement at the within subject level
        % This preservers random effects of subjects
        smp_val = sbj_prev_val(randsample(1:numel(sbj_prev_val), ...
                                          numel(sbj_prev_val), ...
                                          true));
        % Compute subject mean
        strap_val = [strap_val,nanmedian(smp_val)];
    end

    % Gather random sample of group distribution
    all_strap_val = [all_strap_val,nanmedian(strap_val)];

    if(mod(i,1000)==0)
       disp(numel(find(all_strap_val<fnull))/i);
    end

end

%% Test hypothesis: p(mu<.5) one-sided
logger(['  VAL Bootstrap[n=',num2str(Nboot),'], p(mu<.5)=',...
        num2str(numel(find(all_strap_val<fnull))/Nboot)],proj.path.logfile);

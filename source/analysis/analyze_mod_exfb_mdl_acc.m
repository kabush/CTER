%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Analyzing Mdl Accuracy on FB Effect Size  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.mod_exfb_mdl_acc]);
    eval(['! rm -rf ',proj.path.analysis.mod_exfb_mdl_acc]);
    disp(['Creating ',proj.path.analysis.mod_exfb_mdl_acc]);
    eval(['! mkdir ',proj.path.analysis.mod_exfb_mdl_acc]);
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
a_acc = [];

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
        sbj_a_acc = mean(mean(prds.a_cls_acc));

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
        a_acc = [a_acc;repmat(sbj_a_acc,numel(predictions.onset),1)];
        
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

%% Model impact of model accuracy on FB induction performance
[b stats] = robustfit(fb_v_acc,fb_prev_val)
mdl = struct()
mdl.b = b;
mdl.stats = stats;

mdl.Rsqr = corr(fb_prev_val,fb_v_acc*mdl.b(2)+mdl.b(1))^2;

%% Save model
save([proj.path.analysis.mod_exfb_mdl_acc,'fb_prev_mdl.mat'],'mdl');

%% Save data
tbl = table(fb_v_acc,fb_prev_val,'VariableNames',{'val_acc','fb_trig_val'});
writetable(tbl,[proj.path.analysis.mod_exfb_mdl_acc,'mod_exfb_mdl_acc_fb_plot.csv'],'Delimiter','\t');

%% Plotting Valence Trigger Hyperplane as a function of Decoder Accuracy
figure(1)
set(gcf,'color','w');
hold on;

xmin = min(fb_v_acc);
xmax = max(fb_v_acc);
ymin = min(fb_prev_val);
ymax = max(fb_prev_val);

xseq = (xmin-0.05):0.01:(xmax+0.05);
yplot = xseq*mdl.b(2)+mdl.b(1);
Smark = 30;

scatter(fb_v_acc,fb_prev_val,Smark,...
        'MarkerFaceColor',proj.param.plot.very_light_grey,...
        'MarkerEdgeColor',proj.param.plot.light_grey);
plot(xseq,yplot,'r','Linewidth',3);

dax = 0.05
xlim([xmin-dax,xmax+dax]);
ylim([ymin-dax,ymax+dax]);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% Export hi-resolution figure
export_fig mod_exfb_mdl_acc_val_fb_prev -r300
eval(['! mv ',proj.path.code,'mod_exfb_mdl_acc_val_fb_prev.png ',...
      proj.path.fig,'mod_exfb_mdl_acc_val_fb_prev.png']);



%% ----------------------------------------
%% Effect of model accuracy on ex_prev performance
ex_raw_ids = find(trial==0);
ex_prev_raw_ids = find(trial==5);
ex_prev_ids = ex_raw_ids(find(~isnan(pred_valence(ex_prev_raw_ids))));
ex_prev_val = pred_valence(ex_prev_ids);
ex_v_acc = v_acc(ex_prev_ids);

%% Model impact of model accuracy on FB induction performance
[b stats] = robustfit(ex_v_acc,ex_prev_val)
mdl = struct()
mdl.b = b;
mdl.stats = stats;

mdl.Rsqr = corr(ex_prev_val,ex_v_acc*mdl.b(2)+mdl.b(1))^2;

%% Save model
save([proj.path.analysis.mod_exfb_mdl_acc,'ex_prev_mdl.mat'],'mdl');

%% Save data
tbl = table(ex_v_acc,ex_prev_val,'VariableNames',{'val_acc','ex_trig_val'});
writetable(tbl,[proj.path.analysis.mod_exfb_mdl_acc,'mod_exfb_mdl_acc_ex_plot.csv'],'Delimiter','\t');

% %% Plotting Valence Trigger Hyperplane as a function of Decoder Accuracy
figure(2)
set(gcf,'color','w');
hold on;

xmin = min(ex_v_acc);
xmax = max(ex_v_acc);
ymin = min(ex_prev_val);
ymax = max(ex_prev_val);

xseq = (xmin-0.05):0.01:(xmax+0.05);
yplot = xseq*mdl.b(2)+mdl.b(1);
Smark = 30;

scatter(ex_v_acc,ex_prev_val,Smark,...
        'MarkerFaceColor',proj.param.plot.very_light_grey,...
        'MarkerEdgeColor',proj.param.plot.light_grey);
plot(xseq,yplot,'r','Linewidth',3);

dax = 0.05
xlim([xmin-dax,xmax+dax]);
ylim([ymin-dax,ymax+dax]);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% Export hi-resolution figure
export_fig mod_exfb_mdl_acc_val_ex_prev -r300
eval(['! mv ',proj.path.code,'mod_exfb_mdl_acc_val_ex_prev.png ',...
      proj.path.fig,'mod_exfb_mdl_acc_val_ex_prev.png']);


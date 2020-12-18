function in = calc_in(proj,subj_trial,subj_aff)

%% initialize results
in = struct();

%% Separate cue from recall
cue_ids =  find(subj_trial==0);
recall_ids = find(subj_trial==1);

cue_aff = subj_aff(cue_ids);
recall_aff = subj_aff(recall_ids);

%% format 
Nfeel = proj.param.in.Nfeel;
cue_box = repmat(cue_aff,1,Nfeel);
recall_box = reshape(recall_aff',Nfeel,numel(recall_aff)/Nfeel)';
traj_box = repmat(1:Nfeel,numel(cue_aff),1);

%% reformat (vectorize for GLM)
cue_form = reshape(cue_box',1,prod(size(cue_box)))';
recall_form = reshape(recall_box',1,prod(size(recall_box)))';
traj_form = reshape(traj_box',1,prod(size(traj_box)))';

%% Fit GLM
tbl = table(recall_form,cue_form,traj_form,'VariableNames', ...
            {'recall','cue','traj'});

sbj_mdl_fe = fitlme(tbl,['recall ~ 1 + cue + traj']);

%% Extract Fixed effects
[~,~,FE] = fixedEffects(sbj_mdl_fe);
        
in.cue = cue_aff;
in.b1 = FE.Estimate(2); % slope
in.b0 = FE.Estimate(1); % intercept
in.p1 = FE.pValue(2); %slope
in.p0 = FE.pValue(1); %intercept

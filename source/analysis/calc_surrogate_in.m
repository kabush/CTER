function surr = calc_surrogate_in(proj,name,affect_name)
    
surr = struct();

% define range of data to sample
Nrecall = proj.param.in.Nfeel;
Nvol = proj.param.mri.Nvol(3); % 3rd task is rest
rest_start = proj.param.rest.Ntrs_trans+1;
rest_end = Nvol-proj.param.rest.Ntrs_tail-(Nrecall+1);
rest_ids = rest_start:rest_end;

data = load([proj.path.mvpa.rest_affect,'sub-',name,'_',affect_name,'_rest.1D']);
cue_ids = randsample(1:numel(rest_ids),proj.param.rest.Nresample)';
cue_data = data(cue_ids);

if(sum(isnan(data))==0)
    
    %% format rest into IN surrogate trials
    cue_box = repmat(cue_data,1,Nrecall);
    recall_box = [];
    for j = 1:Nrecall
        recall_box = [recall_box,data(cue_ids+j+1)];
    end
    traj_box = repmat(1:Nrecall,proj.param.rest.Nresample,1);
    
    %% reformat (vectorize)
    cue_form = reshape(cue_box',1,prod(size(cue_box)))';
    recall_form = reshape(recall_box',1,prod(size(recall_box)))';
    traj_form = reshape(traj_box',1,prod(size(traj_box)))';
    
    %% build individual subject structures
    surr.name = name;
    
    %% Fit GLM
    tbl = table(recall_form,cue_form,traj_form,'VariableNames', ...
                {'recall','cue','traj'});
    
    sbj_mdl_fe = fitlme(tbl,['recall ~ 1 + cue + traj']);
    
    %% Extract Fixed effects
    [~,~,FE] = fixedEffects(sbj_mdl_fe);
    
    surr.cue = cue_data
    surr.b1 = FE.Estimate(2); % slope
    surr.b0 = FE.Estimate(1); % intercept
    surr.p1 = FE.pValue(2); %slope
    surr.p0 = FE.pValue(1); %intercept
    
end
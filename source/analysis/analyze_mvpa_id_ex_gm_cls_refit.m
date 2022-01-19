%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['********************************************'],proj.path.logfile);
logger([' Analyze MVPA GM Classification (Refit)   '],proj.path.logfile);
logger(['********************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.gm_cls_refit]);
    eval(['! rm -rf ',proj.path.analysis.gm_cls_refit]);
    disp(['Creating ',proj.path.analysis.gm_cls_refit]);
    eval(['! mkdir ',proj.path.analysis.gm_cls_refit]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

all_raw_acc_v = [];
all_raw_acc_a = [];
all_good_ids = []; % Control for missing beta predictions

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    disp(['sub-',name]);

    try

        % Load labels (for sizing)
        filename = [proj.path.betas.fmri_id_ex_beta,'sub-',name, ...
                    '_task-identify_ex_trials.tsv'];
        events = tdfread(filename);
        Nlabels = numel(events.onset);
        
        % Load classification details
        load([proj.path.mvpa.id_ex_gm_cls,'sub-',name,'_prds.mat']);

        prds_v = zeros(size(prds.v_cls_acc,1),Nlabels);
        prds_a = zeros(size(prds.a_cls_acc,1),Nlabels);
        prds_good_ids = zeros(1,Nlabels);
        
        prds_v(:,prds.ids) = prds.v_cls_acc;
        prds_a(:,prds.ids) = prds.a_cls_acc;
        prds_good_ids(1,prds.ids) = 1;
 
        % Take column means (explicit 1 to handle nrows=1 case)
        all_raw_acc_v = [all_raw_acc_v;mean(prds_v,1)]; 
        all_raw_acc_a = [all_raw_acc_a;mean(prds_a,1)]; 
        all_good_ids = [all_good_ids;prds_good_ids];        

    catch
        disp(' no predictions');
    end

end

%% ----------------------------------------
%% VALENCE "Knows-what-it_knows"

% Gather stats pre-refit
mu_acc = [];
for i = 1:size(all_raw_acc_v,1)
    mu_acc = [mu_acc; mean(all_raw_acc_v(i,find(all_good_ids(i,:)==1)))];
end
[h p ci_pre stats] = ttest(mu_acc);

% Refit
[grp_corr_perf_v,corr_ids_v] = refit_acc(all_raw_acc_v,all_good_ids);

% Gather stats post-refit
[h p ci stats] = ttest(grp_corr_perf_v);

logger(['Grp VAL accuracy, pre-refit=', ...
        num2str(mean(mean(all_raw_acc_v,2))),', CI=[',...
        num2str(ci_pre(1)),', ',num2str(ci_pre(2)),'], post-refit=', ...
        num2str(mean(grp_corr_perf_v)), ', CI=[',...
        num2str(ci(1)),', ',num2str(ci(2)),']'],...
        proj.path.logfile);

ncorr = numel(corr_ids_v);
[a b] = binofit(ncorr/2,ncorr,0.05);

rss_sig_subjs = find(grp_corr_perf_v>b(2));
disp(rss_sig_subjs)

sscnt = numel(rss_sig_subjs); 

logger(['Sing. Subj VAL significant post-refit=',num2str(sscnt), ...
        '/',num2str(size(all_raw_acc_v,1))],proj.path.logfile);


%% ----------------------------------------
%% AROUSAL "Knows-what-it_knows"

% Gather stats pre-refit
mu_acc = [];
for i = 1:size(all_raw_acc_a,1)
    mu_acc = [mu_acc; mean(all_raw_acc_a(i,find(all_good_ids(i,:)==1)))];
end
[h p ci_pre stats] = ttest(mu_acc);

% Refit
[grp_corr_perf_a,corr_ids_a] = refit_acc(all_raw_acc_a,all_good_ids);

% Gather stats post-refit
[h p ci stats] = ttest(grp_corr_perf_a);

logger(['Grp ARO accuracy, pre-refit=', ...
        num2str(mean(mean(all_raw_acc_a,2))),', CI=[',...
        num2str(ci_pre(1)),', ',num2str(ci_pre(2)),'], post-refit=', ...
        num2str(mean(grp_corr_perf_a)), ', CI=[',...
        num2str(ci(1)),', ',num2str(ci(2)),']'],...
        proj.path.logfile);

ncorr = numel(corr_ids_a);
[a b] = binofit(ncorr/2,ncorr,0.05);
sscnt = numel(find(grp_corr_perf_a>b(2)));



logger(['Sing. Subj ARO significant post-refit=',num2str(sscnt), ...
        '/',num2str(size(all_raw_acc_a,1))],proj.path.logfile);


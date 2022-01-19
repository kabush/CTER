%% Load in path data
load('proj.mat');

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

logger(['************************************************'],proj.path.logfile);
logger(['Generate Plot of Experiment Design'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Extract subject info.  All subjects receive the same
%% stimuli so we can use one subject to build this.
name = subjs{1}.name;

%% Storage
valence = [];
arousal = [];
trial_type = [];

%% ----------------------------------------
%% Identify scans
task = proj.param.mri.tasks{1};  
Nscans = proj.param.mri.Nscans(1); 
for i = 1:Nscans

    %% Load bids events files for the task
    filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(i),'_events.tsv'];
    events = tdfread(filename);
    
    for j = 1:numel(events.onset)

        type = events.trial_type(j,:);
            
        %% Match event type
        if(strcmp(type,'ex_stim'))
            trial_type = [trial_type,1];
            valence = [valence;str2double(events.valence(j,:))];
            arousal = [arousal;str2double(events.arousal(j,:))];
        end
        
        if(strcmp(type,'in_stim'))
            trial_type = [trial_type,2];
            valence = [valence;str2double(events.valence(j,:))];
            arousal = [arousal;str2double(events.arousal(j,:))];
        end

    end

end

%% ----------------------------------------
%% Modify scans
task = proj.param.mri.tasks{2};  
Nscans = proj.param.mri.Nscans(2); 
for i = 1:Nscans

    %% Load bids events files for the task
    filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(i),'_events.tsv'];
    events = tdfread(filename);
    
    for j = 1:numel(events.onset)

        type = events.trial_type(j,:);
            
        %% Match event type
        if(strcmp(type,'ex_stim'))
            trial_type = [trial_type,3];
            valence = [valence;str2double(events.valence(j,:))];
            arousal = [arousal;str2double(events.arousal(j,:))];
        end
        
        if(strcmp(type,'fb_stim') | strcmp(type,'em_stim'))
            trial_type = [trial_type,4];
            valence = [valence;str2double(events.valence(j,:))];
            arousal = [arousal;str2double(events.arousal(j,:))];
        end

    end

end


%% ----------------------------------------
%% Plot

figure(1);
set(gcf,'color','w');
hold on;

Smark = 40;

% Plot Design Stimuli in Affect Space
id_ex_ids = find(trial_type==1);
scatter(arousal(id_ex_ids),valence(id_ex_ids),Smark,...
        'MarkerFaceColor',proj.param.plot.light_grey,...
        'MarkerEdgeColor',proj.param.plot.light_grey);

id_in_ids = find(trial_type==2);
scatter(arousal(id_in_ids),valence(id_in_ids),Smark,...
        'MarkerFaceColor',proj.param.plot.orange,...
        'MarkerEdgeColor',proj.param.plot.orange);


mod_ex_ids = find(trial_type==3);
scatter(arousal(mod_ex_ids),valence(mod_ex_ids),Smark,...
        'MarkerFaceColor',proj.param.plot.blue,...
        'MarkerEdgeColor',proj.param.plot.blue);


mod_fb_ids = find(trial_type==4);
scatter(arousal(mod_fb_ids),valence(mod_fb_ids),Smark,...
        'MarkerFaceColor',proj.param.plot.red,...
        'MarkerEdgeColor',proj.param.plot.red);

% Rescale limits
xlim([1,8]);
ylim([1,9]);

% Change axis fong
hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

% Export
export_fig stim_design.png -r300
eval(['! mv ',proj.path.code,'stim_design.png ',proj.path.fig,'stim_design.png']);


%% ----------------------------------------
%% Stats of identify EX/IN stimuli (grouped pos/neg val)

% Valence analysis
ex_mu_v = valence(id_ex_ids);
in_mu_v = valence(id_in_ids);

% Arousal analysis
ex_mu_a = arousal(id_ex_ids);
in_mu_a = arousal(id_in_ids);

logger(['Identify EX Stimuli stats'],proj.path.logfile);
logger(['   -Valence: avg v=',num2str(mean(ex_mu_v)),' (',...
        num2str(std(ex_mu_v)),')'],proj.path.logfile);
logger(['   -Arousal: avg a=',num2str(mean(ex_mu_a)),' (',...
        num2str(std(ex_mu_a)),')'],proj.path.logfile);

logger(['Identify IN Stimuli stats'],proj.path.logfile);
logger(['   -Valence: avg v=',num2str(mean(in_mu_v)),' (',...
        num2str(std(ex_mu_v)),')'],proj.path.logfile);
logger(['   -Arousal: avg a=',num2str(mean(in_mu_a)),' (',...
        num2str(std(in_mu_a)),')'],proj.path.logfile);

% Comparison EX vs IN
v_p = ranksum(ex_mu_v,in_mu_v);
a_p = ranksum(ex_mu_a,in_mu_a);

logger(['Identify EX vs IN Stimuli stats'],proj.path.logfile);
logger(['   -Valence Comparison: p=',num2str(v_p)],proj.path.logfile);
logger(['   -Arousal Comparison: p=',num2str(a_p)], ...
       proj.path.logfile);
logger(['  '],proj.path.logfile);


%% ----------------------------------------
%% Stats of modulate EX/FB stimuli (grouped pos/neg val)

% Gather Valence stats of Feedback stimuli
pos_mod_fb_ids = mod_fb_ids(find(valence(mod_fb_ids)>proj.param.mvpa.likert));
neg_mod_fb_ids = mod_fb_ids(find(valence(mod_fb_ids)<proj.param.mvpa.likert));

pos_fb_mu_v = mean(valence(pos_mod_fb_ids));
pos_fb_sig_v = std(valence(pos_mod_fb_ids));
neg_fb_mu_v = mean(valence(neg_mod_fb_ids));
neg_fb_sig_v = std(valence(neg_mod_fb_ids));

logger(['Modulate FB Stimuli Affect stats'],proj.path.logfile);
logger(['   -Pos Valence Cluster: avg v=',num2str(pos_fb_mu_v),' (',...
        num2str(pos_fb_sig_v),')'],proj.path.logfile);
logger(['   -Neg Valence Cluster: avg v=',num2str(neg_fb_mu_v),' (',...
        num2str(neg_fb_sig_v),')'],proj.path.logfile);

% Gather Valence stats of Exstrinsic stimuli
pos_mod_ex_ids = mod_ex_ids(find(valence(mod_ex_ids)>proj.param.mvpa.likert));
neg_mod_ex_ids = mod_ex_ids(find(valence(mod_ex_ids)<proj.param.mvpa.likert));

pos_ex_mu_v = mean(valence(pos_mod_ex_ids));
pos_ex_sig_v = std(valence(pos_mod_ex_ids));
neg_ex_mu_v = mean(valence(neg_mod_ex_ids));
neg_ex_sig_v = std(valence(neg_mod_ex_ids));

logger(['Modulate EX Stimuli Affect stats'],proj.path.logfile);
logger(['   -Pos Valence Cluster: avg v=',num2str(pos_ex_mu_v),' (',...
        num2str(pos_ex_sig_v),')'],proj.path.logfile);
logger(['   -Neg Valence Cluster: avg v=',num2str(neg_ex_mu_v),' (',...
        num2str(neg_ex_sig_v),')'],proj.path.logfile);

% Compare valence properties between FB and EX stimuli (pos/neg VAL
% separately)
pos_v_p = ranksum(valence(pos_mod_fb_ids),valence(pos_mod_ex_ids));
neg_v_p = ranksum(valence(neg_mod_fb_ids),valence(neg_mod_ex_ids));

logger(['Modulate EX vs FB stimuli (Valence comparison)'],proj.path.logfile);
logger(['   -Pos Valence Clusters: p=',num2str(pos_v_p)],proj.path.logfile);
logger(['   -Neg Valence Clusters: p=',num2str(neg_v_p)],proj.path.logfile);

% Compare arousal properties between FB and EX stimuli (pos/neg VAL
% separately)
pos_a_p = ranksum(arousal(pos_mod_fb_ids),arousal(pos_mod_ex_ids));
neg_a_p = ranksum(arousal(neg_mod_fb_ids),arousal(neg_mod_ex_ids));

logger(['Modulate EX vs FB stimuli (Arousal comparison)'],proj.path.logfile);
logger(['   -Pos Valence Clusters: p=',num2str(pos_a_p)],proj.path.logfile);
logger(['   -Neg Valence Clusters: p=',num2str(neg_a_p)],proj.path.logfile);


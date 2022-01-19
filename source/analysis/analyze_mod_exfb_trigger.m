%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Analyze Trigger of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.mod_exfb_trigger]);
    eval(['! rm -rf ',proj.path.analysis.mod_exfb_trigger]);
    disp(['Creating ',proj.path.analysis.mod_exfb_trigger]);
    eval(['! mkdir ',proj.path.analysis.mod_exfb_trigger]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Define the task and type details
task = proj.param.mri.tasks{2};    % modulate
Nscans = proj.param.mri.Nscans(2); % # scans
Nvol = proj.param.mri.Nvol(2);
TR = proj.param.mri.TR;
stim_t = proj.param.mri.stim_t;

trial_type = [];
feedback = [];
avg_feedback = [];
threshold = [];

%% ----------------------------------------
%% Fit beta series for each subject
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% ----------------------------------------
    %% Construct onset times for regression (combined)
    for j = 1:Nscans

        %% Load bids events files for the task
        filename = [proj.path.bids,'sub-',name,'/func/sub-',name,'_task-',task,num2str(j),'_events.tsv'];
        events = tdfread(filename);

        %% Search through events for event type
        for k = 1:numel(events.onset)

            type = events.trial_type(k,:);

            %% Match event type
            if(strcmp(type,'fb_stim') | strcmp(type,'em_stim'))

                feedback = [feedback;str2double(events.feedback(k,:))];
                avg_feedback = [avg_feedback;str2double(events.avg_feedback(k,:))];
                threshold = [threshold; ...
                             str2double(events.threshold(k,:))];

            end

        end

    end
        
end

cln_avg_feedback = avg_feedback(find(~isnan(avg_feedback)));
ftrigger = numel(find(cln_avg_feedback>=.8))/numel(cln_avg_feedback);
logger(['Frac trigger >=.8 = ',num2str(ftrigger)]);
logger(['Median trigger = ',num2str(median(cln_avg_feedback))]);
p = signrank(cln_avg_feedback);
logger(['p = ',num2str(p)]);

% Save model
save([proj.path.analysis.mod_exfb_trigger,'cln_avg_feedback.mat'],'cln_avg_feedback');

% Plot the trigger information
figure(1)
set(gcf,'color','w');
hold on;

hist(cln_avg_feedback);

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

% Export the hi-res figure
export_fig hist_fb_trigger_thresh.png -r300
eval(['! mv ',proj.path.code,'hist_fb_trigger_thresh.png ',proj.path.fig,'hist_fb_trigger_thresh.png']);

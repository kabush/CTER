%% ----------------------------------------
%% BASE PATH
%% ----------------------------------------
base_path = '/home/kabush/';
bids_name = 'CTER';
base_name = 'CTER_rvsn_1';

%% ----------------------------------------
%% Initialize project param structure
proj = struct();

%% ----------------------------------------
%% Link Library Tools
base_lib = [base_path,'lib/'];

proj.path.tools.kablab = [base_lib,'kablab/'];
addpath(genpath(proj.path.tools.kablab));

proj.path.tools.scralyze = [base_lib,'scralyze/'];
addpath(genpath(proj.path.tools.scralyze));

proj.path.tools.export_fig = [base_lib,'export_fig/'];
addpath(genpath(proj.path.tools.export_fig));

proj.path.tools.nifti = [base_lib,'nifti/'];
addpath(genpath(proj.path.tools.nifti));

proj.path.tools.noboxplot = [base_lib,'noboxplot'];
addpath(genpath(proj.path.tools.noboxplot));

%% ----------------------------------------
%% Link Atlases Path
proj.path.atlas = [base_path,'atlas'];

%% ----------------------------------------
%% Project Name
proj.path.name = base_name;
proj.bids.name = bids_name;

%% ----------------------------------------
%% Workspace (code,bids,derivatives,...)
proj.path.home = [base_path,'workspace/'];
proj.path.bids = [proj.path.home,'bids/',proj.bids.name,'/'];
proj.path.data = [proj.path.home,'data/',proj.path.name,'/'];
proj.path.derivs = [proj.path.data,'derivatives/'];
proj.path.code = [proj.path.home,'code/',proj.path.name,'/'];
proj.path.log =[proj.path.code,'log/']; 
proj.path.fig = [proj.path.code,'fig/'];
proj.path.tmp = [proj.path.code,'tmp/'];

%% Subject Lists
proj.path.subj_list = [proj.path.code,'subj_lists/'];

%% Logging (creates a unique time-stampted logfile)
formatOut = 'yyyy_mm_dd_HH:MM:SS';
t = datetime('now');
ds = datestr(t,formatOut);
proj.path.logfile = [proj.path.log,'logfile_',ds,'.txt'];

%% ----------------------------------------
%% Derivatives output directory (All top-level names)
proj.path.mri.name = 'mri/';
proj.path.physio.name = 'physio/';
proj.path.betas.name = 'betas/';
proj.path.trg.name = 'target/';
proj.path.mvpa.name = 'mvpa/';
proj.path.analysis.name = 'analysis/';
proj.path.haufe.name = 'haufe/';

%% ----------------------------------------
%% Specific Data Paths

%% Preprocessing data
proj.path.fmriprep = [proj.path.derivs,'fmriprep/'];
proj.path.kubios = [proj.path.derivs,'kubios/'];
proj.path.mri.mri_clean = [proj.path.derivs,proj.path.mri.name,'mri_clean/'];
proj.path.mri.gm_mask = [proj.path.derivs,proj.path.mri.name,'gm_mask/'];
proj.path.physio.scr_clean = [proj.path.derivs,proj.path.physio.name,'scr_clean/'];
proj.path.physio.emg_clean = [proj.path.derivs,proj.path.physio.name,'emg_clean/'];

%% Beta-Series Paths (feature extraction)
proj.path.betas.fmri_id_ex_beta = [proj.path.derivs,proj.path.betas.name,'fmri_id_ex_beta/'];
proj.path.betas.fmri_id_in_beta = [proj.path.derivs,proj.path.betas.name,'fmri_id_in_beta/'];
proj.path.betas.fmri_mod_exfb_beta = [proj.path.derivs,proj.path.betas.name,'fmri_mod_exfb_beta/'];
proj.path.betas.scr_id_ex_beta = [proj.path.derivs,proj.path.betas.name,'scr_id_ex_beta/'];
proj.path.betas.emg_id_ex_beta = [proj.path.derivs,proj.path.betas.name,'emg_id_ex_beta/'];
proj.path.betas.scr_mod_ex_beta = [proj.path.derivs,proj.path.betas.name,'scr_mod_ex_beta/'];
proj.path.betas.scr_mod_exfb_beta = [proj.path.derivs,proj.path.betas.name,'scr_mod_exfb_beta/'];
proj.path.betas.emg_mod_ex_beta = [proj.path.derivs,proj.path.betas.name,'emg_mod_ex_beta/'];
proj.path.betas.emg_mod_exfb_beta = [proj.path.derivs,proj.path.betas.name,'emg_mod_exfb_beta/'];
proj.path.betas.scr_mod_fb_beta = [proj.path.derivs,proj.path.betas.name,'scr_mod_fb_beta/'];
proj.path.betas.emg_mod_fb_beta = [proj.path.derivs,proj.path.betas.name,'emg_mod_fb_beta/'];

%% MVPA paths (decoding models)
proj.path.mvpa.id_ex_gm_cls = [proj.path.derivs,proj.path.mvpa.name,'id_ex_gm_cls/'];
proj.path.mvpa.id_ex_gm_mdl = [proj.path.derivs,proj.path.mvpa.name,'id_ex_gm_mdl/'];

%% Applying EX models to IN, MOD, and REST trials
proj.path.mvpa.mod_exfb_affect = [proj.path.derivs,proj.path.mvpa.name,'mod_exfb_affect/'];
proj.path.mvpa.id_in_affect = [proj.path.derivs,proj.path.mvpa.name,'id_in_affect/'];
proj.path.mvpa.rest_affect = [proj.path.derivs,proj.path.mvpa.name,'rest_affect/'];

%% Analysis
proj.path.analysis.gm_cls_refit = [proj.path.derivs,proj.path.analysis.name,'gm_cls_refit/'];
proj.path.analysis.rest_in_affect = [proj.path.derivs,proj.path.analysis.name,'rest_in_affect/'];
proj.path.analysis.id_in_affect = [proj.path.derivs,proj.path.analysis.name,'id_in_affect/'];
proj.path.analysis.mod_exfb_trigger = [proj.path.derivs,proj.path.analysis.name,'mod_exfb_trigger/'];

proj.path.analysis.mod_exfb_affect = [proj.path.derivs,proj.path.analysis.name,'mod_exfb_affect/'];

proj.path.analysis.mod_exfb_mdl_acc = [proj.path.derivs,proj.path.analysis.name,'mod_exfb_mdl_acc/'];
proj.path.analysis.mod_exfb_mdl_affect = [proj.path.derivs,proj.path.analysis.name,'mod_exfb_mdl_affect/'];

proj.path.analysis.mod_fb_via_in_perf = [proj.path.derivs,proj.path.analysis.name,'mod_fb_via_in_perf/'];
proj.path.analysis.mod_fb_mdl_via_in_perf = [proj.path.derivs,proj.path.analysis.name,'mod_fb_mdl_via_in_perf/'];

proj.path.analysis.mod_rtfb = [proj.path.derivs,proj.path.analysis.name,'mod_rtfb/'];
proj.path.analysis.scr_id_ex = [proj.path.derivs,proj.path.analysis.name,'scr_id_ex/'];
proj.path.analysis.scr_mod_ex = [proj.path.derivs,proj.path.analysis.name,'scr_mod_ex/'];
proj.path.analysis.scr_mod_fb = [proj.path.derivs,proj.path.analysis.name,'scr_mod_fb/'];
proj.path.analysis.scr_mod_exfb = [proj.path.derivs,proj.path.analysis.name,'scr_mod_exfb/'];

proj.path.analysis.emg_id_ex = [proj.path.derivs,proj.path.analysis.name,'emg_id_ex/'];
proj.path.analysis.emg_mod_ex = [proj.path.derivs,proj.path.analysis.name,'emg_mod_ex/'];
proj.path.analysis.emg_mod_fb = [proj.path.derivs,proj.path.analysis.name,'emg_mod_fb/'];
proj.path.analysis.emg_mod_exfb = [proj.path.derivs,proj.path.analysis.name,'emg_mod_exfb/'];

%% Analysis
proj.path.haufe.id_ex_gm_mdl = [proj.path.derivs,proj.path.haufe.name,'id_ex_gm_mdl/'];

%% ----------------------------------------
%% Project Parameter Definitions

%% Data source
proj.param.studies = {bids_name};

%% Clean build the directories as it runs
proj.flag.clean_build = 1;

%% MRI fmriprep post-processing params
proj.param.mri.tasks = {'identify','modulate','rest'};
proj.param.mri.Nscans = [2,2,1];
proj.param.mri.Nvol = [282,315,225];
proj.param.mri.TR = 2;
proj.param.mri.stim_t = 2;
proj.param.mri.FD_thresh = .5;
proj.param.mri.anat_space = 'MNI152NLin2009cAsym';
proj.param.mri.bold_space = [proj.param.mri.anat_space,'_res-native'];
proj.param.mri.desc = 'preproc_bold';
proj.param.mri.gm_prob = .7;
proj.param.mri.fwhm = 8.0;
proj.param.mri.temp_hz = .0078;

%% Physio Preprocessing parameters
proj.param.physio.chan_scr = 3;      %Signal order from BIOPAC rig
proj.param.physio.chan_emg_corr = 5; %Signal order from BIOPAC rig

proj.param.physio.hz_scr = 2000; %Hz
proj.param.physio.scr.filt_med_samp = 0.01; % (Bach 2015)
proj.param.physio.scr.filt_high = 0.0159; % (Staib 2015)
proj.param.physio.scr.filt_low = 5;
proj.param.physio.scr.filt_type = 2; %% bi-directional

proj.param.physio.hz_emg = 2000; %Hz
proj.param.physio.emg.filt_low = 500.0; %% BIOPAC App. Note #241 (2016)
proj.param.physio.emg.filt_high = 10.0; %% BIOPAC App. Note #241 (2016)
proj.param.physio.emg.filt_type = 2; %% bi-directional

%% Physio Beta parameters
proj.param.betas.hirez = 20;

%% Task specific parameters
proj.param.in.Nfeel = 4;
proj.param.rest.Npseudo = 100;
proj.param.rest.Nmotion = 29;    % Max # of censors allowed
                                 % for resting state beta-series
                                 % to run without error
proj.param.rest.Nresample = 30;  
proj.param.rest.Ntrs_trans = 5;  % To account for startle
proj.param.rest.Ntrs_tail = 10;  % To account for HRF tail
proj.param.mod.fb_prev_t = 2.0;  % Start time of beta extraction
                                 % for brain state at momemnt of
                                 % real-time trigger

%% MVPA parameters
proj.param.mvpa.likert = 5; %Score denoting middle (Lang, 2009)
proj.param.mvpa.pos_cls = 1;
proj.param.mvpa.neg_cls = -1;
proj.param.mvpa.kernel = 'linear';
proj.param.mvpa.Nresample = 30;

%% Bootstrap parameters
proj.param.bootstrap.fnull = .5;
proj.param.bootstrap.Nboot = 10000;

%% Haufe parameters
proj.param.haufe.Npermute = 1000;
proj.param.haufe.Nchunk = 10;

%% Plotting parameters
proj.param.plot.axisLabelFontSize = 18;
proj.param.plot.circleSize = 10;
proj.param.plot.white = [1,1,1];
proj.param.plot.very_light_grey = [.9,.9,.9];
proj.param.plot.light_grey = [.8,.8,.8];
proj.param.plot.dark_grey = [.6,.6,.6];
proj.param.plot.axis_nudge = 0.1;
proj.param.plot.blue = [0,0,1];
proj.param.plot.orange = [1,.6,0];
proj.param.plot.red = [1,0,0];

%% ----------------------------------------
%% Seed random number generator
rng(1,'twister');

%% ----------------------------------------
%% Write out initialized project structure
save('proj.mat','proj');

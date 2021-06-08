%%========================================
%%========================================
%%
%% Keith Bush, PhD (2020)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

tic

%% ----------------------------------------
%% Clean up matlab environment
matlab_reset;

%% ----------------------------------------
%% Link all source code
addpath(genpath('./source/'));

%% ----------------------------------------
%% Initialize the projects directories and parameters.
init_project;

%% ----------------------------------------
%% Set-up top-level data directory structure
clean_project;

%% ----------------------------------------
%% Post-processing fmriprep outputs
fmriprep_process_fmri; % Final process after fmriprep
fmriprep_process_gm;   % Extract grey matter masks

%% ----------------------------------------
%% Analyze Identify EX with MVPA
calc_id_ex_beta; % Extract features

mvpa_id_ex_gm_cls; % intra-subj whole-brain GM MVPA classifications
                   % performs performance estimation using LOOCV
                   % also constructs and saves single model for
                   % application to Mod formats

analyze_mvpa_id_ex_gm_cls_refit; % Estimate true class performance
mvpa_id_ex_gm_mdl;  % Full MVPA model for IN,REST,Mod-EX,Mod-FB

%% ----------------------------------------
%% Analyze REST IN entrainment with MVPA (surrogate)
mvpa_rest_affect;       % Predict affect of rull rest
analyze_rest_in_affect; % Gen. REST-IN surrogates

%% ----------------------------------------
%% Analyze Identify IN with prefit MVPA
calc_id_in_beta;      % Extract features
mvpa_id_in_affect;    % Predict affect of Id-IN trials
analyze_id_in_affect; % Test IN task performance
                    %   1) Test out-sample IN-cue predict
                    %   2) Estimate base Id-IN performance
                    %   3) Control for Rest-IN surrogates

%% ----------------------------------------
%% Check Trigger Threshold
analyze_mod_exfb_trigger;

%% ----------------------------------------
%% Estimate Valence (Neutral versus Val+ State) ***AIM 1***
calc_mod_exfb_beta;         % Extract features
mvpa_mod_exfb_affect;       % Predict affect of Mod-EX & Mod-FB trials
analyze_mod_exfb_affect;    % Estimate Mod-EX vs Mod-FB

%% ------------------------------------------
%% Test for IN skill effects ***Aim 2***
analyze_mod_fb_via_in_perf; % Mod-EX vs Mod-FB effects
analyze_mod_rtfb;           % Stim trigger speed

%% ---------------------------------------
%% Analyze Task Stimuli/Design
analyze_design;

%% ----------------------------------------
%% Preprocess Physiological Signals (EX trials)
preprocess_scr;
preprocess_emg;
 
%% ----------------------------------------
%% Test emotion induction (ID EX)
calc_scr_id_ex_beta;
calc_emg_id_ex_beta;

%% ----------------------------------------
%% Test emotion induction (MOD EX)
calc_scr_mod_ex_beta;
calc_emg_mod_ex_beta;

%% ----------------------------------------
%% Analyze emotion induction (ID EX)
analyze_scr_id_ex;
analyze_emg_id_ex;

%% ----------------------------------------
%% Analyze emotion induction (MOD EX)
analyze_scr_mod_ex;
analyze_emg_mod_ex;

%% ----------------------------------------
%% Encoding evaluation
haufe_id_ex_gm_mdl;

toc

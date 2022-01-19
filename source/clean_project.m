%% Load in path data
load('proj.mat');

if(proj.flag.clean_build)

    %% Create project directory
    disp(['Creating data sub-directories']);

    %% Create all top-level directories
    eval(['! rm -rf ',proj.path.derivs,proj.path.mri.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.mri.name]);

    eval(['! rm -rf ',proj.path.derivs,proj.path.physio.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.physio.name]);

    eval(['! rm -rf ',proj.path.derivs,proj.path.betas.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.betas.name]);

    eval(['! rm -rf ',proj.path.derivs,proj.path.mvpa.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.mvpa.name]);

    eval(['! rm -rf ',proj.path.derivs,proj.path.analysis.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.analysis.name]);

    eval(['! rm -rf ',proj.path.derivs,proj.path.haufe.name]);
    eval(['! mkdir ',proj.path.derivs,proj.path.haufe.name]);

end
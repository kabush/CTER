%% Load in path data
load('proj.mat');

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mri.mri_clean]);
    eval(['! rm -rf ',proj.path.mri.mri_clean]);
    disp(['Creating ',proj.path.mri.mri_clean]);
    eval(['! mkdir ',proj.path.mri.mri_clean]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);

logger(['************************************'],proj.path.logfile);
logger(['fmriprep postprococessing of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************'],proj.path.logfile);

%% Pull top-level scan params
tasks = proj.param.mri.tasks;
Nscans = proj.param.mri.Nscans;

%% Postprocess fmriprep
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Make a directory to store the task and scan data
    eval(['! mkdir ',proj.path.mri.mri_clean,'sub-',name]);

    for j = 1:numel(proj.param.mri.tasks)

        for k = 1:Nscans(j);

            disp(['sub-',name,'_task-',tasks{j},num2str(k)]);
            
            if(strcmp(tasks{j},'rest'))
                infilename = [proj.path.fmriprep,'sub-',name,'/func/sub-',name,...
                            '_task-',tasks{j},'_desc-confounds_regressors.tsv'];
                outfilestub =  ['sub-',name,'_task-',tasks{j}]
            else
                infilename = [proj.path.fmriprep,'sub-',name,'/func/sub-',name,...
                            '_task-',tasks{j},num2str(k),'_desc-confounds_regressors.tsv'];
                outfilestub =  ['sub-',name,'_task-',tasks{j},num2str(k)]
            end

            %%Extract motion information for post-processign
            [censor,fd,motion] = make_motion_regressors(proj,infilename);

            %%write these to subject location in mri_clean
            filename = [outfilestub,'_censor.1D']
            path = [proj.path.mri.mri_clean,'sub-',name,'/',filename]
            dlmwrite(path,censor,' ');
            
            filename = [outfilestub,'_FD.1D']
            path = [proj.path.mri.mri_clean,'sub-',name,'/',filename]
            dlmwrite(path,fd,'');
            
            filename = [outfilestub,'_motion.1D']
            path = [proj.path.mri.mri_clean,'sub-',name,'/',filename]
            dlmwrite(path,motion,' ');

            %% write processing information to tmp for afni code
            dlmwrite([proj.path.tmp,'in_path.txt'],proj.path.fmriprep,'');
            dlmwrite([proj.path.tmp,'out_path.txt'],proj.path.mri.mri_clean,'');
            dlmwrite([proj.path.tmp,'subject.txt'],name,'');
            dlmwrite([proj.path.tmp,'task.txt'],tasks{j},'');
            dlmwrite([proj.path.tmp,'scan.txt'],num2str(k),'');
            dlmwrite([proj.path.tmp,'space.txt'],proj.param.mri.bold_space,'');
            dlmwrite([proj.path.tmp,'desc.txt'],proj.param.mri.desc,'');            
            dlmwrite([proj.path.tmp,'fwhm.txt'],num2str(proj.param.mri.fwhm),'');            
            dlmwrite([proj.path.tmp,'temp_hz.txt'],num2str(proj.param.mri.temp_hz),'');            

            %% Do the preprocessing
            eval(['! ',proj.path.code,'source/process/fmriprep_process_afni ',...
                  proj.path.tmp]);
            
            %% Clean-up
            eval(['! rm ',proj.path.tmp,'*']);

        end

    end

end

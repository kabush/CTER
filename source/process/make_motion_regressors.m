function [censor,fd,motion] = make_motion_regressors(proj,filename)

%% Load the full regressor output from fmriprep
confound = tdfread(filename);

%% Gather BIRC standard motion regressors (Power 24 regressor model)
motion = [];
motion = [motion,confound.trans_x];
motion = [motion,confound.trans_x_power2];
motion = [motion,parse_confound(confound.trans_x_derivative1)];
motion = [motion,parse_confound(confound.trans_x_derivative1_power2)];
motion = [motion,confound.trans_y];
motion = [motion,confound.trans_y_power2];
motion = [motion,parse_confound(confound.trans_y_derivative1)];
motion = [motion,parse_confound(confound.trans_y_derivative1_power2)];
motion = [motion,confound.trans_z];
motion = [motion,confound.trans_z_power2];
motion = [motion,parse_confound(confound.trans_z_derivative1)];
motion = [motion,parse_confound(confound.trans_z_derivative1_power2)];
motion = [motion,confound.rot_x];
motion = [motion,confound.rot_x_power2];
motion = [motion,parse_confound(confound.rot_x_derivative1)];
motion = [motion,parse_confound(confound.rot_x_derivative1_power2)];
motion = [motion,confound.rot_y];
motion = [motion,confound.rot_y_power2];
motion = [motion,parse_confound(confound.rot_y_derivative1)];
motion = [motion,parse_confound(confound.rot_y_derivative1_power2)];
motion = [motion,confound.rot_z];
motion = [motion,confound.rot_z_power2];
motion = [motion,parse_confound(confound.rot_z_derivative1)];
motion = [motion,parse_confound(confound.rot_z_derivative1_power2)];
motion = [motion,confound.csf];
motion = [motion,parse_confound(confound.csf_derivative1)];
motion = [motion,confound.white_matter];
motion = [motion,parse_confound(confound.white_matter_derivative1)];
motion = [motion,confound.global_signal];
motion = [motion,parse_confound(confound.global_signal_derivative1)];

%% Get framewise displacement
fd = parse_confound(confound.framewise_displacement);

%% Censor individual volume
censor = ones(size(fd));
bad = find(fd>proj.param.mri.FD_thresh);
censor(bad)=0;

%% Censor the succeeding volume
if(numel(bad)>0)
    if(bad(end)==numel(censor)) %make sure not to increase length
                                %of censor file ***CHANGE from earlier***
        censor(bad(1:end-1)+1)=0; 
    else
        censor(bad+1)=0;
    end
end

%% Censor singleton volumes
f=find(censor==1);
f_diff=f.*0;
f_diff(2:end)=diff(f);
bad_fs=[];
if isempty(f)==0
    for bad_loop=1:numel(f)
        if bad_loop~=numel(f)
            if f_diff(bad_loop)~=1&&f_diff(bad_loop+1)~=1
                bad_fs=[bad_fs f(bad_loop)];
            end
        elseif bad_loop==numel(f)
            if f_diff(bad_loop)~=1
                bad_fs=[bad_fs f(bad_loop)];
            end
        end
    end
end

censor(bad_fs)=0;

%% This is old way.  Better way is to use FD more wisely
%% during analysis.
%% %censor entire run if less than 50% of data is usable
%% run_length=numel(censor);
%% if numel(find(censor==1))./numel(censor)<.5
%%    censor(1:end)=0;
%% end
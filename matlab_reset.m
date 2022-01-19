%clear command window
clc;

% close all figures
try
    close('all','hidden');
catch % do nothing
end

% clear loaded functions;
clear('all');

%stop and delete timers
alltimer = timerfindall;
if ~isempty(alltimer)
    stop(alltimer);
    delete(alltime);
end

%close open files
fclose('all');

%reset warning status
warning('on','all');
lasterror('reset');
lastwarn('');

%remove <default>properties from root;
prop = get(0,'default');
propname = fieldnames(prop);
for iprop = 1:length(propname)
    set(0,propname{iprop},'remove');
end

%restore original <default> properties of root, load
%default path, run startup.m
matlabrc

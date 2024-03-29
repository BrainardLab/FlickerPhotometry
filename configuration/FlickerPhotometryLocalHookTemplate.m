function FlickerPhotometryLocalHook
% FlickerPhotometryLocalHook
%
% Configure things for working on the FlickerPhotometry project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUseProject('FlickerPhotometry') to set up for
% this project.  You then edit your local copy to match your configuration.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
fprintf('Running FlickerPhotometry local hook\n');

%% Specify project name and location
projectName = 'FlickerPhotometry';
projectBaseDir = tbLocateProject(projectName);
if (ispref(projectName))
    rmpref(projectName);
end

% If we ever needed some user/machine specific preferences, this is one way
% we could do that.
sysInfo = GetComputerInfo();
switch (sysInfo.localHostName)
    case 'eagleray'
        % DHB's desktop
        baseDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
 
    otherwise
        % Some unspecified machine, try user specific customization
        switch(sysInfo.userShortName)
            % Could put user specific things in, but at the moment generic
            % is good enough.
            otherwise
                baseDir = fullfile('/Users/',sysInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
        end
end

%% Set preferences for project output
%
% This will need to be locally configured.
outputBaseDir = '/Users/geoffreyaguirre/Documents/MATLAB/flicker photometry data placeholder'; 

% Set the preferences
setpref(projectName,'outputBaseDir',outputBaseDir);


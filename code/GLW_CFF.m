function GLW_CFF(varargin)
% Heterochromatic flicker photometry experiment using GLWindow
%
% Syntax:
%   GLW_CFF
%
% Description:
%    Implements a heterochromatic flicker photometry experiment using
%    GLWindow. Subjects see a red and green flickering light and use
%    keypresses to adjust a specified cone's contrast (L or M cone) to
%    minimize the flicker. u = adjust up, d = adjust down, and q = quit
%    adjustment. The program saves subjects' data and displays their chosen
%    contrast and adjustments on the screen. It includes optional code for
%    verifying display timing.
%
% Inputs
%    After starting the program, the experimenter is promoted to enter
%    subject ID (character vector) and session number (integer)
%
% Outputs:
%    All variables are automatically saved in a file named
%    'SubjectID_SessionNumber.mat' in a subject-specific folder. Figures
%    are saved separately in the same subject folder.
%
% Optional key/value pairs:
%    'calFile'             - character vector with name of calibration 
%                            file. Default is 'MetropsisCalibration'.
%
%    'adjustCone'          - character vector indicating which cone's
%                            contrast the user can adjust: 'L' or 'M'.
%                            Default is 'M'.
%
%    'steadyConeContrast'  - fixed contrast of non-adjustable cone. Must be
%                            a double between 0 and maxLContrast if
%                            adjusting M or between 0 and maxMContrast if
%                            adjusting L. Default is 0.12.
%
%    'viewDistance'        - double indicating viewing distance in mm.
%                            Default is 400.
%
%    'flickerRate'         - integer indicating flicker frequency in Hz.
%                            Default is 30.
%
%    'timeCheck'           - logical indicating whether to collect and plot
%                            stimulus timing information. Default is false.
%
%    'timeCheckDuration'   - when timeCheck is enabled, integer duration 
%                            (s) for collecting timing data. Default is 10.

% History:
%    06/05/19  dce       Wrote it. Visual angle conversion code from ar
%                        ('ImageSizeFromAngle.m')
%    06/06/19  dce       Minor edits, input error checking
%    06/13/19  dce       Added code to verify timing/find missed frames
%    06/21/19  dce       Added calibration routine
%    06/27/19  dce       Modified user adjustment process to avoid skipped
%                        frames
%    07/10/19  dce       Rewrote stimuli to be in terms of cone contrast
%                        values rather than rgb values.
%    07/19/19  dce       Cleaned up code. Added new key-value pairs,
%                        including options to change frame rate, identity
%                        of adjusted cone, and steady cone contrast

% Examples:
%{
    GLW_CFF
    GLW_CFF('viewDistance', 1000)
    GLW_CFF('adjustCone', 'L', 'steadyConeContrast', 0.16)
    GLW_CFF('timeCheck', true, 'maxFrames', 3600)
    GLW_CFF('calFile', 'EyeTrackerLCDTest')
%}

% Parse key-value pairs
p = inputParser;
p.addParameter('calFile', 'MetropsisCalibration', @(x) (ischar(x)));
p.addParameter('adjustCone', 'M', @(x) (ischar(x) & isscalar(x)));
p.addParameter('steadyConeContrast', 0.12, @(x) (isnumeric(x) & isscalar(x) & x >= 0));
p.addParameter('flickerRate', 30, @(x) (isnumeric(x) & isscalar(x) & x >= 0));
p.addParameter('viewDistance', 400, @(x) (isnumeric(x) & isscalar(x) & x > 0));
p.addParameter('timeCheck', false, @(x) (islogical(x)));
p.addParameter('timeCheckDuration', 10, @(x) (isnumeric(x) & isscalar(x)& x >= 0));
p.parse(varargin{:});

% Key parameters. Maximum cone contrast values vary by monitor and are
% based on a [0.5 0.5 0.5] linear RGB background
bgLinearRGB = [0.25 0.25 0.25]; % Background linear RGB values
nContrastSteps = 21;            % Number of possibilities for adjustable contrast
maxLContrast = 0.154933;        % Maximum L cone contrast for Display++
maxMContrast = 0.165757;        % Maximum M cone contrast for Display++
angle = 2;                      % Visual angle of stimulus (degrees)

% Check if chosen steady cone contrast is within monitor gamut
steadyContrast = p.Results.steadyConeContrast;
if strcmp(p.Results.adjustCone,'M') && (steadyContrast > maxLContrast)
    error('Chosen steady cone contrast %g is greater than max L contrast %g',...
        steadyContrast, maxLContrast);
elseif strcmp(p.Results.adjustCone,'L') && (steadyContrast > maxMContrast)
    error('Chosen steady cone contrast %g is greater than max M contrast %g',...
        steadyContrast, maxMContrast);
end

% Get information on display
disp = mglDescribeDisplays;
last = disp(end);                % We will be using the last display
frameRate = last.refreshRate;
screenSize = last.screenSizeMM;  % Screen dimensions in mm
centerHeight = screenSize(2) / 2;

% Check if chosen flicker rate divides into frame rate. We double the
% flicker rate because each flicker cycle includes a color transition (min 2 frames) 
if mod(frameRate, 2 * p.Results.flickerRate) == 0
    nHoldingFrames = frameRate / (2 * p.Results.flickerRate);
else
    error('Monitor frame rate (%g Hz) is not divisible by chosen flicker rate (2 * %g Hz)',...
        frameRate, p.Results.flickerRate);
end

% Prompt user to enter subject ID and session number
subjectID = input('Enter subject ID: ');
sessionNum = num2str(input('Enter session number: '));

% Create directory named SubjectID for saving data, if it doesn't exist already
outputDir = fullfile(getpref('FlickerPhotometry','outputBaseDir'),subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Create data file with name Subject ID_SessionNumber. Throw error if a
% file already exists with that name
fileName = [subjectID,'_', sessionNum, '.mat'];
fileLoc = fullfile(outputDir,fileName);
if (isfile(fileLoc))
    error('Specified output file %s already exists', fileName);
end

% Load calibration information
[cal,~] = LoadCalFile(p.Results.calFile,[],getpref('BrainardLabToolbox','CalDataFolder'));
load('T_cones_ss2', 'T_cones_ss2', 'S_cones_ss2'); % Cone fundamentals
cal = SetSensorColorSpace(cal,T_cones_ss2, S_cones_ss2);
cal = SetGammaMethod(cal,0);

% Fill table with adjustable cone contrast values. We assume three classes
% of cones where L cones are the first and M cones are the second.
adjustmentTable = zeros(3,nContrastSteps);
if strcmp(p.Results.adjustCone,'M') % Fill table with M cone contrast steps
    adjustmentTable(2,:) = linspace(0,maxMContrast,nContrastSteps);
else                                % Fill table with L cone contrast steps
    adjustmentTable(1,:) = linspace(0,maxLContrast,nContrastSteps);
end

% Convert contrast values to linear RGB values
for i = 1:nContrastSteps
    adjustmentTable(:,i) = contrastTorgb(cal, adjustmentTable(:,i), 'RGB', true);
end
adjustmentTablePos = randi(nContrastSteps); % Random initial position in table

% Create array to store adjustment history (room for 500 adjustments)
dataArray = zeros(3,500);
dataArray(:,1) = adjustmentTable(:,adjustmentTablePos);
dataArrayPos = 2; % Initial position in adjustment history array

try
    % Instructions window
    intro = GLWindow('BackgroundColor', [0 0 0], 'SceneDimensions',...
        screenSize, 'windowID', length(disp));
    intro.addText('A flashing red and green circle will appear on the screen',...
        'Center', [0 0.3 * centerHeight], 'FontSize', 75, 'Color', [1 1 1],...
        'Name', 'line1');
    intro.addText('Try to adjust the green light to minimize the flicker',...
        'Center', [0 0.1 * centerHeight], 'FontSize', 75, 'Color', [1 1 1],...
        'Name', 'line2');
    intro.addText('Press u to increase green intensity and press d to decrease green intensity',...
        'Center', [0 -0.1 * centerHeight], 'FontSize', 75, 'Color', [1 1 1],...
        'Name', 'line3');
    intro.addText('When you are done, press q to quit adjustment',...
        'Center', [0 -0.3 * centerHeight], 'FontSize', 75, 'Color', [1 1 1],...
        'Name', 'line4');
    intro.addText('*Press space bar to begin experiment*',...
        'Center', [0 -0.5 * centerHeight], 'FontSize', 75, 'Color', [0.5 0.5 1],...
        'Name', 'line5');
    
    % Open instructions window and enable character listening
    intro.open;
    mglDisplayCursor(0);
    ListenChar(2);
    FlushEvents;
    
    % Draw until the user presses the space bar 
    while true
        intro.draw;
        if CharAvail
            switch GetChar
                case ' ', break; otherwise, continue;
            end
        end
    end
    intro.close;
    
    % Create stimulus window
    win = GLWindow('BackgroundColor', bgLinearRGB, ...
        'SceneDimensions', screenSize, 'windowID', length(disp));
    
    % Calculate color of circle and diameter in mm. Then add circle.
    % L cone color is set to 12% l contrast
    if strcmp(p.Results.adjustCone,'M')
        steadyConeCol = contrastTorgb(cal, [steadyContrast 0 0], 'RGB', true);
    else
        steadyConeCol = contrastTorgb(cal, [0 steadyContrast 0], 'RGB', true);
    end
    
    diameter = tan(deg2rad(angle/2)) * (2 * p.Results.viewDistance);
    win.addOval([0 0], [diameter diameter], steadyConeCol, 'Name', 'circle');
    
    % Initial parameters.
    % isSteadyCone tells us whether it is the static or adjustable cone
    % stimulus on this frame. How fast we toggle this determines the flicker rate.
    isSteadyCone = true;
    elapsedFrames = 1;     % Tracks total number of elapsed frames
    
    
    % Create array for saving timing data for maxFrames or default duration
    if p.Results.timeCheck
        maxCheckedFrames = p.Results.timeCheckDuration * frameRate;
        timeStamps = zeros(1, maxCheckedFrames);
    end
    
    win.open;
    mglDisplayCursor(0);
    
    % Loop to swich oval color and parse user input
    while true
        % Draw circle
        if isSteadyCone
            color = steadyConeCol;
        else
            color = adjustmentTable(:,adjustmentTablePos)';
        end
        win.setObjectColor('circle', color);
        win.draw;
        
        % Switch color if needed
        if mod(elapsedFrames, nHoldingFrames) == 0
            isSteadyCone = ~isSteadyCone;
        end
        
        % Save timestamp
        if p.Results.timeCheck && elapsedFrames <= maxCheckedFrames
            timeStamps(elapsedFrames) = mglGetSecs;
        end
        elapsedFrames = elapsedFrames + 1;
        
        % Check for user input
        if CharAvail
            switch GetChar
                case 'q'    % Quit adjustment
                    break;
                case 'u'    % Adjust green up
                    adjustmentTablePos = adjustmentTablePos + 1;
                    if adjustmentTablePos > nContrastSteps
                        adjustmentTablePos = nContrastSteps;
                    end
                case 'd'    % Adjust green down
                    adjustmentTablePos = adjustmentTablePos - 1;
                    if adjustmentTablePos < 1
                        adjustmentTablePos = 1;
                    end
            end
            
            % Store new m value in adjustment history table
            dataArray(:,dataArrayPos) = adjustmentTable(:,adjustmentTablePos);
            dataArrayPos = dataArrayPos + 1;
        end
    end
    
    % Clean up once user finishes
    ListenChar(0);
    mglDisplayCursor(1);
    win.close;
    
    if p.Results.timeCheck
        % Plot frame durations
        timeSteps = diff(timeStamps);
        figure(1);
        plot(timeSteps, 'r');      % Actual frame durations
        yline(1/frameRate, 'b');   % Target frame rate
        yline(2/frameRate, 'g');   % Double target frame rate
        yline(0,'g');              %0 time
        axis([0 maxCheckedFrames 0 2.5/frameRate]);
        title('Frame Rate');
        xlabel('Frame');
        ylabel('Duration (s)');
        legend('Measured Frame Rate', 'Target Frame Rate', 'Skipped Frame');
        savefig(fullfile(outputDir, [fileName, '_frameRates.fig']));
        
        % Plot deviations from target frame rate
        figure(2);
        deviation = timeSteps - (1/frameRate);
        plot(deviation, 'r');
        axis([0 maxCheckedFrames -2/frameRate 2/frameRate]);
        title('Deviations from Frame Rate');
        xlabel('Frame');
        ylabel('Difference Between Measured and Target Duration (s)');
        savefig(fullfile(outputDir, [fileName, '_frameDeviations.fig']));
    end
    
    % Reformat adjustment history array
    dataArray = dataArray(dataArray ~= 0);
    col = length(dataArray)/3;
    dataArray = reshape(dataArray, [3 col]);
    for i = 1:col
        dataArray(:,i) = rgbToContrast(cal, dataArray(:,i)', 'RGB', true);
    end
    if strcmp(p.Results.adjustCone,'M')
        dataArray = dataArray(2,:);
    else
        dataArray = dataArray(1,:);
    end
    % Display results and save data
    fprintf('chosen %s cone contrast is %g \n', p.Results.adjustCone, dataArray(end));
    fprintf('adjustment history: ');
    if length(dataArray) > 1
        fprintf('%g, ', dataArray(1:end-1));
    end 
    fprintf('%g', dataArray(end)); 
    fprintf('\n');
    save(fileLoc);
    
    % Handle errors
catch e
    ListenChar(0);
    mglDisplayCursor(1);
    rethrow(e);
end
end

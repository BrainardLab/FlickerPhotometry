function GLW_CFF(fName, varargin)
% Heterochromatic flicker photometry experiment using GLWindow
%
% Syntax:
%   GLW_CFF
%
% Description:
%    Implements a heterochromatic flicker photometry experiment using
%    GLWindow. Subjects see a red and green light flickering at 30Hz and
%    use keypresses to adjust the m cone contrast to minimize the flicker
%    (u = adjust up, d = adjust down, q = quit adjustment). The program
%    saves subjects' adjustment history and displays their chosen m cone
%    contrast and adjustments on the screen. It is designed for displays
%    with 60 or 120 Hz frame rates and includes code for verifying timing.
%    Contrast values are calculated from an [0.5 0.5 0.5] RGB background.
%
% Inputs 
%    None 
%
% Outputs:
%    None
%
% Optional key/value pairs (can only use if you input a filename):
%    'viewDistance'   - double indicating viewing distance in mm. Default
%                       is 400
%
%    'maxFrames'      - double indicating the maximum number of frames,
%                       mostly used for timing verification. Default is Inf
%
%    'calFile'        - string with name of calibration file. Default is
%                       'MetropsisCalibration'

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
%    07/16/19  dce       Changed data saving procedure

% Examples:
%{
    GLW_CFF
    GLW_CFF('Deena.mat')
    GLW_CFF('Deena.mat', 'viewDistance', 1000)
    GLW_CFF('Deena.mat', 'maxFrames', 3600)
    GLW_CFF('Deena.mat', 'calFile', 'MetropsisCalibration')
    GLW_CFF('DHB.mat', 'calFile', 'EyeTrackerLCDTest')
%}

% Parse key-value pairs 
p = inputParser;
p.addParameter('viewDistance', 400, @(x) (isnumeric(x) & isscalar(x)));
p.addParameter('maxFrames', Inf, @(x) (isnumeric(x) & isscalar(x)));
p.addParameter('calFile', 'MetropsisCalibration', @(x) (ischar(x)));
p.parse(varargin{:});

%prompt user for subject ID and session number
subjectID = input('Enter subject ID (string): ');
sessionNum = num2str(input('Enter session number: '));

%Create directory and file for saving data
outputDir = fullfile(getpref('FlickerPhotometry','outputBaseDir'),subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
fileName = [subjectID,'_', sessionNum];
fileLoc = fullfile(outputDir,fileName);
if (exist(fileLoc,'file'))
    error(sprintf('Specified output file %s already exists',fName));
end

% Key parameters
bgLinearRGB = [0.25 0.25 0.25];         % Background linear RGB values
nMContrastSteps = 21;                   % Number of possibilities for m contrast
lConeContrast = 0.12;                   % Contrast for L cone phase of the flicker

% Get information on display
disp = mglDescribeDisplays;
last = disp(end); %we will be using the last display
frameRate = last.refreshRate;
screenSize = last.screenSizeMM; % Screen dimensions in mm
centerHeight = screenSize(2) / 2;

% Load calibration information
[cal,cals] = LoadCalFile(p.Results.calFile,[],getpref('BrainardLabToolbox','CalDataFolder'));
load T_cones_ss2; % Cone fundamentals
cal = SetSensorColorSpace(cal,T_cones_ss2, S_cones_ss2);
cal = SetGammaMethod(cal,0);

% Fill table with m contrast values. Contrast values go from 0 to 16.5757%
% (max contrast on the Metropsis display for a [0.5 0.5 0.5] RGB background)
% We assume three classes of cones and that the m cones are the second.
mArray = zeros(3,nMContrastSteps);
mArray(2,:) = 0:0.00828785:0.165757;

% Convert contrast values to RGB values
for i = 1:nMContrastSteps
    mArray(:,i) = contrastTorgb(cal, mArray(:,i), 'RGB', true);
end
mPosition = randi(nMContrastSteps); % Random initial position in table of values

% Create array to store adjustment history
adjustmentArray = zeros(3,100);
adjustmentArray(:,1) = mArray(:,mPosition);
adjustmentArrayPosition = 2; % Initial position in adjustment history array

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
        intro.addText('*Press any key to begin experiment*',...
        'Center', [0 -0.5 * centerHeight], 'FontSize', 75, 'Color', [0.5 0.5 1],...
        'Name', 'line5');
    
    % Open instructions window and enable character listening
    intro.open;
    mglDisplayCursor(0);
    ListenChar(2);
    FlushEvents;
    
    % Draw until the user presses a key 
    while ~CharAvail
        intro.draw;
    end
    intro.close;
    
    % Create stimulus window
    win = GLWindow('BackgroundColor', bgLinearRGB, ...
        'SceneDimensions', screenSize, 'windowID', length(disp));
    
    % Calculate color of circle and diameter in mm. Then add circle.
    % L cone color is set to 12% l contrast
    lCone = contrastTorgb(cal, [lConeContrast 0 0], 'RGB', true);
    angle = 2; % Visual angle (degrees)
    diameter = tan(deg2rad(angle/2)) * (2 * p.Results.viewDistance);
    win.addOval([0 0], [diameter diameter], lCone, 'Name', 'circle');
    
    % Initial parameters
    %
    % isLCone tells us whether it is L cone or M cone stimulus on this
    % frame. How fast we toggle this determines the flicker rate.
    isLCone = true; 
    
    elapsedFrames = 1; % This tracks total number of elapsed frames
    maxFrames = p.Results.maxFrames;
    if isfinite(maxFrames)
        timeStamps = zeros(1,maxFrames);
    end
    
    win.open;
    mglDisplayCursor(0);
    
    % Loop to swich oval color and parse user input
    while elapsedFrames <= maxFrames
        % Draw circle
        if isLCone
            color = lCone;
        else
            color = mArray(:,mPosition)';
        end
        win.setObjectColor('circle', color);
        win.draw;
        
        % Save timestamp
        if isfinite(maxFrames)
            timeStamps(elapsedFrames) = mglGetSecs;
        end
        elapsedFrames = elapsedFrames + 1;
        
        % Switch color if needed
        if (frameRate == 120 && mod(elapsedFrames, 2) == 1) || (frameRate == 60 && mod(elapsedFrames, 2) == 1)
            isLCone = ~isLCone;
        end
        
        % Check for user input
        if CharAvail
            switch GetChar
                case 'q' % Quit adjustment
                    break;
                case 'u' % Adjust green up
                    mPosition = mPosition + 1;
                    if mPosition > nMContrastSteps
                        mPosition = nMContrastSteps;
                    end
                case 'd' % Adjust green down
                    mPosition = mPosition - 1;
                    if mPosition < 1
                        mPosition = 1;
                    end
            end
            
            % Store new m value in adjustment history table
            adjustmentArray(:,adjustmentArrayPosition) = mArray(:,mPosition);
            adjustmentArrayPosition = adjustmentArrayPosition + 1;
        end
    end
    
    % Clean up once user finishes
    ListenChar(0);
    mglDisplayCursor(1);
    win.close;
    
    if isfinite(maxFrames)
        % Plot frame durations
        timeSteps = diff(timeStamps);
        figure(1);
        plot(timeSteps, 'r'); % Actual frame durations
        yline(1/frameRate, 'b'); % Target frame rate
        yline(2/frameRate, 'g'); % Double target frame rate
        yline(0,'g'); %0 time
        axis([0 maxFrames 0 2.5/frameRate]);
        title('Frame Rate');
        xlabel('Frame');
        ylabel('Duration (s)');
        legend('Measured Frame Rate', 'Target Frame Rate', 'Skipped Frame');
        
        % Plot deviations from target frame rate
        figure(2);
        deviation = timeSteps - (1/frameRate);
        plot(deviation, 'r');
        axis([0 maxFrames -2/frameRate 2/frameRate]);
        title('Deviations from Frame Rate');
        xlabel('Frame');
        ylabel('Difference Between Measured and Target Duration (s)');
    end
    
    % Reformat adjustment history array
    adjustmentArray = adjustmentArray(adjustmentArray ~= 0);
    col = length(adjustmentArray)/3;
    adjustmentArray = reshape(adjustmentArray, [3 col]);
    for i = 1:col
        adjustmentArray(:,i) = rgbToContrast(cal, adjustmentArray(:,i)', 'RGB', true);
    end
    adjustmentArray = adjustmentArray(2,:);
    
    % Display results and save data
    fprintf('chosen m cone contrast is %g \n', adjustmentArray(end));
    fprintf('adjustment history: ');
    fprintf('%g, ', adjustmentArray);
    fprintf('\n');
    save(fileLoc); 
    
catch e % Handle errors
    ListenChar(0);
    mglDisplayCursor(1);
    rethrow(e);
end
end
